#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"

#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"

ClassImp(StTpcEpManager)

//---------------------------------------------------------------------------------
StTpcEpManager::StTpcEpManager(int beamType) : mType(beamType)
{
  // mEnergy = energy;
  clearTpcEp();
}

StTpcEpManager::~StTpcEpManager()
{
  /* */
}
//---------------------------------------------------------------------------------
void StTpcEpManager::clearTpcEp()
{
  mCent9 = -1;
  mRunIndex = -1;
  mVzBin = -1;

  mQCouRawEast      = 0; // TPC EP East
  mQCouReCenterEast = 0;
  v_mQ2RawEast.Set(0.0,0.0);
  v_mQ2ReCenterEast.Set(0.0,0.0);

  mQCouRawWest      = 0; // TPC EP West
  mQCouReCenterWest = 0;
  v_mQ2RawWest.Set(0.0,0.0);
  v_mQ2ReCenterWest.Set(0.0,0.0);
}

void StTpcEpManager::initTpcEpManager(int cent9, int runIndex, int vzBin)
{
  mCent9    = cent9;
  mRunIndex = runIndex;
  mVzBin    = vzBin;
}
//---------------------------------------------------------------------------------
// Utilities
TVector2 StTpcEpManager::calq2Vector(StPicoTrack *picoTrack)
{
  const double phi = picoTrack->pMom().Phi(); // -pi to pi
  TVector2 q2Vector(0.0,0.0);

  const double q2x = TMath::Cos(2.0*phi);
  const double q2y = TMath::Sin(2.0*phi);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

double StTpcEpManager::getWeight(StPicoTrack *picoTrack)
{
  const double pt = picoTrack->pMom().Perp();
  double wgt;
  if(pt <= anaUtils::mPrimPtEpWeight[mType])
  {
    wgt = pt;
  }
  if(pt > anaUtils::mPrimPtEpWeight[mType])
  {
    wgt = anaUtils::mPrimPtEpWeight[mType];
  }

  return wgt;
}

void StTpcEpManager::addTrackRawEast(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2RawEast += wgt*calq2Vector(picoTrack);
  mQCouRawEast++;
}

void StTpcEpManager::addTrackRawWest(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2RawWest += wgt*calq2Vector(picoTrack);
  mQCouRawWest++;
}

void StTpcEpManager::addTrackReCenterEast(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2ReCenterEast += wgt*(calq2Vector(picoTrack) - getReCenterParEast());
  mQCouReCenterEast++;
}

void StTpcEpManager::addTrackReCenterWest(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2ReCenterWest += wgt*(calq2Vector(picoTrack) - getReCenterParWest());
  mQCouReCenterWest++;
}
//---------------------------------------------------------------------------------
// ReCenterPar Correction
void StTpcEpManager::initTpcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQ2ReCenterXEastVz%d",iVz);
    p_mTpcQ2ReCenterXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ2ReCenterYEastVz%d",iVz);
    p_mTpcQ2ReCenterYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ2ReCenterXWestVz%d",iVz);
    p_mTpcQ2ReCenterXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ2ReCenterYWestVz%d",iVz);
    p_mTpcQ2ReCenterYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  }
}

void StTpcEpManager::fillTpcReCenterEast(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);
  TVector2 q2Vector = calq2Vector(picoTrack);
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();

  p_mTpcQ2ReCenterXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mTpcQ2ReCenterYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);
}

void StTpcEpManager::fillTpcReCenterWest(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);
  TVector2 q2Vector = calq2Vector(picoTrack);
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();

  p_mTpcQ2ReCenterXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mTpcQ2ReCenterYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);
}

void StTpcEpManager::writeTpcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mTpcQ2ReCenterXEast[iVz]->Write();
    p_mTpcQ2ReCenterYEast[iVz]->Write();
    p_mTpcQ2ReCenterXWest[iVz]->Write();
    p_mTpcQ2ReCenterYWest[iVz]->Write();
  }
}

void StTpcEpManager::readTpcReCenterCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_TpcReCenterPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mReCenterPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQ2ReCenterXEastVz%d",iVz);
    p_mTpcQ2ReCenterXEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCenterYEastVz%d",iVz);
    p_mTpcQ2ReCenterYEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mTpcQ2ReCenterXWestVz%d",iVz);
    p_mTpcQ2ReCenterXWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCenterYWestVz%d",iVz);
    p_mTpcQ2ReCenterYWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
  }
}

TVector2 StTpcEpManager::getReCenterParEast()
{
  const int binX   = p_mTpcQ2ReCenterXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mTpcQ2ReCenterXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ2ReCenterYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mTpcQ2ReCenterYEast[mVzBin]->GetBinContent(binY);

  TVector2 q2Vector(0.0,0.0);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StTpcEpManager::getReCenterParWest()
{
  const int binX   = p_mTpcQ2ReCenterXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mTpcQ2ReCenterXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ2ReCenterYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mTpcQ2ReCenterYWest[mVzBin]->GetBinContent(binY);

  TVector2 q2Vector(0.0,0.0);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}
//---------------------------------------------------------------------------------
// Shift Correction
void StTpcEpManager::initTpcShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mTpcQShiftCos%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    }
  }
}

void StTpcEpManager::fillTpcShiftEast()
{
  TVector2 Q2Vector = getQ2VecReCenterEast();
  const double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi2Cos = TMath::Cos(2.0*((double)iShift+1.0)*Psi2);
    const double Psi2Sin = TMath::Sin(2.0*((double)iShift+1.0)*Psi2);
    p_mTpcQ2ShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Cos);
    p_mTpcQ2ShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Sin);
  }
}

void StTpcEpManager::fillTpcShiftWest()
{
  TVector2 Q2Vector = getQ2VecReCenterWest();
  const double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi2Cos = TMath::Cos(2.0*((double)iShift+1.0)*Psi2);
    const double Psi2Sin = TMath::Sin(2.0*((double)iShift+1.0)*Psi2);
    p_mTpcQ2ShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Cos);
    p_mTpcQ2ShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Sin);
  }
}

void StTpcEpManager::writeTpcShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mTpcQ2ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ2ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinWest[iVz][iShift]->Write();
    }
  }
}

void StTpcEpManager::readTpcShiftCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftParameter/file_%s_TpcShiftPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mTpcQShiftCos%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

double StTpcEpManager::getPsi2ShiftEast()
{
  TVector2 Q2Vector = getQ2VecReCenterEast();
  const double Psi2ReCenter = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2

  double deltaPsi2 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mTpcQ2ShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mTpcQ2ShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mTpcQ2ShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mTpcQ2ShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi2 += 0.5*(2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*Psi2ReCenter)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*Psi2ReCenter));
  }

  double Psi2ShiftRaw = Psi2ReCenter + deltaPsi2;
  double Psi2Shift    = angleShift(Psi2ShiftRaw);

  return Psi2Shift;
}

double StTpcEpManager::getPsi2ShiftWest()
{
  TVector2 Q2Vector  = getQ2VecReCenterWest();

  const double Psi2ReCenter = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  double deltaPsi2 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mTpcQ2ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mTpcQ2ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mTpcQ2ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mTpcQ2ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi2 += 0.5*(2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*Psi2ReCenter)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*Psi2ReCenter));
  }

  double Psi2ShiftRaw = Psi2ReCenter + deltaPsi2;
  double Psi2Shift    = angleShift(Psi2ShiftRaw);

  return Psi2Shift;
}

double StTpcEpManager::angleShift(double Psi2Raw)
{
  double Psi2Corr = Psi2Raw;
  if(Psi2Raw > 0.5*TMath::Pi())
  {
    Psi2Corr = Psi2Raw - TMath::Pi();
  }
  if(Psi2Raw < -0.5*TMath::Pi())
  {
    Psi2Corr = Psi2Raw + TMath::Pi();
  }

  return Psi2Corr;
}
//---------------------------------------------------------------------------------
// Sub Event Plane Resolution
void StTpcEpManager::initTpcResolution()
{
  p_mTpcSubEp2Res = new TProfile("p_mTpcSubEp2Res","p_mTpcSubEp2Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
}

void StTpcEpManager::fillTpcResolution(double Psi2East, double Psi2West)
{
  double res2Sub = TMath::Cos(2.0*(Psi2East-Psi2West));
  p_mTpcSubEp2Res->Fill((double)mCent9,res2Sub);
}

void StTpcEpManager::writeTpcResolution()
{
  p_mTpcSubEp2Res->Write();
}

void StTpcEpManager::readTpcResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_%s_TpcEpResolution.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mTpcSubEp2Res = (TProfile*)file_mResolution->Get("p_mTpcSubEp2Res");

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    mTpcSubEp2ResVal[iCent] = 0.0;
    mTpcSubEp2ResErr[iCent] = 0.0;
  }

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    double valResSub, errResSub;
    double valResRaw = p_mTpcSubEp2Res->GetBinContent(iCent+1);
    double errResRaw = p_mTpcSubEp2Res->GetBinError(iCent+1);
    if(valResRaw <= 0)
    {
      valResSub = -999.9;
      errResSub = 1.0;
    }
    else
    {
      valResSub = TMath::Sqrt(valResRaw);
      errResSub = errResRaw/(2*valResSub);
    }

    mTpcSubEp2ResVal[iCent] = valResSub;
    mTpcSubEp2ResErr[iCent] = errResSub;
  }
  file_mResolution->Close();
}

double StTpcEpManager::getTpcSubEp2ResVal(int cent9)
{
  return mTpcSubEp2ResVal[cent9];
}

double StTpcEpManager::getTpcSubEp2ResErr(int cent9)
{
  return mTpcSubEp2ResErr[cent9];
}
//---------------------------------------------------------------------------------
// QVector
TVector2 StTpcEpManager::getQ2VecRawEast()
{
  return v_mQ2RawEast;
}

TVector2 StTpcEpManager::getQ2VecRawWest()
{
  return v_mQ2RawWest;
}

TVector2 StTpcEpManager::getQ2VecReCenterEast()
{
  return v_mQ2ReCenterEast;
}

TVector2 StTpcEpManager::getQ2VecReCenterWest()
{
  return v_mQ2ReCenterWest;
}

double StTpcEpManager::getPsi2RawEast()
{
  TVector2 Q2Vector = getQ2VecRawEast();
  double Psi2       = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2Raw    = angleShift(Psi2);

  return Psi2Raw;
}

double StTpcEpManager::getPsi2RawWest()
{
  TVector2 Q2Vector = getQ2VecRawWest();
  double Psi2       = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2Raw    = angleShift(Psi2);

  return Psi2Raw;
}

double StTpcEpManager::getPsi2ReCenterEast()
{
  TVector2 Q2Vector   = getQ2VecReCenterEast();
  double Psi2         = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2ReCenter = angleShift(Psi2);

  return Psi2ReCenter;
}

double StTpcEpManager::getPsi2ReCenterWest()
{
  TVector2 Q2Vector   = getQ2VecReCenterWest();
  double Psi2         = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2ReCenter = angleShift(Psi2);

  return Psi2ReCenter;
}

int StTpcEpManager::getNumTrkRawEast()
{
  return mQCouRawEast;
}

int StTpcEpManager::getNumTrkRawWest()
{
  return mQCouRawWest;
}

int StTpcEpManager::getNumTrkReCenterEast()
{
  return mQCouReCenterEast;
}

int StTpcEpManager::getNumTrkReCenterWest()
{
  return mQCouReCenterWest;
}
//---------------------------------------------------------------------------------
// raw EP
void StTpcEpManager::initTpcSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp2RawEastCent%d",iCent);
    h_mTpcEp2RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2RawWestCent%d",iCent);
    h_mTpcEp2RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2RawCorrCent%d",iCent);
    h_mTpcEp2RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StTpcEpManager::fillTpcSubEpRaw(double Psi2East, double Psi2West)
{
  h_mTpcEp2RawEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mTpcEp2RawWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mTpcEp2RawCorr[mCent9]->Fill(Psi2East,Psi2West);
}

void StTpcEpManager::writeTpcSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2RawEast[iCent]->Write();
    h_mTpcEp2RawWest[iCent]->Write();
    h_mTpcEp2RawCorr[iCent]->Write();
  }
}

// recenter EP
void StTpcEpManager::initTpcSubEpReCenter()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp2ReCenterEastCent%d",iCent);
    h_mTpcEp2ReCenterEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2ReCenterWestCent%d",iCent);
    h_mTpcEp2ReCenterWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2ReCenterCorrCent%d",iCent);
    h_mTpcEp2ReCenterCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StTpcEpManager::fillTpcSubEpReCenter(double Psi2East, double Psi2West)
{
  h_mTpcEp2ReCenterEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mTpcEp2ReCenterWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mTpcEp2ReCenterCorr[mCent9]->Fill(Psi2East,Psi2West);
}

void StTpcEpManager::writeTpcSubEpReCenter()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2ReCenterEast[iCent]->Write();
    h_mTpcEp2ReCenterWest[iCent]->Write();
    h_mTpcEp2ReCenterCorr[iCent]->Write();
  }
}

// shift EP
void StTpcEpManager::initTpcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp2ShiftEastCent%d",iCent);
    h_mTpcEp2ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2ShiftWestCent%d",iCent);
    h_mTpcEp2ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2ShiftCorrCent%d",iCent);
    h_mTpcEp2ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StTpcEpManager::fillTpcSubEpShift(double Psi2East, double Psi2West)
{
  h_mTpcEp2ShiftEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mTpcEp2ShiftWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mTpcEp2ShiftCorr[mCent9]->Fill(Psi2East,Psi2West);
}

void StTpcEpManager::writeTpcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2ShiftEast[iCent]->Write();
    h_mTpcEp2ShiftWest[iCent]->Write();
    h_mTpcEp2ShiftCorr[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
