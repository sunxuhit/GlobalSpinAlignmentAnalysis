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
#include "Utility/include/StSpinAlignmentFunctions.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"
// #include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

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

  mQCouEastRaw      = 0; // TPC EP East
  mQCouEastReCenter = 0;
  v_mQ2EastRaw.Set(0.0,0.0);
  v_mQ2EastReCenter.Set(0.0,0.0);

  mQCouWestRaw      = 0; // TPC EP West
  mQCouWestReCenter = 0;
  v_mQ2WestRaw.Set(0.0,0.0);
  v_mQ2WestReCenter.Set(0.0,0.0);
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

void StTpcEpManager::addTrackEastRaw(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2EastRaw += wgt*calq2Vector(picoTrack);
  mQCouEastRaw++;
}

void StTpcEpManager::addTrackWestRaw(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2WestRaw += wgt*calq2Vector(picoTrack);
  mQCouWestRaw++;
}

void StTpcEpManager::addTrackEastReCenter(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2EastReCenter += wgt*(calq2Vector(picoTrack) - getReCenterParEast());
  mQCouEastReCenter++;
}

void StTpcEpManager::addTrackWestReCenter(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2WestReCenter += wgt*(calq2Vector(picoTrack) - getReCenterParWest());
  mQCouWestReCenter++;
}
//---------------------------------------------------------------------------------
// ReCenterPar Correction
void StTpcEpManager::initTpcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQReCenterXEastVz%d",iVz);
    p_mTpcQReCenterXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    proName = Form("p_mTpcQReCenterYEastVz%d",iVz);
    p_mTpcQReCenterYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);

    proName = Form("p_mTpcQReCenterXWestVz%d",iVz);
    p_mTpcQReCenterXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    proName = Form("p_mTpcQReCenterYWestVz%d",iVz);
    p_mTpcQReCenterYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
  }
}

void StTpcEpManager::fillTpcReCenterEast(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);
  TVector2 q2Vector = calq2Vector(picoTrack);
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();

  p_mTpcQReCenterXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mTpcQReCenterYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);
}

void StTpcEpManager::fillTpcReCenterWest(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);
  TVector2 q2Vector = calq2Vector(picoTrack);
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();

  p_mTpcQReCenterXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mTpcQReCenterYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);
}

void StTpcEpManager::writeTpcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mTpcQReCenterXEast[iVz]->Write();
    p_mTpcQReCenterYEast[iVz]->Write();
    p_mTpcQReCenterXWest[iVz]->Write();
    p_mTpcQReCenterYWest[iVz]->Write();
  }
}

void StTpcEpManager::readTpcReCenterCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_TpcReCenterPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mReCenterPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQReCenterXEastVz%d",iVz);
    p_mTpcQReCenterXEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mTpcQReCenterYEastVz%d",iVz);
    p_mTpcQReCenterYEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mTpcQReCenterXWestVz%d",iVz);
    p_mTpcQReCenterXWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mTpcQReCenterYWestVz%d",iVz);
    p_mTpcQReCenterYWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
  }
}

TVector2 StTpcEpManager::getReCenterParEast()
{
  const int binX   = p_mTpcQReCenterXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mTpcQReCenterXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQReCenterYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mTpcQReCenterYEast[mVzBin]->GetBinContent(binY);

  TVector2 q2Vector(0.0,0.0);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StTpcEpManager::getReCenterParWest()
{
  const int binX   = p_mTpcQReCenterXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mTpcQReCenterXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQReCenterYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mTpcQReCenterYWest[mVzBin]->GetBinContent(binY);

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
      p_mTpcQShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
      proName = Form("p_mTpcQShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);

      proName = Form("p_mTpcQShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
      proName = Form("p_mTpcQShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    }
  }
}

void StTpcEpManager::fillTpcShiftEast()
{
  TVector2 Q2Vector = getQVectorEastReCenter();
  const double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi2Cos = TMath::Cos(2.0*((double)iShift+1.0)*Psi2);
    const double Psi2Sin = TMath::Sin(2.0*((double)iShift+1.0)*Psi2);
    p_mTpcQShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Cos);
    p_mTpcQShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Sin);
  }
}

void StTpcEpManager::fillTpcShiftWest()
{
  TVector2 Q2Vector = getQVectorWestReCenter();
  const double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi2Cos = TMath::Cos(2.0*((double)iShift+1.0)*Psi2);
    const double Psi2Sin = TMath::Sin(2.0*((double)iShift+1.0)*Psi2);
    p_mTpcQShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Cos);
    p_mTpcQShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Sin);
  }
}

void StTpcEpManager::writeTpcShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mTpcQShiftCosEast[iVz][iShift]->Write();
      p_mTpcQShiftSinEast[iVz][iShift]->Write();
      p_mTpcQShiftCosWest[iVz][iShift]->Write();
      p_mTpcQShiftSinWest[iVz][iShift]->Write();
    }
  }
}

void StTpcEpManager::readTpcShiftCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftParameter/file_%s_TpcShiftPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(InPutFile_Shift.Data());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mTpcQShiftCos%dEastVz%d",iShift,iVz);
      p_mTpcQShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName);
      proName = Form("p_mTpcQShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName);

      proName = Form("p_mTpcQShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName);
      proName = Form("p_mTpcQShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName);
    }
  }
}

double StTpcEpManager::calShiftAngle2East()
{
  TVector2 Q2Vector = getQVectorEastReCenter();
  const double PsiReCenter = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2

  double deltaPsi = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mTpcQShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mTpcQShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mTpcQShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mTpcQShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += 0.5*(2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*PsiReCenter));
  }

  double PsiShiftRaw = PsiReCenter + deltaPsi;
  double PsiShift    = angleShift(PsiShiftRaw);

  return PsiShift;
}

double StTpcEpManager::calShiftAngle2West()
{
  TVector2 Q2Vector  = getQVectorWestReCenter();

  const double PsiReCenter = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  double deltaPsi = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mTpcQShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mTpcQShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mTpcQShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mTpcQShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += 0.5*(2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*PsiReCenter));
  }

  double PsiShiftRaw = PsiReCenter + deltaPsi;
  double PsiShift = angleShift(PsiShiftRaw);

  return PsiShift;
}

double StTpcEpManager::angleShift(double PsiRaw)
{
  double PsiCorr = PsiRaw;
  if(PsiRaw > 0.5*TMath::Pi())
  {
    PsiCorr = PsiRaw - TMath::Pi();
  }
  if(PsiRaw < -0.5*TMath::Pi())
  {
    PsiCorr = PsiRaw + TMath::Pi();
  }

  return PsiCorr;
}
//---------------------------------------------------------------------------------
// Sub Event Plane Resolution
void StTpcEpManager::readTpcResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_%s_TpcEpResolution.root",globCons::str_mBeamType[mType].c_str(),g  lobCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.Data());
  p_mTpcEpResolutionSub = (TProfile*)file_mResolution->Get("p_mRes2_Sub");

  for(int iCent = 0; iCent < 9; ++iCent)
  {
    double valResSub, errResSub;
    double valResRaw = p_mTpcEpResolutionSub->GetBinContent(iCent+1);
    double errResRaw = p_mTpcEpResolutionSub->GetBinError(iCent+1);
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

    mTpcResSubVal[iCent] = valResSub;
    mTpcResSubErr[iCent] = errResSub;
  }
}

double StTpcEpManager::getTpcResSubVal(int cent9)
{
  return mTpcResSubVal[cent9];
}

double StTpcEpManager::getTpcResSubErr(int cent9)
{
  return mTpcResSubErr[cent9];
}
//---------------------------------------------------------------------------------
// QVector
TVector2 StTpcEpManager::getQVectorEastRaw()
{
  return v_mQ2EastRaw;
}

TVector2 StTpcEpManager::getQVectorWestRaw()
{
  return v_mQ2WestRaw;
}

TVector2 StTpcEpManager::getQVectorEastReCenter()
{
  return v_mQ2EastReCenter;
}

TVector2 StTpcEpManager::getQVectorWestReCenter()
{
  return v_mQ2WestReCenter;
}

int StTpcEpManager::getNumTrackEastRaw()
{
  return mQCouEastRaw;
}

int StTpcEpManager::getNumTrackWestRaw()
{
  return mQCouWestRaw;
}

int StTpcEpManager::getNumTrackEastReCenter()
{
  return mQCouEastReCenter;
}

int StTpcEpManager::getNumTrackWestReCenter()
{
  return mQCouWestReCenter;
}
//---------------------------------------------------------------------------------
