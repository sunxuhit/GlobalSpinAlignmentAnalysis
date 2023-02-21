#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TH2F.h"

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
  clearTpcEpManager();
}

StTpcEpManager::~StTpcEpManager()
{
  /* */
}
//---------------------------------------------------------------------------------
void StTpcEpManager::clearTpcEpManager()
{
  mCent9 = -1;
  mRunIndex = -1;
  mVzBin = -1;

  mQCouRawEast = 0; // raw EP
  mQCouRawWest = 0;
  v_mQ2RawWest.Set(0.0,0.0); // 2nd raw EP
  v_mQ2RawEast.Set(0.0,0.0);
  v_mQ3RawEast.Set(0.0,0.0); // 3rd raw EP
  v_mQ3RawWest.Set(0.0,0.0);

  mQCouReCenterEast = 0; // recenter EP
  mQCouReCenterWest = 0;
  v_mQ2ReCenterEast.Set(0.0,0.0); // 2nd recenter EP
  v_mQ2ReCenterWest.Set(0.0,0.0);
  v_mQ3ReCenterEast.Set(0.0,0.0);
  v_mQ3ReCenterWest.Set(0.0,0.0); // 3rd recenter EP
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

TVector2 StTpcEpManager::calq3Vector(StPicoTrack *picoTrack)
{
  const double phi = picoTrack->pMom().Phi(); // -pi to pi
  TVector2 q3Vector(0.0,0.0);

  const double q3x = TMath::Cos(3.0*phi);
  const double q3y = TMath::Sin(3.0*phi);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
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
  v_mQ3RawEast += wgt*calq3Vector(picoTrack);
  mQCouRawEast++;
}

void StTpcEpManager::addTrackRawWest(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2RawWest += wgt*calq2Vector(picoTrack);
  v_mQ3RawWest += wgt*calq3Vector(picoTrack);
  mQCouRawWest++;
}

void StTpcEpManager::addTrackReCenterEast(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2ReCenterEast += wgt*(calq2Vector(picoTrack) - getq2VecCtrEast());
  v_mQ3ReCenterEast += wgt*(calq3Vector(picoTrack) - getq3VecCtrEast());
  mQCouReCenterEast++;
}

void StTpcEpManager::addTrackReCenterWest(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2ReCenterWest += wgt*(calq2Vector(picoTrack) - getq2VecCtrWest());
  v_mQ3ReCenterWest += wgt*(calq3Vector(picoTrack) - getq3VecCtrWest());
  mQCouReCenterWest++;
}
//---------------------------------------------------------------------------------
// ReCenterPar Correction
void StTpcEpManager::initTpcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQ2ReCenterXEastVz%d",iVz); // 2nd EP
    p_mTpcQ2ReCenterXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ2ReCenterYEastVz%d",iVz);
    p_mTpcQ2ReCenterYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ2ReCenterXWestVz%d",iVz);
    p_mTpcQ2ReCenterXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ2ReCenterYWestVz%d",iVz);
    p_mTpcQ2ReCenterYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ3ReCenterXEastVz%d",iVz); // 3rd EP
    p_mTpcQ3ReCenterXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ3ReCenterYEastVz%d",iVz);
    p_mTpcQ3ReCenterYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ3ReCenterXWestVz%d",iVz);
    p_mTpcQ3ReCenterXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ3ReCenterYWestVz%d",iVz);
    p_mTpcQ3ReCenterYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
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

  TVector2 q3Vector = calq3Vector(picoTrack);
  const double q3x  = q3Vector.X();
  const double q3y  = q3Vector.Y();
  p_mTpcQ3ReCenterXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3x,wgt);
  p_mTpcQ3ReCenterYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3y,wgt);
}

void StTpcEpManager::fillTpcReCenterWest(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);

  TVector2 q2Vector = calq2Vector(picoTrack);
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();
  p_mTpcQ2ReCenterXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mTpcQ2ReCenterYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);

  TVector2 q3Vector = calq3Vector(picoTrack);
  const double q3x  = q3Vector.X();
  const double q3y  = q3Vector.Y();
  p_mTpcQ3ReCenterXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3x,wgt);
  p_mTpcQ3ReCenterYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3y,wgt);
}

void StTpcEpManager::writeTpcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mTpcQ2ReCenterXEast[iVz]->Write();
    p_mTpcQ2ReCenterYEast[iVz]->Write();
    p_mTpcQ2ReCenterXWest[iVz]->Write();
    p_mTpcQ2ReCenterYWest[iVz]->Write();

    p_mTpcQ3ReCenterXEast[iVz]->Write();
    p_mTpcQ3ReCenterYEast[iVz]->Write();
    p_mTpcQ3ReCenterXWest[iVz]->Write();
    p_mTpcQ3ReCenterYWest[iVz]->Write();
  }
}

void StTpcEpManager::readTpcReCenter()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_TpcReCenterPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mReCenterPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQ2ReCenterXEastVz%d",iVz); // 2nd EP
    p_mTpcQ2ReCenterXEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCenterYEastVz%d",iVz);
    p_mTpcQ2ReCenterYEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mTpcQ2ReCenterXWestVz%d",iVz);
    p_mTpcQ2ReCenterXWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCenterYWestVz%d",iVz);
    p_mTpcQ2ReCenterYWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mTpcQ3ReCenterXEastVz%d",iVz); // 3rd EP
    p_mTpcQ3ReCenterXEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mTpcQ3ReCenterYEastVz%d",iVz);
    p_mTpcQ3ReCenterYEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mTpcQ3ReCenterXWestVz%d",iVz);
    p_mTpcQ3ReCenterXWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mTpcQ3ReCenterYWestVz%d",iVz);
    p_mTpcQ3ReCenterYWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
  }
}

TVector2 StTpcEpManager::getq2VecCtrEast()
{
  const int binX   = p_mTpcQ2ReCenterXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mTpcQ2ReCenterXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ2ReCenterYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mTpcQ2ReCenterYEast[mVzBin]->GetBinContent(binY);

  TVector2 q2Vector(0.0,0.0);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StTpcEpManager::getq2VecCtrWest()
{
  const int binX   = p_mTpcQ2ReCenterXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mTpcQ2ReCenterXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ2ReCenterYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mTpcQ2ReCenterYWest[mVzBin]->GetBinContent(binY);

  TVector2 q2Vector(0.0,0.0);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StTpcEpManager::getq3VecCtrEast()
{
  const int binX   = p_mTpcQ3ReCenterXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3x = p_mTpcQ3ReCenterXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ3ReCenterYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3y = p_mTpcQ3ReCenterYEast[mVzBin]->GetBinContent(binY);

  TVector2 q3Vector(0.0,0.0);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
}

TVector2 StTpcEpManager::getq3VecCtrWest()
{
  const int binX   = p_mTpcQ3ReCenterXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3x = p_mTpcQ3ReCenterXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ3ReCenterYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3y = p_mTpcQ3ReCenterYWest[mVzBin]->GetBinContent(binY);

  TVector2 q3Vector(0.0,0.0);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
}
//---------------------------------------------------------------------------------
// Shift Correction
void StTpcEpManager::initTpcShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mTpcQ2ShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mTpcQ2ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ2ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ2ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ2ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ3ShiftCos%dEastVz%d",iShift,iVz); // 3rd EP
      p_mTpcQ3ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ3ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ3ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ3ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
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

  TVector2 Q3Vector = getQ3VecReCenterEast();
  const double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi3Cos = TMath::Cos(3.0*((double)iShift+1.0)*Psi3);
    const double Psi3Sin = TMath::Sin(3.0*((double)iShift+1.0)*Psi3);
    p_mTpcQ3ShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Cos);
    p_mTpcQ3ShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Sin);
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

  TVector2 Q3Vector = getQ3VecReCenterWest();
  const double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi3Cos = TMath::Cos(3.0*((double)iShift+1.0)*Psi3);
    const double Psi3Sin = TMath::Sin(3.0*((double)iShift+1.0)*Psi3);
    p_mTpcQ3ShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Cos);
    p_mTpcQ3ShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Sin);
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

      p_mTpcQ3ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ3ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinWest[iVz][iShift]->Write();
    }
  }
}

void StTpcEpManager::readTpcShift()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftParameter/file_%s_TpcShiftPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mTpcQ2ShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mTpcQ2ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ2ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dEastVz%d",iShift,iVz); // 3rd EP
      p_mTpcQ3ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
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

    deltaPsi2 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*Psi2ReCenter)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*Psi2ReCenter));
  }

  double Psi2ShiftRaw = Psi2ReCenter + deltaPsi2/2.0;
  double Psi2Shift    = transPsi2(Psi2ShiftRaw);

  return Psi2Shift;
}

double StTpcEpManager::getPsi2ShiftWest()
{
  TVector2 Q2Vector = getQ2VecReCenterWest();

  const double Psi2ReCenter = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  double deltaPsi2 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mTpcQ2ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mTpcQ2ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mTpcQ2ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mTpcQ2ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi2 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*Psi2ReCenter)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*Psi2ReCenter));
  }

  double Psi2ShiftRaw = Psi2ReCenter + deltaPsi2/2.0;
  double Psi2Shift    = transPsi2(Psi2ShiftRaw);

  return Psi2Shift;
}

double StTpcEpManager::transPsi2(double Psi2)
{
  double Psi2Corr = Psi2;
  if(Psi2 >  TMath::Pi()/2.0) Psi2Corr = Psi2 - TMath::Pi();
  if(Psi2 < -TMath::Pi()/2.0) Psi2Corr = Psi2 + TMath::Pi();

  return Psi2Corr;
}

double StTpcEpManager::getPsi3ShiftEast()
{
  TVector2 Q3Vector = getQ3VecReCenterEast();
  const double Psi3ReCenter = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3

  double deltaPsi3 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mTpcQ3ShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mTpcQ3ShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mTpcQ3ShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mTpcQ3ShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi3 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(3.0*((double)iShift+1.0)*Psi3ReCenter)+meanCos*TMath::Sin(3.0*((double)iShift+1.0)*Psi3ReCenter));
  }

  double Psi3ShiftRaw = Psi3ReCenter + deltaPsi3/3.0;
  double Psi3Shift    = transPsi3(Psi3ShiftRaw);

  return Psi3Shift;
}

double StTpcEpManager::getPsi3ShiftWest()
{
  TVector2 Q3Vector = getQ3VecReCenterWest();

  const double Psi3ReCenter = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0;
  double deltaPsi3 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mTpcQ3ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mTpcQ3ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mTpcQ3ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mTpcQ3ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi3 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(3.0*((double)iShift+1.0)*Psi3ReCenter)+meanCos*TMath::Sin(3.0*((double)iShift+1.0)*Psi3ReCenter));
  }

  double Psi3ShiftRaw = Psi3ReCenter + deltaPsi3/3.0;
  double Psi3Shift    = transPsi3(Psi3ShiftRaw);

  return Psi3Shift;
}

double StTpcEpManager::transPsi3(double Psi3)
{
  double Psi3Corr = Psi3;
  if(Psi3 >  TMath::Pi()/3.0) Psi3Corr = Psi3 - TMath::TwoPi()/3.0;
  if(Psi3 < -TMath::Pi()/3.0) Psi3Corr = Psi3 + TMath::TwoPi()/3.0;

  return Psi3Corr;
}
//---------------------------------------------------------------------------------
// Sub Event Plane Resolution
void StTpcEpManager::initTpcResolution()
{
  p_mTpcSubEp2Res = new TProfile("p_mTpcSubEp2Res","p_mTpcSubEp2Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  p_mTpcSubEp3Res = new TProfile("p_mTpcSubEp3Res","p_mTpcSubEp3Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
}

void StTpcEpManager::fillTpcResolution(double Psi2East, double Psi2West, double Psi3East, double Psi3West)
{
  double res2Sub = TMath::Cos(2.0*(Psi2East-Psi2West));
  p_mTpcSubEp2Res->Fill((double)mCent9,res2Sub);

  double res3Sub = TMath::Cos(3.0*(Psi3East-Psi3West));
  p_mTpcSubEp3Res->Fill((double)mCent9,res3Sub);
}

void StTpcEpManager::writeTpcResolution()
{
  p_mTpcSubEp2Res->Write();
  p_mTpcSubEp3Res->Write();
}

void StTpcEpManager::readTpcResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_%s_TpcEpResolution.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mTpcSubEp2Res = (TProfile*)file_mResolution->Get("p_mTpcSubEp2Res");
  p_mTpcSubEp3Res = (TProfile*)file_mResolution->Get("p_mTpcSubEp3Res");
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    mTpcSubEp2ResVal[iCent] = 0.0;
    mTpcSubEp2ResErr[iCent] = 0.0;
    mTpcSubEp3ResVal[iCent] = 0.0;
    mTpcSubEp3ResErr[iCent] = 0.0;
  }

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    double valRes2Sub = -999.9;
    double errRes2Sub = 1.0;
    double valRes3Sub = -999.9;
    double errRes3Sub = 1.0;
    double valRes2Raw = p_mTpcSubEp2Res->GetBinContent(iCent+1);
    double errRes2Raw = p_mTpcSubEp2Res->GetBinError(iCent+1);
    double valRes3Raw = p_mTpcSubEp3Res->GetBinContent(iCent+1);
    double errRes3Raw = p_mTpcSubEp3Res->GetBinError(iCent+1);
    if(valRes2Raw > 0)
    {
      valRes2Sub = TMath::Sqrt(valRes2Raw);
      errRes2Sub = errRes2Raw/(2.0*valRes2Sub);
    }

    if(valRes3Raw > 0)
    {
      valRes3Sub = TMath::Sqrt(valRes3Raw);
      errRes3Sub = errRes3Raw/(2.0*valRes3Sub);
    }

    mTpcSubEp2ResVal[iCent] = valRes2Sub;
    mTpcSubEp2ResErr[iCent] = errRes2Sub;
    mTpcSubEp3ResVal[iCent] = valRes3Sub;
    mTpcSubEp3ResErr[iCent] = errRes3Sub;
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

double StTpcEpManager::getTpcSubEp3ResVal(int cent9)
{
  return mTpcSubEp3ResVal[cent9];
}

double StTpcEpManager::getTpcSubEp3ResErr(int cent9)
{
  return mTpcSubEp3ResErr[cent9];
}
//---------------------------------------------------------------------------------
// Charged Hadron Elliptic and Triangular FLow
void StTpcEpManager::initTpcSubEpFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    std::string proName = Form("p_mTpcSubEpEFlowCent%d",i_cent);
    p_mTpcSubEpEFlow[i_cent] = new TProfile(proName.c_str(),proName.c_str(),50,0.0,10.0);
    proName = Form("p_mTpcSubEpTFlowCent%d",i_cent);
    p_mTpcSubEpTFlow[i_cent] = new TProfile(proName.c_str(),proName.c_str(),50,0.0,10.0);
  }
}

void StTpcEpManager::fillTpcSubEpEFlow(double pt, double v2, double reweight)
{
    p_mTpcSubEpEFlow[mCent9]->Fill(pt, v2, reweight);
}

void StTpcEpManager::fillTpcSubEpTFlow(double pt, double v3, double reweight)
{
    p_mTpcSubEpTFlow[mCent9]->Fill(pt, v3, reweight);
}

void StTpcEpManager::writeTpcSubEpFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    p_mTpcSubEpEFlow[i_cent]->Write();
    p_mTpcSubEpTFlow[i_cent]->Write();
  }
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
  double Psi2Raw    = transPsi2(Psi2);

  return Psi2Raw;
}

double StTpcEpManager::getPsi2RawWest()
{
  TVector2 Q2Vector = getQ2VecRawWest();
  double Psi2       = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2Raw    = transPsi2(Psi2);

  return Psi2Raw;
}

double StTpcEpManager::getPsi2ReCenterEast()
{
  TVector2 Q2Vector   = getQ2VecReCenterEast();
  double Psi2         = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2ReCenter = transPsi2(Psi2);

  return Psi2ReCenter;
}

double StTpcEpManager::getPsi2ReCenterWest()
{
  TVector2 Q2Vector   = getQ2VecReCenterWest();
  double Psi2         = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2ReCenter = transPsi2(Psi2);

  return Psi2ReCenter;
}

TVector2 StTpcEpManager::getQ3VecRawEast()
{
  return v_mQ3RawEast;
}

TVector2 StTpcEpManager::getQ3VecRawWest()
{
  return v_mQ3RawWest;
}

TVector2 StTpcEpManager::getQ3VecReCenterEast()
{
  return v_mQ3ReCenterEast;
}

TVector2 StTpcEpManager::getQ3VecReCenterWest()
{
  return v_mQ3ReCenterWest;
}

double StTpcEpManager::getPsi3RawEast()
{
  TVector2 Q3Vector = getQ3VecRawEast();
  double Psi3       = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
  double Psi3Raw    = transPsi3(Psi3);

  return Psi3Raw;
}

double StTpcEpManager::getPsi3RawWest()
{
  TVector2 Q3Vector = getQ3VecRawWest();
  double Psi3       = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
  double Psi3Raw    = transPsi3(Psi3);

  return Psi3Raw;
}

double StTpcEpManager::getPsi3ReCenterEast()
{
  TVector2 Q3Vector   = getQ3VecReCenterEast();
  double Psi3         = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
  double Psi3ReCenter = transPsi3(Psi3);

  return Psi3ReCenter;
}

double StTpcEpManager::getPsi3ReCenterWest()
{
  TVector2 Q3Vector   = getQ3VecReCenterWest();
  double Psi3         = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
  double Psi3ReCenter = transPsi3(Psi3);

  return Psi3ReCenter;
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
    std::string histName = Form("h_mTpcEp2RawEastCent%d",iCent); // 2nd EP
    h_mTpcEp2RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2RawWestCent%d",iCent);
    h_mTpcEp2RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2RawCorrCent%d",iCent);
    h_mTpcEp2RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());

    histName = Form("h_mTpcEp3RawEastCent%d",iCent); // 3rd EP
    h_mTpcEp3RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp3RawWestCent%d",iCent);
    h_mTpcEp3RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp3RawCorrCent%d",iCent);
    h_mTpcEp3RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StTpcEpManager::fillTpcSubEpRaw(double Psi2East, double Psi2West, double Psi3East, double Psi3West)
{
  h_mTpcEp2RawEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mTpcEp2RawWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mTpcEp2RawCorr[mCent9]->Fill(Psi2East,Psi2West);

  h_mTpcEp3RawEast[mCent9]->Fill(mRunIndex,Psi3East);
  h_mTpcEp3RawWest[mCent9]->Fill(mRunIndex,Psi3West);
  h_mTpcEp3RawCorr[mCent9]->Fill(Psi3East,Psi3West);
}

void StTpcEpManager::writeTpcSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2RawEast[iCent]->Write();
    h_mTpcEp2RawWest[iCent]->Write();
    h_mTpcEp2RawCorr[iCent]->Write();

    h_mTpcEp3RawEast[iCent]->Write();
    h_mTpcEp3RawWest[iCent]->Write();
    h_mTpcEp3RawCorr[iCent]->Write();
  }
}

// recenter EP
void StTpcEpManager::initTpcSubEpReCenter()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp2ReCenterEastCent%d",iCent); // 2nd EP
    h_mTpcEp2ReCenterEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2ReCenterWestCent%d",iCent);
    h_mTpcEp2ReCenterWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2ReCenterCorrCent%d",iCent);
    h_mTpcEp2ReCenterCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());

    histName = Form("h_mTpcEp3ReCenterEastCent%d",iCent); // 3rd EP
    h_mTpcEp3ReCenterEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp3ReCenterWestCent%d",iCent);
    h_mTpcEp3ReCenterWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp3ReCenterCorrCent%d",iCent);
    h_mTpcEp3ReCenterCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StTpcEpManager::fillTpcSubEpReCenter(double Psi2East, double Psi2West, double Psi3East, double Psi3West)
{
  h_mTpcEp2ReCenterEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mTpcEp2ReCenterWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mTpcEp2ReCenterCorr[mCent9]->Fill(Psi2East,Psi2West);

  h_mTpcEp3ReCenterEast[mCent9]->Fill(mRunIndex,Psi3East);
  h_mTpcEp3ReCenterWest[mCent9]->Fill(mRunIndex,Psi3West);
  h_mTpcEp3ReCenterCorr[mCent9]->Fill(Psi3East,Psi3West);
}

void StTpcEpManager::writeTpcSubEpReCenter()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2ReCenterEast[iCent]->Write();
    h_mTpcEp2ReCenterWest[iCent]->Write();
    h_mTpcEp2ReCenterCorr[iCent]->Write();

    h_mTpcEp3ReCenterEast[iCent]->Write();
    h_mTpcEp3ReCenterWest[iCent]->Write();
    h_mTpcEp3ReCenterCorr[iCent]->Write();
  }
}

// shift EP
void StTpcEpManager::initTpcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp2ShiftEastCent%d",iCent); // 2nd EP
    h_mTpcEp2ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2ShiftWestCent%d",iCent);
    h_mTpcEp2ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp2ShiftCorrCent%d",iCent);
    h_mTpcEp2ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());

    histName = Form("h_mTpcEp3ShiftEastCent%d",iCent); // 3rd EP
    h_mTpcEp3ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp3ShiftWestCent%d",iCent);
    h_mTpcEp3ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mTpcEp3ShiftCorrCent%d",iCent);
    h_mTpcEp3ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StTpcEpManager::fillTpcSubEpShift(double Psi2East, double Psi2West, double Psi3East, double Psi3West)
{
  h_mTpcEp2ShiftEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mTpcEp2ShiftWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mTpcEp2ShiftCorr[mCent9]->Fill(Psi2East,Psi2West);

  h_mTpcEp3ShiftEast[mCent9]->Fill(mRunIndex,Psi3East);
  h_mTpcEp3ShiftWest[mCent9]->Fill(mRunIndex,Psi3West);
  h_mTpcEp3ShiftCorr[mCent9]->Fill(Psi3East,Psi3West);
}

void StTpcEpManager::writeTpcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2ShiftEast[iCent]->Write();
    h_mTpcEp2ShiftWest[iCent]->Write();
    h_mTpcEp2ShiftCorr[iCent]->Write();

    h_mTpcEp3ShiftEast[iCent]->Write();
    h_mTpcEp3ShiftWest[iCent]->Write();
    h_mTpcEp3ShiftCorr[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
