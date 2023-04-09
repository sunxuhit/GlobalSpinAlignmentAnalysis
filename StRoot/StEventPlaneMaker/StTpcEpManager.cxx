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

#include <iostream>

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
  mQCouRawFull = 0;
  v_mQ1RawWest.Set(0.0,0.0); // 1st raw EP
  v_mQ1RawEast.Set(0.0,0.0);
  v_mQ1RawFull.Set(0.0,0.0);
  v_mQ2RawWest.Set(0.0,0.0); // 2nd raw EP
  v_mQ2RawEast.Set(0.0,0.0);
  v_mQ2RawFull.Set(0.0,0.0);
  v_mQ3RawEast.Set(0.0,0.0); // 3rd raw EP
  v_mQ3RawWest.Set(0.0,0.0);
  v_mQ3RawFull.Set(0.0,0.0);

  mQCouReCtrEast = 0; // recenter EP
  mQCouReCtrWest = 0;
  mQCouReCtrFull = 0;
  v_mQ1ReCtrEast.Set(0.0,0.0); // 1st recenter EP
  v_mQ1ReCtrWest.Set(0.0,0.0);
  v_mQ1ReCtrFull.Set(0.0,0.0);
  v_mQ2ReCtrEast.Set(0.0,0.0); // 2nd recenter EP
  v_mQ2ReCtrWest.Set(0.0,0.0);
  v_mQ2ReCtrFull.Set(0.0,0.0);
  v_mQ3ReCtrEast.Set(0.0,0.0); // 3rd recenter EP
  v_mQ3ReCtrWest.Set(0.0,0.0);
  v_mQ3ReCtrFull.Set(0.0,0.0);
}

void StTpcEpManager::initTpcEpManager(int cent9, int runIndex, int vzBin)
{
  mCent9    = cent9;
  mRunIndex = runIndex;
  mVzBin    = vzBin;
}
//---------------------------------------------------------------------------------
// Utilities
TVector2 StTpcEpManager::calq1Vector(StPicoTrack *picoTrack)
{
  const double phi = picoTrack->pMom().Phi(); // -pi to pi
  TVector2 q1Vector(0.0,0.0);

  const double q1x = TMath::Cos(1.0*phi);
  const double q1y = TMath::Sin(1.0*phi);
  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

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
  double wgt = (pt > anaUtils::mPrimPtEpWeight[mType]) ? anaUtils::mPrimPtEpWeight[mType] : pt;
  if(pt < anaUtils::mPrimPtEpMin[mType]) wgt = 0.0;

  // double wgt = pt;
  // if(pt > anaUtils::mPrimPtEpWeight[mType])
  // {
  //   wgt = anaUtils::mPrimPtEpWeight[mType];
  // }
  // if(pt <= anaUtils::mPrimPtEpWeight[mType])
  // {
  //   wgt = pt;
  // }
  // if(pt < anaUtils::mPrimPtEpMin[mType])
  // {
  //   wgt = 0.0;
  // }

  return wgt;
}

void StTpcEpManager::addTrackRawEast(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ1RawEast += wgt*calq1Vector(picoTrack);
  v_mQ2RawEast += wgt*calq2Vector(picoTrack);
  v_mQ3RawEast += wgt*calq3Vector(picoTrack);
  mQCouRawEast++;
}

void StTpcEpManager::addTrackRawWest(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ1RawWest += wgt*calq1Vector(picoTrack);
  v_mQ2RawWest += wgt*calq2Vector(picoTrack);
  v_mQ3RawWest += wgt*calq3Vector(picoTrack);
  mQCouRawWest++;
}

void StTpcEpManager::addTrackRawFull(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ1RawFull += wgt*calq1Vector(picoTrack);
  v_mQ2RawFull += wgt*calq2Vector(picoTrack);
  v_mQ3RawFull += wgt*calq3Vector(picoTrack);
  mQCouRawFull++;
}

void StTpcEpManager::addTrackReCtrEast(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ1ReCtrEast += wgt*(calq1Vector(picoTrack) - getq1VecCtrEast());
  v_mQ2ReCtrEast += wgt*(calq2Vector(picoTrack) - getq2VecCtrEast());
  v_mQ3ReCtrEast += wgt*(calq3Vector(picoTrack) - getq3VecCtrEast());
  mQCouReCtrEast++;
}

void StTpcEpManager::addTrackReCtrWest(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ1ReCtrWest += wgt*(calq1Vector(picoTrack) - getq1VecCtrWest());
  v_mQ2ReCtrWest += wgt*(calq2Vector(picoTrack) - getq2VecCtrWest());
  v_mQ3ReCtrWest += wgt*(calq3Vector(picoTrack) - getq3VecCtrWest());
  mQCouReCtrWest++;
}

void StTpcEpManager::addTrackReCtrFull(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ1ReCtrFull += wgt*(calq1Vector(picoTrack) - getq1VecCtrFull());
  v_mQ2ReCtrFull += wgt*(calq2Vector(picoTrack) - getq2VecCtrFull());
  v_mQ3ReCtrFull += wgt*(calq3Vector(picoTrack) - getq3VecCtrFull());
  mQCouReCtrFull++;
}
//---------------------------------------------------------------------------------
// ReCenterPar Correction
void StTpcEpManager::initTpcReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQ1ReCtrXEastVz%d",iVz); // 1st EP
    p_mTpcQ1ReCtrXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ1ReCtrYEastVz%d",iVz);
    p_mTpcQ1ReCtrYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ1ReCtrXWestVz%d",iVz);
    p_mTpcQ1ReCtrXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ1ReCtrYWestVz%d",iVz);
    p_mTpcQ1ReCtrYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ1ReCtrXFullVz%d",iVz);
    p_mTpcQ1ReCtrXFull[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ1ReCtrYFullVz%d",iVz);
    p_mTpcQ1ReCtrYFull[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ2ReCtrXEastVz%d",iVz); // 2nd EP
    p_mTpcQ2ReCtrXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ2ReCtrYEastVz%d",iVz);
    p_mTpcQ2ReCtrYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ2ReCtrXWestVz%d",iVz);
    p_mTpcQ2ReCtrXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ2ReCtrYWestVz%d",iVz);
    p_mTpcQ2ReCtrYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ2ReCtrXFullVz%d",iVz);
    p_mTpcQ2ReCtrXFull[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ2ReCtrYFullVz%d",iVz);
    p_mTpcQ2ReCtrYFull[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ3ReCtrXEastVz%d",iVz); // 3rd EP
    p_mTpcQ3ReCtrXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ3ReCtrYEastVz%d",iVz);
    p_mTpcQ3ReCtrYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ3ReCtrXWestVz%d",iVz);
    p_mTpcQ3ReCtrXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ3ReCtrYWestVz%d",iVz);
    p_mTpcQ3ReCtrYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mTpcQ3ReCtrXFullVz%d",iVz);
    p_mTpcQ3ReCtrXFull[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mTpcQ3ReCtrYFullVz%d",iVz);
    p_mTpcQ3ReCtrYFull[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  }
}

void StTpcEpManager::fillTpcReCtrEast(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);

  TVector2 q1Vector = calq1Vector(picoTrack); // 1st EP
  const double q1x  = q1Vector.X();
  const double q1y  = q1Vector.Y();
  p_mTpcQ1ReCtrXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1x,wgt);
  p_mTpcQ1ReCtrYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1y,wgt);

  TVector2 q2Vector = calq2Vector(picoTrack); // 2nd EP
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();
  p_mTpcQ2ReCtrXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mTpcQ2ReCtrYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);

  TVector2 q3Vector = calq3Vector(picoTrack); // 3rd EP
  const double q3x  = q3Vector.X();
  const double q3y  = q3Vector.Y();
  p_mTpcQ3ReCtrXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3x,wgt);
  p_mTpcQ3ReCtrYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3y,wgt);
}

void StTpcEpManager::fillTpcReCtrWest(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);

  TVector2 q1Vector = calq1Vector(picoTrack); // 1st EP
  const double q1x  = q1Vector.X();
  const double q1y  = q1Vector.Y();
  p_mTpcQ1ReCtrXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1x,wgt);
  p_mTpcQ1ReCtrYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1y,wgt);

  TVector2 q2Vector = calq2Vector(picoTrack); // 2nd EP
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();
  p_mTpcQ2ReCtrXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mTpcQ2ReCtrYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);

  TVector2 q3Vector = calq3Vector(picoTrack); // 3rd EP
  const double q3x  = q3Vector.X();
  const double q3y  = q3Vector.Y();
  p_mTpcQ3ReCtrXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3x,wgt);
  p_mTpcQ3ReCtrYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3y,wgt);
}

void StTpcEpManager::fillTpcReCtrFull(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);

  TVector2 q1Vector = calq1Vector(picoTrack); // 1st EP
  const double q1x  = q1Vector.X();
  const double q1y  = q1Vector.Y();
  p_mTpcQ1ReCtrXFull[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1x,wgt);
  p_mTpcQ1ReCtrYFull[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1y,wgt);

  TVector2 q2Vector = calq2Vector(picoTrack); // 2nd EP
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();
  p_mTpcQ2ReCtrXFull[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mTpcQ2ReCtrYFull[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);

  TVector2 q3Vector = calq3Vector(picoTrack); // 3rd EP
  const double q3x  = q3Vector.X();
  const double q3y  = q3Vector.Y();
  p_mTpcQ3ReCtrXFull[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3x,wgt);
  p_mTpcQ3ReCtrYFull[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3y,wgt);
}

void StTpcEpManager::writeTpcReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mTpcQ1ReCtrXEast[iVz]->Write();
    p_mTpcQ1ReCtrYEast[iVz]->Write();
    p_mTpcQ1ReCtrXWest[iVz]->Write();
    p_mTpcQ1ReCtrYWest[iVz]->Write();
    p_mTpcQ1ReCtrXFull[iVz]->Write();
    p_mTpcQ1ReCtrYFull[iVz]->Write();

    p_mTpcQ2ReCtrXEast[iVz]->Write();
    p_mTpcQ2ReCtrYEast[iVz]->Write();
    p_mTpcQ2ReCtrXWest[iVz]->Write();
    p_mTpcQ2ReCtrYWest[iVz]->Write();
    p_mTpcQ2ReCtrXFull[iVz]->Write();
    p_mTpcQ2ReCtrYFull[iVz]->Write();

    p_mTpcQ3ReCtrXEast[iVz]->Write();
    p_mTpcQ3ReCtrYEast[iVz]->Write();
    p_mTpcQ3ReCtrXWest[iVz]->Write();
    p_mTpcQ3ReCtrYWest[iVz]->Write();
    p_mTpcQ3ReCtrXFull[iVz]->Write();
    p_mTpcQ3ReCtrYFull[iVz]->Write();
  }
}

void StTpcEpManager::readTpcReCtr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_TpcReCenterPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mReCtrPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQ1ReCtrXEastVz%d",iVz); // 1st EP
    p_mTpcQ1ReCtrXEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mTpcQ1ReCtrYEastVz%d",iVz);
    p_mTpcQ1ReCtrYEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mTpcQ1ReCtrXWestVz%d",iVz);
    p_mTpcQ1ReCtrXWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mTpcQ1ReCtrYWestVz%d",iVz);
    p_mTpcQ1ReCtrYWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mTpcQ1ReCtrXFullVz%d",iVz);
    p_mTpcQ1ReCtrXFull[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mTpcQ1ReCtrYFullVz%d",iVz);
    p_mTpcQ1ReCtrYFull[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mTpcQ2ReCtrXEastVz%d",iVz); // 2nd EP
    p_mTpcQ2ReCtrXEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCtrYEastVz%d",iVz);
    p_mTpcQ2ReCtrYEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mTpcQ2ReCtrXWestVz%d",iVz);
    p_mTpcQ2ReCtrXWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCtrYWestVz%d",iVz);
    p_mTpcQ2ReCtrYWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mTpcQ2ReCtrXFullVz%d",iVz);
    p_mTpcQ2ReCtrXFull[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCtrYFullVz%d",iVz);
    p_mTpcQ2ReCtrYFull[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mTpcQ3ReCtrXEastVz%d",iVz); // 3rd EP
    p_mTpcQ3ReCtrXEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mTpcQ3ReCtrYEastVz%d",iVz);
    p_mTpcQ3ReCtrYEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mTpcQ3ReCtrXWestVz%d",iVz);
    p_mTpcQ3ReCtrXWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mTpcQ3ReCtrYWestVz%d",iVz);
    p_mTpcQ3ReCtrYWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mTpcQ3ReCtrXFullVz%d",iVz);
    p_mTpcQ3ReCtrXFull[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mTpcQ3ReCtrYFullVz%d",iVz);
    p_mTpcQ3ReCtrYFull[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
  }
}

TVector2 StTpcEpManager::getq1VecCtrEast()
{
  TVector2 q1Vector(0.0,0.0);

  const int binX   = p_mTpcQ1ReCtrXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1x = p_mTpcQ1ReCtrXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ1ReCtrYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1y = p_mTpcQ1ReCtrYEast[mVzBin]->GetBinContent(binY);

  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

TVector2 StTpcEpManager::getq1VecCtrWest()
{
  TVector2 q1Vector(0.0,0.0);

  const int binX   = p_mTpcQ1ReCtrXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1x = p_mTpcQ1ReCtrXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ1ReCtrYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1y = p_mTpcQ1ReCtrYWest[mVzBin]->GetBinContent(binY);

  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

TVector2 StTpcEpManager::getq1VecCtrFull()
{
  TVector2 q1Vector(0.0,0.0);

  const int binX   = p_mTpcQ1ReCtrXFull[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1x = p_mTpcQ1ReCtrXFull[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ1ReCtrYFull[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1y = p_mTpcQ1ReCtrYFull[mVzBin]->GetBinContent(binY);

  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

TVector2 StTpcEpManager::getq2VecCtrEast()
{
  TVector2 q2Vector(0.0,0.0);

  const int binX   = p_mTpcQ2ReCtrXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mTpcQ2ReCtrXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ2ReCtrYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mTpcQ2ReCtrYEast[mVzBin]->GetBinContent(binY);

  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StTpcEpManager::getq2VecCtrWest()
{
  TVector2 q2Vector(0.0,0.0);

  const int binX   = p_mTpcQ2ReCtrXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mTpcQ2ReCtrXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ2ReCtrYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mTpcQ2ReCtrYWest[mVzBin]->GetBinContent(binY);

  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StTpcEpManager::getq2VecCtrFull()
{
  TVector2 q2Vector(0.0,0.0);

  const int binX   = p_mTpcQ2ReCtrXFull[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mTpcQ2ReCtrXFull[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ2ReCtrYFull[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mTpcQ2ReCtrYFull[mVzBin]->GetBinContent(binY);

  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StTpcEpManager::getq3VecCtrEast()
{
  TVector2 q3Vector(0.0,0.0);

  const int binX   = p_mTpcQ3ReCtrXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3x = p_mTpcQ3ReCtrXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ3ReCtrYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3y = p_mTpcQ3ReCtrYEast[mVzBin]->GetBinContent(binY);

  q3Vector.Set(q3x,q3y);

  return q3Vector;
}

TVector2 StTpcEpManager::getq3VecCtrWest()
{
  TVector2 q3Vector(0.0,0.0);

  const int binX   = p_mTpcQ3ReCtrXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3x = p_mTpcQ3ReCtrXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ3ReCtrYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3y = p_mTpcQ3ReCtrYWest[mVzBin]->GetBinContent(binY);

  q3Vector.Set(q3x,q3y);

  return q3Vector;
}

TVector2 StTpcEpManager::getq3VecCtrFull()
{
  TVector2 q3Vector(0.0,0.0);

  const int binX   = p_mTpcQ3ReCtrXFull[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3x = p_mTpcQ3ReCtrXFull[mVzBin]->GetBinContent(binX);

  const int binY   = p_mTpcQ3ReCtrYFull[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3y = p_mTpcQ3ReCtrYFull[mVzBin]->GetBinContent(binY);

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
      std::string proName = Form("p_mTpcQ1ShiftCos%dEastVz%d",iShift,iVz); // 1st EP
      p_mTpcQ1ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ1ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ1ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ1ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ1ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ1ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ1ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ1ShiftCos%dFullVz%d",iShift,iVz);
      p_mTpcQ1ShiftCosFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ1ShiftSin%dFullVz%d",iShift,iVz);
      p_mTpcQ1ShiftSinFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ2ShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mTpcQ2ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ2ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ2ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ2ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ2ShiftCos%dFullVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ2ShiftSin%dFullVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ3ShiftCos%dEastVz%d",iShift,iVz); // 3rd EP
      p_mTpcQ3ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ3ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ3ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ3ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mTpcQ3ShiftCos%dFullVz%d",iShift,iVz);
      p_mTpcQ3ShiftCosFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mTpcQ3ShiftSin%dFullVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    }
  }
}

void StTpcEpManager::fillTpcShiftEast()
{
  TVector2 Q1Vector = getQ1VecReCtrEast(); // 1st EP
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi1Cos = TMath::Cos(1.0*((double)iShift+1.0)*Psi1);
      const double Psi1Sin = TMath::Sin(1.0*((double)iShift+1.0)*Psi1);
      p_mTpcQ1ShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
      p_mTpcQ1ShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
    }
  }

  TVector2 Q2Vector = getQ2VecReCtrEast(); // 2nd EP
  if(Q2Vector.Mod() > 0.0)
  {
    const double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi2Cos = TMath::Cos(2.0*((double)iShift+1.0)*Psi2);
      const double Psi2Sin = TMath::Sin(2.0*((double)iShift+1.0)*Psi2);
      p_mTpcQ2ShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Cos);
      p_mTpcQ2ShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Sin);
    }
  }

  TVector2 Q3Vector = getQ3VecReCtrEast(); // 3rd EP
  if(Q3Vector.Mod() > 0.0)
  {
    const double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi3Cos = TMath::Cos(3.0*((double)iShift+1.0)*Psi3);
      const double Psi3Sin = TMath::Sin(3.0*((double)iShift+1.0)*Psi3);
      p_mTpcQ3ShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Cos);
      p_mTpcQ3ShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Sin);
    }
  }
}

void StTpcEpManager::fillTpcShiftWest()
{
  TVector2 Q1Vector = getQ1VecReCtrWest(); // 1st EP
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi1Cos = TMath::Cos(1.0*((double)iShift+1.0)*Psi1);
      const double Psi1Sin = TMath::Sin(1.0*((double)iShift+1.0)*Psi1);
      p_mTpcQ1ShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
      p_mTpcQ1ShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
    }
  }

  TVector2 Q2Vector = getQ2VecReCtrWest(); // 2nd EP
  if(Q2Vector.Mod() > 0.0)
  {
    const double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi2Cos = TMath::Cos(2.0*((double)iShift+1.0)*Psi2);
      const double Psi2Sin = TMath::Sin(2.0*((double)iShift+1.0)*Psi2);
      p_mTpcQ2ShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Cos);
      p_mTpcQ2ShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Sin);
    }
  }

  TVector2 Q3Vector = getQ3VecReCtrWest(); // 3rd EP
  if(Q3Vector.Mod() > 0.0)
  {
    const double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi3Cos = TMath::Cos(3.0*((double)iShift+1.0)*Psi3);
      const double Psi3Sin = TMath::Sin(3.0*((double)iShift+1.0)*Psi3);
      p_mTpcQ3ShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Cos);
      p_mTpcQ3ShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Sin);
    }
  }
}

void StTpcEpManager::fillTpcShiftFull()
{
  TVector2 Q1Vector = getQ1VecReCtrFull(); // 1st EP
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi1Cos = TMath::Cos(1.0*((double)iShift+1.0)*Psi1);
      const double Psi1Sin = TMath::Sin(1.0*((double)iShift+1.0)*Psi1);
      p_mTpcQ1ShiftCosFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
      p_mTpcQ1ShiftSinFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
    }
  }

  TVector2 Q2Vector = getQ2VecReCtrFull(); // 2nd EP
  if(Q2Vector.Mod() > 0.0)
  {
    const double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi2Cos = TMath::Cos(2.0*((double)iShift+1.0)*Psi2);
      const double Psi2Sin = TMath::Sin(2.0*((double)iShift+1.0)*Psi2);
      p_mTpcQ2ShiftCosFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Cos);
      p_mTpcQ2ShiftSinFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Sin);
    }
  }

  TVector2 Q3Vector = getQ3VecReCtrFull(); // 3rd EP
  if(Q3Vector.Mod() > 0.0)
  {
    const double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi3Cos = TMath::Cos(3.0*((double)iShift+1.0)*Psi3);
      const double Psi3Sin = TMath::Sin(3.0*((double)iShift+1.0)*Psi3);
      p_mTpcQ3ShiftCosFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Cos);
      p_mTpcQ3ShiftSinFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Sin);
    }
  }
}

void StTpcEpManager::writeTpcShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mTpcQ1ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ1ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ1ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ1ShiftSinWest[iVz][iShift]->Write();
      p_mTpcQ1ShiftCosFull[iVz][iShift]->Write();
      p_mTpcQ1ShiftSinFull[iVz][iShift]->Write();

      p_mTpcQ2ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ2ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinWest[iVz][iShift]->Write();
      p_mTpcQ2ShiftCosFull[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinFull[iVz][iShift]->Write();

      p_mTpcQ3ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ3ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinWest[iVz][iShift]->Write();
      p_mTpcQ3ShiftCosFull[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinFull[iVz][iShift]->Write();
    }
  }
}

void StTpcEpManager::readTpcShift()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftPar/file_TpcShiftPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mTpcQ1ShiftCos%dEastVz%d",iShift,iVz); // 1st EP
      p_mTpcQ1ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ1ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ1ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ1ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ1ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ1ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ1ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ1ShiftCos%dFullVz%d",iShift,iVz);
      p_mTpcQ1ShiftCosFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ1ShiftSin%dFullVz%d",iShift,iVz);
      p_mTpcQ1ShiftSinFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ2ShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mTpcQ2ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ2ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ2ShiftCos%dFullVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dFullVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dEastVz%d",iShift,iVz); // 3rd EP
      p_mTpcQ3ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dFullVz%d",iShift,iVz);
      p_mTpcQ3ShiftCosFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dFullVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

double StTpcEpManager::getPsi1ShiftEast(TVector2 Q1Vector)
{
  // TVector2 Q1Vector = getQ1VecReCtrEast();
  double Psi1Shift = -999.9;
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mTpcQ1ShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mTpcQ1ShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mTpcQ1ShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mTpcQ1ShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(1.0*((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(1.0*((double)iShift+1.0)*Psi1ReCtr));
    }

    double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StTpcEpManager::getPsi1ShiftWest(TVector2 Q1Vector)
{
  // TVector2 Q1Vector = getQ1VecReCtrWest();
  double Psi1Shift = -999.9;
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mTpcQ1ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mTpcQ1ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mTpcQ1ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mTpcQ1ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(1.0*((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(1.0*((double)iShift+1.0)*Psi1ReCtr));
    }

    double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StTpcEpManager::getPsi1ShiftFull(TVector2 Q1Vector)
{
  // TVector2 Q1Vector = getQ1VecReCtrFull();
  double Psi1Shift = -999.9;
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mTpcQ1ShiftCosFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mTpcQ1ShiftCosFull[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mTpcQ1ShiftSinFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mTpcQ1ShiftSinFull[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(1.0*((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(1.0*((double)iShift+1.0)*Psi1ReCtr));
    }

    double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StTpcEpManager::transPsi1(double Psi1)
{
  double Psi1Corr = Psi1;
  if(Psi1 >  TMath::Pi()) Psi1Corr = Psi1 - TMath::TwoPi();
  if(Psi1 < -TMath::Pi()) Psi1Corr = Psi1 + TMath::TwoPi();

  return Psi1Corr;
}

bool StTpcEpManager::isPsi1InRange(double Psi1)
{
  if(Psi1 < -TMath::Pi() || Psi1 > TMath::Pi())
  {
    return false;
  }

  return true;
}

double StTpcEpManager::getPsi2ShiftEast(TVector2 Q2Vector)
{
  // TVector2 Q2Vector = getQ2VecReCtrEast();
  double Psi2Shift = -999.9;
  if(Q2Vector.Mod() > 0.0)
  {
    const double Psi2ReCtr = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    double deltaPsi2 = 0.0;

    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mTpcQ2ShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mTpcQ2ShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mTpcQ2ShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mTpcQ2ShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi2 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*Psi2ReCtr)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*Psi2ReCtr));
    }

    double Psi2ShiftRaw = Psi2ReCtr + deltaPsi2/2.0;
    Psi2Shift = transPsi2(Psi2ShiftRaw);
  }

  return Psi2Shift;
}

double StTpcEpManager::getPsi2ShiftWest(TVector2 Q2Vector)
{
  // TVector2 Q2Vector = getQ2VecReCtrWest();
  double Psi2Shift = -999.9;
  if(Q2Vector.Mod() > 0.0)
  {
    const double Psi2ReCtr = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
    double deltaPsi2 = 0.0;

    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mTpcQ2ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mTpcQ2ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mTpcQ2ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mTpcQ2ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi2 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*Psi2ReCtr)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*Psi2ReCtr));
    }

    double Psi2ShiftRaw = Psi2ReCtr + deltaPsi2/2.0;
    Psi2Shift = transPsi2(Psi2ShiftRaw);
  }

  return Psi2Shift;
}

double StTpcEpManager::getPsi2ShiftFull(TVector2 Q2Vector)
{
  // TVector2 Q2Vector = getQ2VecReCtrFull();
  double Psi2Shift = -999.9;
  if(Q2Vector.Mod() > 0.0)
  {
    const double Psi2ReCtr = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
    double deltaPsi2 = 0.0;

    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mTpcQ2ShiftCosFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mTpcQ2ShiftCosFull[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mTpcQ2ShiftSinFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mTpcQ2ShiftSinFull[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi2 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*Psi2ReCtr)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*Psi2ReCtr));
    }

    double Psi2ShiftRaw = Psi2ReCtr + deltaPsi2/2.0;
    Psi2Shift = transPsi2(Psi2ShiftRaw);
  }

  return Psi2Shift;
}

double StTpcEpManager::transPsi2(double Psi2)
{
  double Psi2Corr = Psi2;
  if(Psi2 >  TMath::Pi()/2.0) Psi2Corr = Psi2 - TMath::Pi();
  if(Psi2 < -TMath::Pi()/2.0) Psi2Corr = Psi2 + TMath::Pi();

  return Psi2Corr;
}

bool StTpcEpManager::isPsi2InRange(double Psi2)
{
  if(Psi2 < -TMath::Pi()/2.0 || Psi2 > TMath::Pi()/2.0)
  {
    return false;
  }

  return true;
}

double StTpcEpManager::getPsi3ShiftEast(TVector2 Q3Vector)
{
  // TVector2 Q3Vector = getQ3VecReCtrEast();
  double Psi3Shift = -999.9;
  if(Q3Vector.Mod() > 0.0)
  {
    const double Psi3ReCtr = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    double deltaPsi3 = 0.0;

    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mTpcQ3ShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mTpcQ3ShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mTpcQ3ShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mTpcQ3ShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi3 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(3.0*((double)iShift+1.0)*Psi3ReCtr)+meanCos*TMath::Sin(3.0*((double)iShift+1.0)*Psi3ReCtr));
    }

    double Psi3ShiftRaw = Psi3ReCtr + deltaPsi3/3.0;
    Psi3Shift = transPsi3(Psi3ShiftRaw);
  }

  return Psi3Shift;
}

double StTpcEpManager::getPsi3ShiftWest(TVector2 Q3Vector)
{
  // TVector2 Q3Vector = getQ3VecReCtrWest();
  double Psi3Shift = -999.9;
  if(Q3Vector.Mod() > 0.0)
  {
    const double Psi3ReCtr = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0;
    double deltaPsi3 = 0.0;

    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mTpcQ3ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mTpcQ3ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mTpcQ3ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mTpcQ3ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi3 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(3.0*((double)iShift+1.0)*Psi3ReCtr)+meanCos*TMath::Sin(3.0*((double)iShift+1.0)*Psi3ReCtr));
    }

    double Psi3ShiftRaw = Psi3ReCtr + deltaPsi3/3.0;
    Psi3Shift = transPsi3(Psi3ShiftRaw);
  }

  return Psi3Shift;
}

double StTpcEpManager::getPsi3ShiftFull(TVector2 Q3Vector)
{
  // TVector2 Q3Vector = getQ3VecReCtrFull();
  double Psi3Shift = -999.9;
  if(Q3Vector.Mod() > 0.0)
  {
    const double Psi3ReCtr = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0;
    double deltaPsi3 = 0.0;

    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mTpcQ3ShiftCosFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mTpcQ3ShiftCosFull[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mTpcQ3ShiftSinFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mTpcQ3ShiftSinFull[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi3 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(3.0*((double)iShift+1.0)*Psi3ReCtr)+meanCos*TMath::Sin(3.0*((double)iShift+1.0)*Psi3ReCtr));
    }

    double Psi3ShiftRaw = Psi3ReCtr + deltaPsi3/3.0;
    Psi3Shift = transPsi3(Psi3ShiftRaw);
  }

  return Psi3Shift;
}

double StTpcEpManager::transPsi3(double Psi3)
{
  double Psi3Corr = Psi3;
  if(Psi3 >  TMath::Pi()/3.0) Psi3Corr = Psi3 - TMath::TwoPi()/3.0;
  if(Psi3 < -TMath::Pi()/3.0) Psi3Corr = Psi3 + TMath::TwoPi()/3.0;

  return Psi3Corr;
}

bool StTpcEpManager::isPsi3InRange(double Psi3)
{
  if(Psi3 < -TMath::Pi()/3.0 || Psi3 > TMath::Pi()/3.0)
  {
    return false;
  }

  return true;
}
//---------------------------------------------------------------------------------
// Sub Event Plane Resolution
void StTpcEpManager::initTpcResolution()
{
  p_mTpcSubEp1Res = new TProfile("p_mTpcSubEp1Res","p_mTpcSubEp1Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  p_mTpcSubEp2Res = new TProfile("p_mTpcSubEp2Res","p_mTpcSubEp2Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  p_mTpcSubEp3Res = new TProfile("p_mTpcSubEp3Res","p_mTpcSubEp3Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
}

void StTpcEpManager::fillTpcResolution(double Psi1East, double Psi1West, double Psi2East, double Psi2West, double Psi3East, double Psi3West)
{
  double res1Sub = TMath::Cos(1.0*(Psi1West-Psi1East));
  p_mTpcSubEp1Res->Fill((double)mCent9,res1Sub);

  double res2Sub = TMath::Cos(2.0*(Psi2West-Psi2East));
  p_mTpcSubEp2Res->Fill((double)mCent9,res2Sub);

  double res3Sub = TMath::Cos(3.0*(Psi3West-Psi3East));
  p_mTpcSubEp3Res->Fill((double)mCent9,res3Sub);
}

void StTpcEpManager::writeTpcResolution()
{
  p_mTpcSubEp1Res->Write();
  p_mTpcSubEp2Res->Write();
  p_mTpcSubEp3Res->Write();
}

void StTpcEpManager::readTpcResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_TpcEpResolution_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mTpcSubEp1Res = (TProfile*)file_mResolution->Get("p_mTpcSubEp1Res");
  p_mTpcSubEp2Res = (TProfile*)file_mResolution->Get("p_mTpcSubEp2Res");
  p_mTpcSubEp3Res = (TProfile*)file_mResolution->Get("p_mTpcSubEp3Res");
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    mTpcSubEp1ResVal[iCent] = 0.0;
    mTpcSubEp1ResErr[iCent] = 0.0;
    mTpcSubEp2ResVal[iCent] = 0.0;
    mTpcSubEp2ResErr[iCent] = 0.0;
    mTpcSubEp3ResVal[iCent] = 0.0;
    mTpcSubEp3ResErr[iCent] = 0.0;
  }

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    double valRes1Sub = -999.9; // 1st EP
    double errRes1Sub = 1.0;
    double valRes1Raw = p_mTpcSubEp1Res->GetBinContent(iCent+1);
    double errRes1Raw = p_mTpcSubEp1Res->GetBinError(iCent+1);
    if(valRes1Raw > 0)
    {
      valRes1Sub = TMath::Sqrt(valRes1Raw);
      errRes1Sub = errRes1Raw/(2.0*valRes1Sub);
    }
    mTpcSubEp1ResVal[iCent] = valRes1Sub;
    mTpcSubEp1ResErr[iCent] = errRes1Sub;

    double valRes2Sub = -999.9; // 2nd EP
    double errRes2Sub = 1.0;
    double valRes2Raw = p_mTpcSubEp2Res->GetBinContent(iCent+1);
    double errRes2Raw = p_mTpcSubEp2Res->GetBinError(iCent+1);
    if(valRes2Raw > 0)
    {
      valRes2Sub = TMath::Sqrt(valRes2Raw);
      errRes2Sub = errRes2Raw/(2.0*valRes2Sub);
    }
    mTpcSubEp2ResVal[iCent] = valRes2Sub;
    mTpcSubEp2ResErr[iCent] = errRes2Sub;

    double valRes3Sub = -999.9; // 3rd EP
    double errRes3Sub = 1.0;
    double valRes3Raw = p_mTpcSubEp3Res->GetBinContent(iCent+1);
    double errRes3Raw = p_mTpcSubEp3Res->GetBinError(iCent+1);
    if(valRes3Raw > 0)
    {
      valRes3Sub = TMath::Sqrt(valRes3Raw);
      errRes3Sub = errRes3Raw/(2.0*valRes3Sub);
    }
    mTpcSubEp3ResVal[iCent] = valRes3Sub;
    mTpcSubEp3ResErr[iCent] = errRes3Sub;
  }
  file_mResolution->Close();
}

double StTpcEpManager::getTpcSubEp1ResVal(int cent9)
{
  return mTpcSubEp1ResVal[cent9];
}

double StTpcEpManager::getTpcSubEp1ResErr(int cent9)
{
  return mTpcSubEp1ResErr[cent9];
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
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string proName = Form("p_mTpcSubEpV1Cent%d",iCent);
    p_mTpcSubEpV1[iCent] = new TProfile(proName.c_str(),proName.c_str(),50,0.0,10.0);
    proName = Form("p_mTpcSubEpV2Cent%d",iCent);
    p_mTpcSubEpV2[iCent] = new TProfile(proName.c_str(),proName.c_str(),50,0.0,10.0);
    proName = Form("p_mTpcSubEpV3Cent%d",iCent);
    p_mTpcSubEpV3[iCent] = new TProfile(proName.c_str(),proName.c_str(),50,0.0,10.0);
  }
}

void StTpcEpManager::fillTpcSubEpV1(double pt, double v1, double reweight)
{
    p_mTpcSubEpV1[mCent9]->Fill(pt, v1, reweight);
}

void StTpcEpManager::fillTpcSubEpV2(double pt, double v2, double reweight)
{
    p_mTpcSubEpV2[mCent9]->Fill(pt, v2, reweight);
}

void StTpcEpManager::fillTpcSubEpV3(double pt, double v3, double reweight)
{
    p_mTpcSubEpV3[mCent9]->Fill(pt, v3, reweight);
}

void StTpcEpManager::writeTpcSubEpFlow()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    p_mTpcSubEpV1[iCent]->Write();
    p_mTpcSubEpV2[iCent]->Write();
    p_mTpcSubEpV3[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
// QVector
TVector2 StTpcEpManager::getQ1VecRawEast()
{
  return -1.0*v_mQ1RawEast; // flip the sign for Q1VecEast
}

TVector2 StTpcEpManager::getQ1VecRawWest()
{
  return v_mQ1RawWest;
}

TVector2 StTpcEpManager::getQ1VecRawFull()
{
  return v_mQ1RawFull;
}

TVector2 StTpcEpManager::getQ1VecReCtrEast()
{
  return -1.0*v_mQ1ReCtrEast; // flip the sign for Q1VecEast
}

TVector2 StTpcEpManager::getQ1VecReCtrWest()
{
  return v_mQ1ReCtrWest;
}

TVector2 StTpcEpManager::getQ1VecReCtrFull()
{
  return v_mQ1ReCtrFull;
}

double StTpcEpManager::getPsi1RawEast()
{
  double Psi1Raw = -999.9;
  TVector2 Q1Vector = getQ1VecRawEast();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Raw = transPsi1(Psi1);
  }

  return Psi1Raw;
}

double StTpcEpManager::getPsi1RawWest()
{
  double Psi1Raw = -999.9;
  TVector2 Q1Vector = getQ1VecRawWest();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Raw = transPsi1(Psi1);
  }

  return Psi1Raw;
}

double StTpcEpManager::getPsi1RawFull()
{
  double Psi1Raw = -999.9;
  TVector2 Q1Vector = getQ1VecRawFull();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Raw = transPsi1(Psi1);
  }

  return Psi1Raw;
}

double StTpcEpManager::getPsi1ReCtrEast()
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector   = getQ1VecReCtrEast();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

double StTpcEpManager::getPsi1ReCtrWest()
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector   = getQ1VecReCtrWest();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

double StTpcEpManager::getPsi1ReCtrFull()
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector   = getQ1VecReCtrFull();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

TVector2 StTpcEpManager::getQ2VecRawEast()
{
  return v_mQ2RawEast;
}

TVector2 StTpcEpManager::getQ2VecRawWest()
{
  return v_mQ2RawWest;
}

TVector2 StTpcEpManager::getQ2VecRawFull()
{
  return v_mQ2RawFull;
}

TVector2 StTpcEpManager::getQ2VecReCtrEast()
{
  return v_mQ2ReCtrEast;
}

TVector2 StTpcEpManager::getQ2VecReCtrWest()
{
  return v_mQ2ReCtrWest;
}

TVector2 StTpcEpManager::getQ2VecReCtrFull()
{
  return v_mQ2ReCtrFull;
}

double StTpcEpManager::getPsi2RawEast()
{
  double Psi2Raw = -999.9;
  TVector2 Q2Vector = getQ2VecRawEast();
  if(Q2Vector.Mod() > 0.0)
  {
    double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    Psi2Raw = transPsi2(Psi2);
  }

  return Psi2Raw;
}

double StTpcEpManager::getPsi2RawWest()
{
  double Psi2Raw = -999.9;
  TVector2 Q2Vector = getQ2VecRawWest();
  if(Q2Vector.Mod() > 0.0)
  {
    double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    Psi2Raw = transPsi2(Psi2);
  }

  return Psi2Raw;
}

double StTpcEpManager::getPsi2RawFull()
{
  double Psi2Raw = -999.9;
  TVector2 Q2Vector = getQ2VecRawFull();
  if(Q2Vector.Mod() > 0.0)
  {
    double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    Psi2Raw = transPsi2(Psi2);
  }

  return Psi2Raw;
}

double StTpcEpManager::getPsi2ReCtrEast()
{
  double Psi2ReCtr = -999.9;
  TVector2 Q2Vector   = getQ2VecReCtrEast();
  if(Q2Vector.Mod() > 0.0)
  {
    double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    Psi2ReCtr = transPsi2(Psi2);
  }

  return Psi2ReCtr;
}

double StTpcEpManager::getPsi2ReCtrWest()
{
  double Psi2ReCtr = -999.9;
  TVector2 Q2Vector   = getQ2VecReCtrWest();
  if(Q2Vector.Mod() > 0.0)
  {
    double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    Psi2ReCtr = transPsi2(Psi2);
  }

  return Psi2ReCtr;
}

double StTpcEpManager::getPsi2ReCtrFull()
{
  double Psi2ReCtr = -999.9;
  TVector2 Q2Vector   = getQ2VecReCtrFull();
  if(Q2Vector.Mod() > 0.0)
  {
    double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
    Psi2ReCtr = transPsi2(Psi2);
  }

  return Psi2ReCtr;
}

TVector2 StTpcEpManager::getQ3VecRawEast()
{
  return v_mQ3RawEast;
}

TVector2 StTpcEpManager::getQ3VecRawWest()
{
  return v_mQ3RawWest;
}

TVector2 StTpcEpManager::getQ3VecRawFull()
{
  return v_mQ3RawFull;
}

TVector2 StTpcEpManager::getQ3VecReCtrEast()
{
  return v_mQ3ReCtrEast;
}

TVector2 StTpcEpManager::getQ3VecReCtrWest()
{
  return v_mQ3ReCtrWest;
}

TVector2 StTpcEpManager::getQ3VecReCtrFull()
{
  return v_mQ3ReCtrFull;
}

double StTpcEpManager::getPsi3RawEast()
{
  double Psi3Raw = -999.9;
  TVector2 Q3Vector = getQ3VecRawEast();
  if(Q3Vector.Mod() > 0.0)
  {
    double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    Psi3Raw = transPsi3(Psi3);
  }

  return Psi3Raw;
}

double StTpcEpManager::getPsi3RawWest()
{
  double Psi3Raw = -999.9;
  TVector2 Q3Vector = getQ3VecRawWest();
  if(Q3Vector.Mod() > 0.0)
  {
    double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    Psi3Raw = transPsi3(Psi3);
  }

  return Psi3Raw;
}

double StTpcEpManager::getPsi3RawFull()
{
  double Psi3Raw = -999.9;
  TVector2 Q3Vector = getQ3VecRawFull();
  if(Q3Vector.Mod() > 0.0)
  {
    double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    Psi3Raw = transPsi3(Psi3);
  }

  return Psi3Raw;
}

double StTpcEpManager::getPsi3ReCtrEast()
{
  double Psi3ReCtr = -999.9;
  TVector2 Q3Vector   = getQ3VecReCtrEast();
  if(Q3Vector.Mod() > 0.0)
  {
    double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    Psi3ReCtr = transPsi3(Psi3);
  }

  return Psi3ReCtr;
}

double StTpcEpManager::getPsi3ReCtrWest()
{
  double Psi3ReCtr = -999.9;
  TVector2 Q3Vector   = getQ3VecReCtrWest();
  if(Q3Vector.Mod() > 0.0)
  {
    double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    Psi3ReCtr = transPsi3(Psi3);
  }

  return Psi3ReCtr;
}

double StTpcEpManager::getPsi3ReCtrFull()
{
  double Psi3ReCtr = -999.9;
  TVector2 Q3Vector   = getQ3VecReCtrFull();
  if(Q3Vector.Mod() > 0.0)
  {
    double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
    Psi3ReCtr = transPsi3(Psi3);
  }

  return Psi3ReCtr;
}

int StTpcEpManager::getNumTrkRawEast()
{
  return mQCouRawEast;
}

int StTpcEpManager::getNumTrkRawWest()
{
  return mQCouRawWest;
}

int StTpcEpManager::getNumTrkRawFull()
{
  return mQCouRawFull;
}

int StTpcEpManager::getNumTrkReCtrEast()
{
  return mQCouReCtrEast;
}

int StTpcEpManager::getNumTrkReCtrWest()
{
  return mQCouReCtrWest;
}

int StTpcEpManager::getNumTrkReCtrFull()
{
  return mQCouReCtrFull;
}
//---------------------------------------------------------------------------------
// raw EP
void StTpcEpManager::initTpcSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp1RawEastCent%d",iCent); // 1st EP
    h_mTpcEp1RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1RawWestCent%d",iCent);
    h_mTpcEp1RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1RawFullCent%d",iCent);
    h_mTpcEp1RawFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1RawCorrCent%d",iCent);
    h_mTpcEp1RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());

    histName = Form("h_mTpcEp2RawEastCent%d",iCent); // 2nd EP
    h_mTpcEp2RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,720,-1.0*TMath::Pi(),1.0*TMath::Pi());
    histName = Form("h_mTpcEp2RawWestCent%d",iCent);
    h_mTpcEp2RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,720,-1.0*TMath::Pi(),1.0*TMath::Pi());
    histName = Form("h_mTpcEp2RawFullCent%d",iCent);
    h_mTpcEp2RawFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,720,-1.0*TMath::Pi(),1.0*TMath::Pi());
    histName = Form("h_mTpcEp2RawCorrCent%d",iCent);
    h_mTpcEp2RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),180,-1.0*TMath::Pi(),1.0*TMath::Pi(),180,-1.0*TMath::Pi(),1.0*TMath::Pi());

    histName = Form("h_mTpcEp3RawEastCent%d",iCent); // 3rd EP
    h_mTpcEp3RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-0.5*TMath::Pi(),0.5*TMath::Pi());
    histName = Form("h_mTpcEp3RawWestCent%d",iCent);
    h_mTpcEp3RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-0.5*TMath::Pi(),0.5*TMath::Pi());
    histName = Form("h_mTpcEp3RawFullCent%d",iCent);
    h_mTpcEp3RawFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-0.5*TMath::Pi(),0.5*TMath::Pi());
    histName = Form("h_mTpcEp3RawCorrCent%d",iCent);
    h_mTpcEp3RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-0.5*TMath::Pi(),0.5*TMath::Pi(),135,-0.5*TMath::Pi(),0.5*TMath::Pi());
  }
}

void StTpcEpManager::fillTpcSubEpRaw(double Psi1East, double Psi1West, double Psi1Full, double Psi2East, double Psi2West, double Psi2Full, double Psi3East, double Psi3West, double Psi3Full)
{
  h_mTpcEp1RawEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mTpcEp1RawWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mTpcEp1RawFull[mCent9]->Fill(mRunIndex,Psi1Full);
  h_mTpcEp1RawCorr[mCent9]->Fill(Psi1East,Psi1West);

  h_mTpcEp2RawEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mTpcEp2RawWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mTpcEp2RawFull[mCent9]->Fill(mRunIndex,Psi2Full);
  h_mTpcEp2RawCorr[mCent9]->Fill(Psi2East,Psi2West);

  h_mTpcEp3RawEast[mCent9]->Fill(mRunIndex,Psi3East);
  h_mTpcEp3RawWest[mCent9]->Fill(mRunIndex,Psi3West);
  h_mTpcEp3RawFull[mCent9]->Fill(mRunIndex,Psi3Full);
  h_mTpcEp3RawCorr[mCent9]->Fill(Psi3East,Psi3West);
}

void StTpcEpManager::writeTpcSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp1RawEast[iCent]->Write();
    h_mTpcEp1RawWest[iCent]->Write();
    h_mTpcEp1RawFull[iCent]->Write();
    h_mTpcEp1RawCorr[iCent]->Write();

    h_mTpcEp2RawEast[iCent]->Write();
    h_mTpcEp2RawWest[iCent]->Write();
    h_mTpcEp2RawFull[iCent]->Write();
    h_mTpcEp2RawCorr[iCent]->Write();

    h_mTpcEp3RawEast[iCent]->Write();
    h_mTpcEp3RawWest[iCent]->Write();
    h_mTpcEp3RawFull[iCent]->Write();
    h_mTpcEp3RawCorr[iCent]->Write();
  }
}

// recenter EP
void StTpcEpManager::initTpcSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp1ReCtrEastCent%d",iCent); // 1st EP
    h_mTpcEp1ReCtrEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1ReCtrWestCent%d",iCent);
    h_mTpcEp1ReCtrWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1ReCtrFullCent%d",iCent);
    h_mTpcEp1ReCtrFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1ReCtrCorrCent%d",iCent);
    h_mTpcEp1ReCtrCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());

    histName = Form("h_mTpcEp2ReCtrEastCent%d",iCent); // 2nd EP
    h_mTpcEp2ReCtrEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,720,-1.0*TMath::Pi(),1.0*TMath::Pi());
    histName = Form("h_mTpcEp2ReCtrWestCent%d",iCent);
    h_mTpcEp2ReCtrWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,720,-1.0*TMath::Pi(),1.0*TMath::Pi());
    histName = Form("h_mTpcEp2ReCtrFullCent%d",iCent);
    h_mTpcEp2ReCtrFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,720,-1.0*TMath::Pi(),1.0*TMath::Pi());
    histName = Form("h_mTpcEp2ReCtrCorrCent%d",iCent);
    h_mTpcEp2ReCtrCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),180,-1.0*TMath::Pi(),1.0*TMath::Pi(),180,-1.0*TMath::Pi(),1.0*TMath::Pi());

    histName = Form("h_mTpcEp3ReCtrEastCent%d",iCent); // 3rd EP
    h_mTpcEp3ReCtrEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-0.5*TMath::Pi(),0.5*TMath::Pi());
    histName = Form("h_mTpcEp3ReCtrWestCent%d",iCent);
    h_mTpcEp3ReCtrWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-0.5*TMath::Pi(),0.5*TMath::Pi());
    histName = Form("h_mTpcEp3ReCtrFullCent%d",iCent);
    h_mTpcEp3ReCtrFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-0.5*TMath::Pi(),0.5*TMath::Pi());
    histName = Form("h_mTpcEp3ReCtrCorrCent%d",iCent);
    h_mTpcEp3ReCtrCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-0.5*TMath::Pi(),0.5*TMath::Pi(),135,-0.5*TMath::Pi(),0.5*TMath::Pi());
  }
}

void StTpcEpManager::fillTpcSubEpReCtr(double Psi1East, double Psi1West, double Psi1Full, double Psi2East, double Psi2West, double Psi2Full, double Psi3East, double Psi3West, double Psi3Full)
{
  h_mTpcEp1ReCtrEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mTpcEp1ReCtrWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mTpcEp1ReCtrFull[mCent9]->Fill(mRunIndex,Psi1Full);
  h_mTpcEp1ReCtrCorr[mCent9]->Fill(Psi1East,Psi1West);

  h_mTpcEp2ReCtrEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mTpcEp2ReCtrWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mTpcEp2ReCtrFull[mCent9]->Fill(mRunIndex,Psi2Full);
  h_mTpcEp2ReCtrCorr[mCent9]->Fill(Psi2East,Psi2West);

  h_mTpcEp3ReCtrEast[mCent9]->Fill(mRunIndex,Psi3East);
  h_mTpcEp3ReCtrWest[mCent9]->Fill(mRunIndex,Psi3West);
  h_mTpcEp3ReCtrFull[mCent9]->Fill(mRunIndex,Psi3Full);
  h_mTpcEp3ReCtrCorr[mCent9]->Fill(Psi3East,Psi3West);
}

void StTpcEpManager::writeTpcSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp1ReCtrEast[iCent]->Write();
    h_mTpcEp1ReCtrWest[iCent]->Write();
    h_mTpcEp1ReCtrFull[iCent]->Write();
    h_mTpcEp1ReCtrCorr[iCent]->Write();

    h_mTpcEp2ReCtrEast[iCent]->Write();
    h_mTpcEp2ReCtrWest[iCent]->Write();
    h_mTpcEp2ReCtrFull[iCent]->Write();
    h_mTpcEp2ReCtrCorr[iCent]->Write();

    h_mTpcEp3ReCtrEast[iCent]->Write();
    h_mTpcEp3ReCtrWest[iCent]->Write();
    h_mTpcEp3ReCtrFull[iCent]->Write();
    h_mTpcEp3ReCtrCorr[iCent]->Write();
  }
}

// shift EP
void StTpcEpManager::initTpcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp1ShiftEastCent%d",iCent); // 1st EP
    h_mTpcEp1ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1ShiftWestCent%d",iCent);
    h_mTpcEp1ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1ShiftFullCent%d",iCent);
    h_mTpcEp1ShiftFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1ShiftCorrCent%d",iCent);
    h_mTpcEp1ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());

    histName = Form("h_mTpcEp2ShiftEastCent%d",iCent); // 2nd EP
    h_mTpcEp2ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,720,-1.0*TMath::Pi(),1.0*TMath::Pi());
    histName = Form("h_mTpcEp2ShiftWestCent%d",iCent);
    h_mTpcEp2ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,720,-1.0*TMath::Pi(),1.0*TMath::Pi());
    histName = Form("h_mTpcEp2ShiftFullCent%d",iCent);
    h_mTpcEp2ShiftFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,720,-1.0*TMath::Pi(),1.0*TMath::Pi());
    histName = Form("h_mTpcEp2ShiftCorrCent%d",iCent);
    h_mTpcEp2ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),180,-1.0*TMath::Pi(),1.0*TMath::Pi(),180,-1.0*TMath::Pi(),1.0*TMath::Pi());

    histName = Form("h_mTpcEp3ShiftEastCent%d",iCent); // 3rd EP
    h_mTpcEp3ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-0.5*TMath::Pi(),0.5*TMath::Pi());
    histName = Form("h_mTpcEp3ShiftWestCent%d",iCent);
    h_mTpcEp3ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-0.5*TMath::Pi(),0.5*TMath::Pi());
    histName = Form("h_mTpcEp3ShiftFullCent%d",iCent);
    h_mTpcEp3ShiftFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-0.5*TMath::Pi(),0.5*TMath::Pi());
    histName = Form("h_mTpcEp3ShiftCorrCent%d",iCent);
    h_mTpcEp3ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-0.5*TMath::Pi(),0.5*TMath::Pi(),135,-0.5*TMath::Pi(),0.5*TMath::Pi());
  }
}

void StTpcEpManager::fillTpcSubEpShift(double Psi1East, double Psi1West, double Psi1Full, double Psi2East, double Psi2West, double Psi2Full, double Psi3East, double Psi3West, double Psi3Full)
{
  h_mTpcEp1ShiftEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mTpcEp1ShiftWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mTpcEp1ShiftFull[mCent9]->Fill(mRunIndex,Psi1Full);
  h_mTpcEp1ShiftCorr[mCent9]->Fill(Psi1East,Psi1West);

  h_mTpcEp2ShiftEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mTpcEp2ShiftWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mTpcEp2ShiftFull[mCent9]->Fill(mRunIndex,Psi2Full);
  h_mTpcEp2ShiftCorr[mCent9]->Fill(Psi2East,Psi2West);

  h_mTpcEp3ShiftEast[mCent9]->Fill(mRunIndex,Psi3East);
  h_mTpcEp3ShiftWest[mCent9]->Fill(mRunIndex,Psi3West);
  h_mTpcEp3ShiftFull[mCent9]->Fill(mRunIndex,Psi3Full);
  h_mTpcEp3ShiftCorr[mCent9]->Fill(Psi3East,Psi3West);
}

void StTpcEpManager::writeTpcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp1ShiftEast[iCent]->Write();
    h_mTpcEp1ShiftWest[iCent]->Write();
    h_mTpcEp1ShiftFull[iCent]->Write();
    h_mTpcEp1ShiftCorr[iCent]->Write();

    h_mTpcEp2ShiftEast[iCent]->Write();
    h_mTpcEp2ShiftWest[iCent]->Write();
    h_mTpcEp2ShiftFull[iCent]->Write();
    h_mTpcEp2ShiftCorr[iCent]->Write();

    h_mTpcEp3ShiftEast[iCent]->Write();
    h_mTpcEp3ShiftWest[iCent]->Write();
    h_mTpcEp3ShiftFull[iCent]->Write();
    h_mTpcEp3ShiftCorr[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
