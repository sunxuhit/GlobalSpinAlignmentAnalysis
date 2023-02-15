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
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StEventPlaneMaker/StEpdEpManager.h"

ClassImp(StEpdEpManager)

//---------------------------------------------------------------------------------
StEpdEpManager::StEpdEpManager(int beamType) : mType(beamType)
{
  // mEnergy = energy;
  clearEpdEpManager();
  mEpdGeom = new StEpdGeom();
}

StEpdEpManager::~StEpdEpManager()
{
  /* */
}
//---------------------------------------------------------------------------------
void StEpdEpManager::clearEpdEpManager()
{
  mCent9 = -1;
  mRunIndex = -1;
  mVzBin = -1;

  // mQCouRawEast = 0; // raw EP
  // mQCouRawWest = 0;
  // v_mQ2RawWest.Set(0.0,0.0); // 2nd raw EP
  // v_mQ2RawEast.Set(0.0,0.0);
  // v_mQ3RawEast.Set(0.0,0.0); // 3rd raw EP
  // v_mQ3RawWest.Set(0.0,0.0);

  // mQCouReCenterEast = 0; // recenter EP
  // mQCouReCenterWest = 0;
  // v_mQ2ReCenterEast.Set(0.0,0.0); // 2nd recenter EP
  // v_mQ2ReCenterWest.Set(0.0,0.0);
  // v_mQ3ReCenterEast.Set(0.0,0.0);
  // v_mQ3ReCenterWest.Set(0.0,0.0); // 3rd recenter EP
}

void StEpdEpManager::initEpdEpManager(int cent9, int runIndex, int vzBin)
{
  mCent9    = cent9;
  mRunIndex = runIndex;
  mVzBin    = vzBin;
}
//---------------------------------------------------------------------------------
// Utilities
TVector2 StEpdEpManager::calq1Vector(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  TVector2 q1Vector(0.0,0.0);
  TVector3 EpdPoint = mEpdGeom->RandomPointOnTile(EpdHit->id()); // get a random position within the tile
  TVector3 StraightLine = EpdPoint-primVtx;
  const double phi = StraightLine.Phi(); // -pi to pi
					 //
  const double q1x = TMath::Cos(1.0*phi);
  const double q1y = TMath::Sin(1.0*phi);
  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

double StEpdEpManager::getTileWeight(StPicoEpdHit *picoEpdHit)
{
  const double nMip = picoEpdHit->nMIP();
  double tileWgt = (nMip < anaUtils::mMipEpdEpMax[mType]) ? nMip : mMipEpdEpMax[mType];
  if(nMip < anaUtils::mMipEpdEpMin[mType]) tileWgt = 0.0;

  return tileWgt;
}

double StEpdEpManager::getPhiWeight(StPicoEpdHit *picoEpdHit)
{
  double phiWgt = 1.0; // TODO: add real phi weight

  return phiWgt;
}

double StEpdEpManager::getEtaWeight(StPicoEpdHit *picoEpdHit)
{
  double etaWgt = 1.0; // TODO: add real eta weight

  return etaWgt;
}

#if 0
void StEpdEpManager::addHitRawEast(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2RawEast += wgt*calq2Vector(picoTrack);
  v_mQ3RawEast += wgt*calq3Vector(picoTrack);
  mQCouRawEast++;
}

void StEpdEpManager::addHitRawWest(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2RawWest += wgt*calq2Vector(picoTrack);
  v_mQ3RawWest += wgt*calq3Vector(picoTrack);
  mQCouRawWest++;
}

void StEpdEpManager::addHitReCenterEast(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2ReCenterEast += wgt*(calq2Vector(picoTrack) - getq2VecReCenterEast());
  v_mQ3ReCenterEast += wgt*(calq3Vector(picoTrack) - getq3VecReCenterEast());
  mQCouReCenterEast++;
}

void StEpdEpManager::addHitReCenterWest(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2ReCenterWest += wgt*(calq2Vector(picoTrack) - getq2VecReCenterWest());
  v_mQ3ReCenterWest += wgt*(calq3Vector(picoTrack) - getq3VecReCenterWest());
  mQCouReCenterWest++;
}
//---------------------------------------------------------------------------------
// ReCenterPar Correction
void StEpdEpManager::initEpdReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mEpdQ2ReCenterXEastVz%d",iVz); // 2nd EP
    p_mEpdQ2ReCenterXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mEpdQ2ReCenterYEastVz%d",iVz);
    p_mEpdQ2ReCenterYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mEpdQ2ReCenterXWestVz%d",iVz);
    p_mEpdQ2ReCenterXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mEpdQ2ReCenterYWestVz%d",iVz);
    p_mEpdQ2ReCenterYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mEpdQ3ReCenterXEastVz%d",iVz); // 3rd EP
    p_mEpdQ3ReCenterXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mEpdQ3ReCenterYEastVz%d",iVz);
    p_mEpdQ3ReCenterYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mEpdQ3ReCenterXWestVz%d",iVz);
    p_mEpdQ3ReCenterXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mEpdQ3ReCenterYWestVz%d",iVz);
    p_mEpdQ3ReCenterYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  }
}

void StEpdEpManager::fillEpdReCenterEast(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);

  TVector2 q2Vector = calq2Vector(picoTrack);
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();
  p_mEpdQ2ReCenterXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mEpdQ2ReCenterYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);

  TVector2 q3Vector = calq3Vector(picoTrack);
  const double q3x  = q3Vector.X();
  const double q3y  = q3Vector.Y();
  p_mEpdQ3ReCenterXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3x,wgt);
  p_mEpdQ3ReCenterYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3y,wgt);
}

void StEpdEpManager::fillEpdReCenterWest(StPicoTrack* picoTrack)
{
  const double wgt  = getWeight(picoTrack);

  TVector2 q2Vector = calq2Vector(picoTrack);
  const double q2x  = q2Vector.X();
  const double q2y  = q2Vector.Y();
  p_mEpdQ2ReCenterXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2x,wgt);
  p_mEpdQ2ReCenterYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q2y,wgt);

  TVector2 q3Vector = calq3Vector(picoTrack);
  const double q3x  = q3Vector.X();
  const double q3y  = q3Vector.Y();
  p_mEpdQ3ReCenterXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3x,wgt);
  p_mEpdQ3ReCenterYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q3y,wgt);
}

void StEpdEpManager::writeEpdReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mEpdQ2ReCenterXEast[iVz]->Write();
    p_mEpdQ2ReCenterYEast[iVz]->Write();
    p_mEpdQ2ReCenterXWest[iVz]->Write();
    p_mEpdQ2ReCenterYWest[iVz]->Write();

    p_mEpdQ3ReCenterXEast[iVz]->Write();
    p_mEpdQ3ReCenterYEast[iVz]->Write();
    p_mEpdQ3ReCenterXWest[iVz]->Write();
    p_mEpdQ3ReCenterYWest[iVz]->Write();
  }
}

void StEpdEpManager::readEpdReCenter()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_EpdReCenterPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mReCenterPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mEpdQ2ReCenterXEastVz%d",iVz); // 2nd EP
    p_mEpdQ2ReCenterXEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mEpdQ2ReCenterYEastVz%d",iVz);
    p_mEpdQ2ReCenterYEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mEpdQ2ReCenterXWestVz%d",iVz);
    p_mEpdQ2ReCenterXWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mEpdQ2ReCenterYWestVz%d",iVz);
    p_mEpdQ2ReCenterYWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mEpdQ3ReCenterXEastVz%d",iVz); // 3rd EP
    p_mEpdQ3ReCenterXEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mEpdQ3ReCenterYEastVz%d",iVz);
    p_mEpdQ3ReCenterYEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mEpdQ3ReCenterXWestVz%d",iVz);
    p_mEpdQ3ReCenterXWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mEpdQ3ReCenterYWestVz%d",iVz);
    p_mEpdQ3ReCenterYWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
  }
}

TVector2 StEpdEpManager::getq2VecReCenterEast()
{
  const int binX   = p_mEpdQ2ReCenterXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mEpdQ2ReCenterXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mEpdQ2ReCenterYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mEpdQ2ReCenterYEast[mVzBin]->GetBinContent(binY);

  TVector2 q2Vector(0.0,0.0);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StEpdEpManager::getq2VecReCenterWest()
{
  const int binX   = p_mEpdQ2ReCenterXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2x = p_mEpdQ2ReCenterXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mEpdQ2ReCenterYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q2y = p_mEpdQ2ReCenterYWest[mVzBin]->GetBinContent(binY);

  TVector2 q2Vector(0.0,0.0);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StEpdEpManager::getq3VecReCenterEast()
{
  const int binX   = p_mEpdQ3ReCenterXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3x = p_mEpdQ3ReCenterXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mEpdQ3ReCenterYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3y = p_mEpdQ3ReCenterYEast[mVzBin]->GetBinContent(binY);

  TVector2 q3Vector(0.0,0.0);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
}

TVector2 StEpdEpManager::getq3VecReCenterWest()
{
  const int binX   = p_mEpdQ3ReCenterXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3x = p_mEpdQ3ReCenterXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mEpdQ3ReCenterYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q3y = p_mEpdQ3ReCenterYWest[mVzBin]->GetBinContent(binY);

  TVector2 q3Vector(0.0,0.0);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
}
//---------------------------------------------------------------------------------
// Shift Correction
void StEpdEpManager::initEpdShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ2ShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mEpdQ2ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ2ShiftSin%dEastVz%d",iShift,iVz);
      p_mEpdQ2ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mEpdQ2ShiftCos%dWestVz%d",iShift,iVz);
      p_mEpdQ2ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ2ShiftSin%dWestVz%d",iShift,iVz);
      p_mEpdQ2ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mEpdQ3ShiftCos%dEastVz%d",iShift,iVz); // 3rd EP
      p_mEpdQ3ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ3ShiftSin%dEastVz%d",iShift,iVz);
      p_mEpdQ3ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mEpdQ3ShiftCos%dWestVz%d",iShift,iVz);
      p_mEpdQ3ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ3ShiftSin%dWestVz%d",iShift,iVz);
      p_mEpdQ3ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    }
  }
}

void StEpdEpManager::fillEpdShiftEast()
{
  TVector2 Q2Vector = getQ2VecReCenterEast();
  const double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi2Cos = TMath::Cos(2.0*((double)iShift+1.0)*Psi2);
    const double Psi2Sin = TMath::Sin(2.0*((double)iShift+1.0)*Psi2);
    p_mEpdQ2ShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Cos);
    p_mEpdQ2ShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Sin);
  }

  TVector2 Q3Vector = getQ3VecReCenterEast();
  const double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi3Cos = TMath::Cos(3.0*((double)iShift+1.0)*Psi3);
    const double Psi3Sin = TMath::Sin(3.0*((double)iShift+1.0)*Psi3);
    p_mEpdQ3ShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Cos);
    p_mEpdQ3ShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Sin);
  }
}

void StEpdEpManager::fillEpdShiftWest()
{
  TVector2 Q2Vector = getQ2VecReCenterWest();
  const double Psi2 = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi2Cos = TMath::Cos(2.0*((double)iShift+1.0)*Psi2);
    const double Psi2Sin = TMath::Sin(2.0*((double)iShift+1.0)*Psi2);
    p_mEpdQ2ShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Cos);
    p_mEpdQ2ShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi2Sin);
  }

  TVector2 Q3Vector = getQ3VecReCenterWest();
  const double Psi3 = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0;
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi3Cos = TMath::Cos(3.0*((double)iShift+1.0)*Psi3);
    const double Psi3Sin = TMath::Sin(3.0*((double)iShift+1.0)*Psi3);
    p_mEpdQ3ShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Cos);
    p_mEpdQ3ShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi3Sin);
  }
}

void StEpdEpManager::writeEpdShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mEpdQ2ShiftCosEast[iVz][iShift]->Write();
      p_mEpdQ2ShiftSinEast[iVz][iShift]->Write();
      p_mEpdQ2ShiftCosWest[iVz][iShift]->Write();
      p_mEpdQ2ShiftSinWest[iVz][iShift]->Write();

      p_mEpdQ3ShiftCosEast[iVz][iShift]->Write();
      p_mEpdQ3ShiftSinEast[iVz][iShift]->Write();
      p_mEpdQ3ShiftCosWest[iVz][iShift]->Write();
      p_mEpdQ3ShiftSinWest[iVz][iShift]->Write();
    }
  }
}

void StEpdEpManager::readEpdShift()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftParameter/file_%s_EpdShiftPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ2ShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mEpdQ2ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ2ShiftSin%dEastVz%d",iShift,iVz);
      p_mEpdQ2ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mEpdQ2ShiftCos%dWestVz%d",iShift,iVz);
      p_mEpdQ2ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ2ShiftSin%dWestVz%d",iShift,iVz);
      p_mEpdQ2ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mEpdQ3ShiftCos%dEastVz%d",iShift,iVz); // 3rd EP
      p_mEpdQ3ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ3ShiftSin%dEastVz%d",iShift,iVz);
      p_mEpdQ3ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mEpdQ3ShiftCos%dWestVz%d",iShift,iVz);
      p_mEpdQ3ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ3ShiftSin%dWestVz%d",iShift,iVz);
      p_mEpdQ3ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

double StEpdEpManager::getPsi2ShiftEast()
{
  TVector2 Q2Vector = getQ2VecReCenterEast();
  const double Psi2ReCenter = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2

  double deltaPsi2 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mEpdQ2ShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mEpdQ2ShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mEpdQ2ShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mEpdQ2ShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi2 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*Psi2ReCenter)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*Psi2ReCenter));
  }

  double Psi2ShiftRaw = Psi2ReCenter + deltaPsi2/2.0;
  double Psi2Shift    = transPsi2(Psi2ShiftRaw);

  return Psi2Shift;
}

double StEpdEpManager::getPsi2ShiftWest()
{
  TVector2 Q2Vector = getQ2VecReCenterWest();

  const double Psi2ReCenter = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0;
  double deltaPsi2 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mEpdQ2ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mEpdQ2ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mEpdQ2ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mEpdQ2ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi2 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(2.0*((double)iShift+1.0)*Psi2ReCenter)+meanCos*TMath::Sin(2.0*((double)iShift+1.0)*Psi2ReCenter));
  }

  double Psi2ShiftRaw = Psi2ReCenter + deltaPsi2/2.0;
  double Psi2Shift    = transPsi2(Psi2ShiftRaw);

  return Psi2Shift;
}

double StEpdEpManager::transPsi2(double Psi2)
{
  double Psi2Corr = Psi2;
  if(Psi2 >  TMath::Pi()/2.0) Psi2Corr = Psi2 - TMath::Pi();
  if(Psi2 < -TMath::Pi()/2.0) Psi2Corr = Psi2 + TMath::Pi();

  return Psi2Corr;
}

double StEpdEpManager::getPsi3ShiftEast()
{
  TVector2 Q3Vector = getQ3VecReCenterEast();
  const double Psi3ReCenter = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3

  double deltaPsi3 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mEpdQ3ShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mEpdQ3ShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mEpdQ3ShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mEpdQ3ShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi3 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(3.0*((double)iShift+1.0)*Psi3ReCenter)+meanCos*TMath::Sin(3.0*((double)iShift+1.0)*Psi3ReCenter));
  }

  double Psi3ShiftRaw = Psi3ReCenter + deltaPsi3/3.0;
  double Psi3Shift    = transPsi3(Psi3ShiftRaw);

  return Psi3Shift;
}

double StEpdEpManager::getPsi3ShiftWest()
{
  TVector2 Q3Vector = getQ3VecReCenterWest();

  const double Psi3ReCenter = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0;
  double deltaPsi3 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mEpdQ3ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mEpdQ3ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mEpdQ3ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mEpdQ3ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi3 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(3.0*((double)iShift+1.0)*Psi3ReCenter)+meanCos*TMath::Sin(3.0*((double)iShift+1.0)*Psi3ReCenter));
  }

  double Psi3ShiftRaw = Psi3ReCenter + deltaPsi3/3.0;
  double Psi3Shift    = transPsi3(Psi3ShiftRaw);

  return Psi3Shift;
}

double StEpdEpManager::transPsi3(double Psi3)
{
  double Psi3Corr = Psi3;
  if(Psi3 >  TMath::Pi()/3.0) Psi3Corr = Psi3 - TMath::TwoPi()/3.0;
  if(Psi3 < -TMath::Pi()/3.0) Psi3Corr = Psi3 + TMath::TwoPi()/3.0;

  return Psi3Corr;
}
//---------------------------------------------------------------------------------
// Sub Event Plane Resolution
void StEpdEpManager::initEpdResolution()
{
  p_mEpdSubEp2Res = new TProfile("p_mEpdSubEp2Res","p_mEpdSubEp2Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  p_mEpdSubEp3Res = new TProfile("p_mEpdSubEp3Res","p_mEpdSubEp3Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
}

void StEpdEpManager::fillEpdResolution(double Psi2East, double Psi2West, double Psi3East, double Psi3West)
{
  double res2Sub = TMath::Cos(2.0*(Psi2East-Psi2West));
  p_mEpdSubEp2Res->Fill((double)mCent9,res2Sub);

  double res3Sub = TMath::Cos(3.0*(Psi3East-Psi3West));
  p_mEpdSubEp3Res->Fill((double)mCent9,res3Sub);
}

void StEpdEpManager::writeEpdResolution()
{
  p_mEpdSubEp2Res->Write();
  p_mEpdSubEp3Res->Write();
}

void StEpdEpManager::readEpdResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_%s_EpdEpResolution.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mEpdSubEp2Res = (TProfile*)file_mResolution->Get("p_mEpdSubEp2Res");
  p_mEpdSubEp3Res = (TProfile*)file_mResolution->Get("p_mEpdSubEp3Res");
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    mEpdSubEp2ResVal[iCent] = 0.0;
    mEpdSubEp2ResErr[iCent] = 0.0;
    mEpdSubEp3ResVal[iCent] = 0.0;
    mEpdSubEp3ResErr[iCent] = 0.0;
  }

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    double valRes2Sub = -999.9;
    double errRes2Sub = 1.0;
    double valRes3Sub = -999.9;
    double errRes3Sub = 1.0;
    double valRes2Raw = p_mEpdSubEp2Res->GetBinContent(iCent+1);
    double errRes2Raw = p_mEpdSubEp2Res->GetBinError(iCent+1);
    double valRes3Raw = p_mEpdSubEp3Res->GetBinContent(iCent+1);
    double errRes3Raw = p_mEpdSubEp3Res->GetBinError(iCent+1);
    if(valRes2Raw > 0)
    {
      valRes2Sub = TMath::Sqrt(valRes2Raw);
      errRes2Sub = errRes2Raw/(2*valRes2Sub);
    }

    if(valRes3Raw > 0)
    {
      valRes3Sub = TMath::Sqrt(valRes3Raw);
      errRes3Sub = errRes3Raw/(2*valRes3Sub);
    }

    mEpdSubEp2ResVal[iCent] = valRes2Sub;
    mEpdSubEp2ResErr[iCent] = errRes2Sub;
    mEpdSubEp3ResVal[iCent] = valRes3Sub;
    mEpdSubEp3ResErr[iCent] = errRes3Sub;
  }
  file_mResolution->Close();
}

double StEpdEpManager::getEpdSubEp2ResVal(int cent9)
{
  return mEpdSubEp2ResVal[cent9];
}

double StEpdEpManager::getEpdSubEp2ResErr(int cent9)
{
  return mEpdSubEp2ResErr[cent9];
}

double StEpdEpManager::getEpdSubEp3ResVal(int cent9)
{
  return mEpdSubEp3ResVal[cent9];
}

double StEpdEpManager::getEpdSubEp3ResErr(int cent9)
{
  return mEpdSubEp3ResErr[cent9];
}
//---------------------------------------------------------------------------------
// QVector
TVector2 StEpdEpManager::getQ2VecRawEast()
{
  return v_mQ2RawEast;
}

TVector2 StEpdEpManager::getQ2VecRawWest()
{
  return v_mQ2RawWest;
}

TVector2 StEpdEpManager::getQ2VecReCenterEast()
{
  return v_mQ2ReCenterEast;
}

TVector2 StEpdEpManager::getQ2VecReCenterWest()
{
  return v_mQ2ReCenterWest;
}

double StEpdEpManager::getPsi2RawEast()
{
  TVector2 Q2Vector = getQ2VecRawEast();
  double Psi2       = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2Raw    = transPsi2(Psi2);

  return Psi2Raw;
}

double StEpdEpManager::getPsi2RawWest()
{
  TVector2 Q2Vector = getQ2VecRawWest();
  double Psi2       = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2Raw    = transPsi2(Psi2);

  return Psi2Raw;
}

double StEpdEpManager::getPsi2ReCenterEast()
{
  TVector2 Q2Vector   = getQ2VecReCenterEast();
  double Psi2         = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2ReCenter = transPsi2(Psi2);

  return Psi2ReCenter;
}

double StEpdEpManager::getPsi2ReCenterWest()
{
  TVector2 Q2Vector   = getQ2VecReCenterWest();
  double Psi2         = TMath::ATan2(Q2Vector.Y(),Q2Vector.X())/2.0; // -pi/2 to pi/2
  double Psi2ReCenter = transPsi2(Psi2);

  return Psi2ReCenter;
}

TVector2 StEpdEpManager::getQ3VecRawEast()
{
  return v_mQ3RawEast;
}

TVector2 StEpdEpManager::getQ3VecRawWest()
{
  return v_mQ3RawWest;
}

TVector2 StEpdEpManager::getQ3VecReCenterEast()
{
  return v_mQ3ReCenterEast;
}

TVector2 StEpdEpManager::getQ3VecReCenterWest()
{
  return v_mQ3ReCenterWest;
}

double StEpdEpManager::getPsi3RawEast()
{
  TVector2 Q3Vector = getQ3VecRawEast();
  double Psi3       = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
  double Psi3Raw    = transPsi3(Psi3);

  return Psi3Raw;
}

double StEpdEpManager::getPsi3RawWest()
{
  TVector2 Q3Vector = getQ3VecRawWest();
  double Psi3       = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
  double Psi3Raw    = transPsi3(Psi3);

  return Psi3Raw;
}

double StEpdEpManager::getPsi3ReCenterEast()
{
  TVector2 Q3Vector   = getQ3VecReCenterEast();
  double Psi3         = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
  double Psi3ReCenter = transPsi3(Psi3);

  return Psi3ReCenter;
}

double StEpdEpManager::getPsi3ReCenterWest()
{
  TVector2 Q3Vector   = getQ3VecReCenterWest();
  double Psi3         = TMath::ATan2(Q3Vector.Y(),Q3Vector.X())/3.0; // -pi/3 to pi/3
  double Psi3ReCenter = transPsi3(Psi3);

  return Psi3ReCenter;
}

int StEpdEpManager::getNumTrkRawEast()
{
  return mQCouRawEast;
}

int StEpdEpManager::getNumTrkRawWest()
{
  return mQCouRawWest;
}

int StEpdEpManager::getNumTrkReCenterEast()
{
  return mQCouReCenterEast;
}

int StEpdEpManager::getNumTrkReCenterWest()
{
  return mQCouReCenterWest;
}
//---------------------------------------------------------------------------------
// raw EP
void StEpdEpManager::initEpdSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp2RawEastCent%d",iCent); // 2nd EP
    h_mEpdEp2RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp2RawWestCent%d",iCent);
    h_mEpdEp2RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp2RawCorrCent%d",iCent);
    h_mEpdEp2RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());

    histName = Form("h_mEpdEp3RawEastCent%d",iCent); // 3rd EP
    h_mEpdEp3RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp3RawWestCent%d",iCent);
    h_mEpdEp3RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp3RawCorrCent%d",iCent);
    h_mEpdEp3RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpRaw(double Psi2East, double Psi2West, double Psi3East, double Psi3West)
{
  h_mEpdEp2RawEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mEpdEp2RawWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mEpdEp2RawCorr[mCent9]->Fill(Psi2East,Psi2West);

  h_mEpdEp3RawEast[mCent9]->Fill(mRunIndex,Psi3East);
  h_mEpdEp3RawWest[mCent9]->Fill(mRunIndex,Psi3West);
  h_mEpdEp3RawCorr[mCent9]->Fill(Psi3East,Psi3West);
}

void StEpdEpManager::writeEpdSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp2RawEast[iCent]->Write();
    h_mEpdEp2RawWest[iCent]->Write();
    h_mEpdEp2RawCorr[iCent]->Write();

    h_mEpdEp3RawEast[iCent]->Write();
    h_mEpdEp3RawWest[iCent]->Write();
    h_mEpdEp3RawCorr[iCent]->Write();
  }
}

// recenter EP
void StEpdEpManager::initEpdSubEpReCenter()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp2ReCenterEastCent%d",iCent); // 2nd EP
    h_mEpdEp2ReCenterEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp2ReCenterWestCent%d",iCent);
    h_mEpdEp2ReCenterWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp2ReCenterCorrCent%d",iCent);
    h_mEpdEp2ReCenterCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());

    histName = Form("h_mEpdEp3ReCenterEastCent%d",iCent); // 3rd EP
    h_mEpdEp3ReCenterEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp3ReCenterWestCent%d",iCent);
    h_mEpdEp3ReCenterWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp3ReCenterCorrCent%d",iCent);
    h_mEpdEp3ReCenterCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpReCenter(double Psi2East, double Psi2West, double Psi3East, double Psi3West)
{
  h_mEpdEp2ReCenterEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mEpdEp2ReCenterWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mEpdEp2ReCenterCorr[mCent9]->Fill(Psi2East,Psi2West);

  h_mEpdEp3ReCenterEast[mCent9]->Fill(mRunIndex,Psi3East);
  h_mEpdEp3ReCenterWest[mCent9]->Fill(mRunIndex,Psi3West);
  h_mEpdEp3ReCenterCorr[mCent9]->Fill(Psi3East,Psi3West);
}

void StEpdEpManager::writeEpdSubEpReCenter()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp2ReCenterEast[iCent]->Write();
    h_mEpdEp2ReCenterWest[iCent]->Write();
    h_mEpdEp2ReCenterCorr[iCent]->Write();

    h_mEpdEp3ReCenterEast[iCent]->Write();
    h_mEpdEp3ReCenterWest[iCent]->Write();
    h_mEpdEp3ReCenterCorr[iCent]->Write();
  }
}

// shift EP
void StEpdEpManager::initEpdSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp2ShiftEastCent%d",iCent); // 2nd EP
    h_mEpdEp2ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp2ShiftWestCent%d",iCent);
    h_mEpdEp2ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp2ShiftCorrCent%d",iCent);
    h_mEpdEp2ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());

    histName = Form("h_mEpdEp3ShiftEastCent%d",iCent); // 3rd EP
    h_mEpdEp3ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp3ShiftWestCent%d",iCent);
    h_mEpdEp3ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mEpdEp3ShiftCorrCent%d",iCent);
    h_mEpdEp3ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpShift(double Psi2East, double Psi2West, double Psi3East, double Psi3West)
{
  h_mEpdEp2ShiftEast[mCent9]->Fill(mRunIndex,Psi2East);
  h_mEpdEp2ShiftWest[mCent9]->Fill(mRunIndex,Psi2West);
  h_mEpdEp2ShiftCorr[mCent9]->Fill(Psi2East,Psi2West);

  h_mEpdEp3ShiftEast[mCent9]->Fill(mRunIndex,Psi3East);
  h_mEpdEp3ShiftWest[mCent9]->Fill(mRunIndex,Psi3West);
  h_mEpdEp3ShiftCorr[mCent9]->Fill(Psi3East,Psi3West);
}

void StEpdEpManager::writeEpdSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp2ShiftEast[iCent]->Write();
    h_mEpdEp2ShiftWest[iCent]->Write();
    h_mEpdEp2ShiftCorr[iCent]->Write();

    h_mEpdEp3ShiftEast[iCent]->Write();
    h_mEpdEp3ShiftWest[iCent]->Write();
    h_mEpdEp3ShiftCorr[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
#endif
