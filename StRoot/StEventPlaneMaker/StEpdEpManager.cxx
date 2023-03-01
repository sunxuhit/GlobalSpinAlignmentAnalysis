#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TF1.h"

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

double funcEpdEpResFull(double *x_val, double *par)
{
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  double resFull = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return resFull;
}

ClassImp(StEpdEpManager)

//---------------------------------------------------------------------------------
StEpdEpManager::StEpdEpManager(int beamType) : mType(beamType)
{
  // mEnergy = energy;
  clearEpdEpManager();
  mEpdGeom = new StEpdGeom();
  mUsePhiWgt = false;
  mUseEtaWgt = false;
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

  mQ1WgtSideRawEast = 0.0;
  mQ1WgtSideRawWest = 0.0;
  v_mQ1SideRawEast.Set(0.0,0.0);
  v_mQ1SideRawWest.Set(0.0,0.0);

  mQ1WgtSideWgtEast = 0.0;
  mQ1WgtSideWgtWest = 0.0;
  v_mQ1SideWgtEast.Set(0.0,0.0);
  v_mQ1SideWgtWest.Set(0.0,0.0);

  mQ1WgtSideReCtrEast = 0.0;
  mQ1WgtSideReCtrWest = 0.0;
  v_mQ1SideReCtrEast.Set(0.0,0.0);
  v_mQ1SideReCtrWest.Set(0.0,0.0);
  /*
  for(int iRing = 0; iRing < mNumRings; ++iRing)
  {
    mQ1WgtRingRawEast[iRing] = 0.0
    mQ1WgtRingRawWest[iRing] = 0.0
    v_mQ1RingRawEast[iRing].Set(0.0,0.0);
    v_mQ1RingRawWest[iRing].Set(0.0,0.0);

    mQ1WgtRingWgtEast[iRing] = 0.0
    mQ1WgtRingWgtWest[iRing] = 0.0
    v_mQ1RingWgtEast[iRing].Set(0.0,0.0);
    v_mQ1RingWgtWest[iRing].Set(0.0,0.0);

    mQ1WgtRingReCtrEast[iRing] = 0.0
    mQ1WgtRingReCtrWest[iRing] = 0.0
    v_mQ1RingReCtrEast[iRing].Set(0.0,0.0);
    v_mQ1RingReCtrWest[iRing].Set(0.0,0.0);
  }
  */
}

void StEpdEpManager::initEpdEpManager(int cent9, int runIndex, int vzBin)
{
  mCent9    = cent9;
  mRunIndex = runIndex;
  mVzBin    = vzBin;
}
//---------------------------------------------------------------------------------
// Utilities
TVector3 StEpdEpManager::getEpdRanVec(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  TVector3 EpdPoint  = mEpdGeom->RandomPointOnTile(picoEpdHit->id()); // get a random position within the tile
  TVector3 EpdVector = EpdPoint - primVtx;

  return EpdVector;
}

TVector2 StEpdEpManager::calq1Vector(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  TVector2 q1Vector(0.0,0.0);
  TVector3 EpdVector = getEpdRanVec(picoEpdHit, primVtx);
  const double phi   = EpdVector.Phi(); // -pi to pi
  const double q1x   = TMath::Cos(1.0*phi);
  const double q1y   = TMath::Sin(1.0*phi);
  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

double StEpdEpManager::getTileWeight(StPicoEpdHit *picoEpdHit)
{
  const double nMip = picoEpdHit->nMIP();
  double tileWgt = (nMip < anaUtils::mMipEpdEpMax[mType]) ? nMip : anaUtils::mMipEpdEpMax[mType];
  if(nMip < anaUtils::mMipEpdEpMin[mType]) tileWgt = 0.0;

  return tileWgt;
}

double StEpdEpManager::getPhiWeight(StPicoEpdHit *picoEpdHit)
{
  double phiWgt = 1.0;
  if(mUsePhiWgt)
  {
    const int tileId      = picoEpdHit->id();
    const int secId       = picoEpdHit->position()-1; // convert to 0-11
    const int tileIdOnSec = picoEpdHit->tile()-1;     // conver to 0-30
    if(tileId < 0) // East EPD
    {
      int binWgt = h_mEpdPhiWgtEast[mCent9]->FindBin((double)secId,(double)tileIdOnSec);
      phiWgt     = 1.0/h_mEpdPhiWgtEast[mCent9]->GetBinContent(binWgt);
    }
    if(tileId > 0) // West EPD
    {
      int binWgt = h_mEpdPhiWgtWest[mCent9]->FindBin((double)secId,(double)tileIdOnSec);
      phiWgt     = 1.0/h_mEpdPhiWgtWest[mCent9]->GetBinContent(binWgt);
    }
  }

  return phiWgt;
}

double StEpdEpManager::getEtaWeight(StPicoEpdHit *picoEpdHit)
{
  double etaWgt = 1.0; // TODO: add real eta weight

  return etaWgt;
}

void StEpdEpManager::addHitRawEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId < mNumRingsUsed)
  {
    v_mQ1SideRawEast  += tileWgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtSideRawEast += tileWgt;
  }

  // v_mQ1RingRawEast[ringId]  += tileWgt*calq1Vector(picoEpdHit,primVtx);
  // mQ1WgtRingRawEast[ringId] += tileWgt;
}

void StEpdEpManager::addHitRawWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId < mNumRingsUsed)
  {
    v_mQ1SideRawWest  += tileWgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtSideRawWest += tileWgt;
  }

  // v_mQ1RingRawWest[ringId]  += tileWgt*calq1Vector(picoEpdHit,primVtx);
  // mQ1WgtRingRawWest[ringId] += tileWgt;
}

void StEpdEpManager::addHitWgtEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId < mNumRingsUsed)
  {
    v_mQ1SideWgtEast  += wgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtSideWgtEast += wgt;
  }

  // v_mQ1RingWgtEast[ringId]  += wgt*calq1Vector(picoEpdHit,primVtx);
  // mQ1WgtRingWgtEast[ringId] += wgt;
}

void StEpdEpManager::addHitWgtWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId < mNumRingsUsed)
  {
    v_mQ1SideWgtWest  += wgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtSideWgtWest += wgt;
  }

  // v_mQ1RingWgtWest[ringId]  += wgt*calq1Vector(picoEpdHit,primVtx);
  // mQ1WgtRingWgtWest[ringId] += wgt;
}

void StEpdEpManager::addHitReCtrEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId < mNumRingsUsed)
  {
    v_mQ1SideReCtrEast  += wgt*(calq1Vector(picoEpdHit,primVtx) - getq1VecCtrEast());
    mQ1WgtSideReCtrEast += wgt;
  }
}

void StEpdEpManager::addHitReCtrWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId < mNumRingsUsed)
  {
    v_mQ1SideReCtrWest  += wgt*(calq1Vector(picoEpdHit,primVtx) - getq1VecCtrWest());
    mQ1WgtSideReCtrWest += wgt;
  }
}
//---------------------------------------------------------------------------------
// phi Weight Correction
void StEpdEpManager::initEpdPhiWgt()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdPhiWgtEastCent%d",iCent);
    h_mEpdPhiWgtEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),mNumSectors,-0.5,(double)mNumSectors-0.5,mNumTiles,-0.5,(double)mNumTiles-0.5);
    histName = Form("h_mEpdPhiAveEastCent%d",iCent);
    h_mEpdPhiAveEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),mNumSectors,-0.5,(double)mNumSectors-0.5,mNumTiles,-0.5,(double)mNumTiles-0.5);

    histName = Form("h_mEpdPhiWgtWestCent%d",iCent);
    h_mEpdPhiWgtWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),mNumSectors,-0.5,(double)mNumSectors-0.5,mNumTiles,-0.5,(double)mNumTiles-0.5);
    histName = Form("h_mEpdPhiAveWestCent%d",iCent);
    h_mEpdPhiAveWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),mNumSectors,-0.5,(double)mNumSectors-0.5,mNumTiles,-0.5,(double)mNumTiles-0.5);
  }
}

void StEpdEpManager::fillEpdPhiWgtEast(StPicoEpdHit* picoEpdHit)
{
  const int secId       = picoEpdHit->position()-1; // convert to 0-11
  const int tileIdOnSec = picoEpdHit->tile()-1;     // conver to 0-30
  const int ringId      = picoEpdHit->row()-1;      // conver to 0-15
  const double tileWgt  = getTileWeight(picoEpdHit);
  h_mEpdPhiWgtEast[mCent9]->Fill((double)secId,(double)tileIdOnSec,tileWgt);

  for(int iSec = 0; iSec < mNumSectors; ++iSec)
  {
    if(tileIdOnSec == 0)
    {
      h_mEpdPhiAveEast[mCent9]->Fill((double)iSec,0.0,tileWgt/12.0); // average over 12 tiles in the most inner ring
    }
    else
    {
      for(int iTile = 2*ringId-1; iTile < 2*ringId+1; ++iTile)
      {
	h_mEpdPhiAveEast[mCent9]->Fill((double)iSec,(double)iTile,tileWgt/24.0); // average over 24 tiles
      }
    }
  }
}

void StEpdEpManager::fillEpdPhiWgtWest(StPicoEpdHit* picoEpdHit)
{
  const int secId       = picoEpdHit->position()-1; // convert to 0-11
  const int tileIdOnSec = picoEpdHit->tile()-1;     // conver to 0-30
  const int ringId      = picoEpdHit->row()-1;      // conver to 0-15
  const double tileWgt  = getTileWeight(picoEpdHit);
  h_mEpdPhiWgtWest[mCent9]->Fill((double)secId,(double)tileIdOnSec,tileWgt);

  for(int iSec = 0; iSec < mNumSectors; ++iSec)
  {
    if(tileIdOnSec == 0)
    {
      h_mEpdPhiAveWest[mCent9]->Fill((double)iSec,0.0,tileWgt/12.0); // average over 12 tiles in the most inner ring
    }
    else
    {
      for(int iTile = 2*ringId-1; iTile < 2*ringId+1; ++iTile)
      {
	h_mEpdPhiAveWest[mCent9]->Fill((double)iSec,(double)iTile,tileWgt/24.0); // average over 24 tiles
      }
    }
  }
}

void StEpdEpManager::writeEpdPhiWgt()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdPhiWgtEast[iCent]->Write();
    h_mEpdPhiAveEast[iCent]->Write();
    h_mEpdPhiWgtWest[iCent]->Write();
    h_mEpdPhiAveWest[iCent]->Write();
  }
}

void StEpdEpManager::readEpdPhiWgt()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/GainCorrPar/file_EpdPhiWgtPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mPhiWgtPar = TFile::Open(inputFile.c_str());
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdPhiWgtEastCent%d",iCent);
    h_mEpdPhiWgtEast[iCent] = (TH2F*)file_mPhiWgtPar->Get(histName.c_str());

    histName = Form("h_mEpdPhiWgtWestCent%d",iCent);
    h_mEpdPhiWgtWest[iCent] = (TH2F*)file_mPhiWgtPar->Get(histName.c_str());
  }
  mUsePhiWgt = true;
}
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// ReCtrPar Correction
void StEpdEpManager::initEpdReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mEpdQ1ReCtrXEastVz%d",iVz); // 1st EP
    p_mEpdQ1ReCtrXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mEpdQ1ReCtrYEastVz%d",iVz);
    p_mEpdQ1ReCtrYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mEpdQ1ReCtrXWestVz%d",iVz);
    p_mEpdQ1ReCtrXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mEpdQ1ReCtrYWestVz%d",iVz);
    p_mEpdQ1ReCtrYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  }
}

void StEpdEpManager::fillEpdReCtrEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;

  const TVector2 q1Vector = calq1Vector(picoEpdHit,primVtx);
  const double q1x = q1Vector.X();
  const double q1y = q1Vector.Y();

  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId < mNumRingsUsed)
  {
    p_mEpdQ1ReCtrXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1x,wgt);
    p_mEpdQ1ReCtrYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1y,wgt);
  }
}

void StEpdEpManager::fillEpdReCtrWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;

  const TVector2 q1Vector = calq1Vector(picoEpdHit,primVtx);
  const double q1x = q1Vector.X();
  const double q1y = q1Vector.Y();

  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId < mNumRingsUsed)
  {
    p_mEpdQ1ReCtrXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1x,wgt);
    p_mEpdQ1ReCtrYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1y,wgt);
  }
}

void StEpdEpManager::writeEpdReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mEpdQ1ReCtrXEast[iVz]->Write();
    p_mEpdQ1ReCtrYEast[iVz]->Write();
    p_mEpdQ1ReCtrXWest[iVz]->Write();
    p_mEpdQ1ReCtrYWest[iVz]->Write();
  }
}

void StEpdEpManager::readEpdReCtr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_EpdReCenterPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mReCtrPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mEpdQ1ReCtrXEastVz%d",iVz); // 1st EP
    p_mEpdQ1ReCtrXEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mEpdQ1ReCtrYEastVz%d",iVz);
    p_mEpdQ1ReCtrYEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mEpdQ1ReCtrXWestVz%d",iVz);
    p_mEpdQ1ReCtrXWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mEpdQ1ReCtrYWestVz%d",iVz);
    p_mEpdQ1ReCtrYWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
  }
}

TVector2 StEpdEpManager::getq1VecCtrEast()
{
  const int binX   = p_mEpdQ1ReCtrXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1x = p_mEpdQ1ReCtrXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mEpdQ1ReCtrYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1y = p_mEpdQ1ReCtrYEast[mVzBin]->GetBinContent(binY);

  TVector2 q1Vector(0.0,0.0);
  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

TVector2 StEpdEpManager::getq1VecCtrWest()
{
  const int binX   = p_mEpdQ1ReCtrXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1x = p_mEpdQ1ReCtrXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mEpdQ1ReCtrYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1y = p_mEpdQ1ReCtrYWest[mVzBin]->GetBinContent(binY);

  TVector2 q1Vector(0.0,0.0);
  q1Vector.Set(q1x,q1y);

  return q1Vector;
}
//---------------------------------------------------------------------------------
// Shift Correction
void StEpdEpManager::initEpdShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1ShiftCos%dEastVz%d",iShift,iVz); // 1st EP
      p_mEpdQ1ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1ShiftSin%dEastVz%d",iShift,iVz);
      p_mEpdQ1ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mEpdQ1ShiftCos%dWestVz%d",iShift,iVz);
      p_mEpdQ1ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1ShiftSin%dWestVz%d",iShift,iVz);
      p_mEpdQ1ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    }
  }
}

void StEpdEpManager::fillEpdShiftEast()
{
  TVector2 Q1Vector = getQ1VecReCtrEast();
  const double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X());
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
    const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
    p_mEpdQ1ShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
    p_mEpdQ1ShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
  }
}

void StEpdEpManager::fillEpdShiftWest()
{
  TVector2 Q1Vector = getQ1VecReCtrWest();
  const double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X());
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
    const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
    p_mEpdQ1ShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
    p_mEpdQ1ShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
  }
}

void StEpdEpManager::writeEpdShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mEpdQ1ShiftCosEast[iVz][iShift]->Write();
      p_mEpdQ1ShiftSinEast[iVz][iShift]->Write();
      p_mEpdQ1ShiftCosWest[iVz][iShift]->Write();
      p_mEpdQ1ShiftSinWest[iVz][iShift]->Write();
    }
  }
}

void StEpdEpManager::readEpdShift()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftPar/file_EpdShiftPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1ShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mEpdQ1ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1ShiftSin%dEastVz%d",iShift,iVz);
      p_mEpdQ1ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mEpdQ1ShiftCos%dWestVz%d",iShift,iVz);
      p_mEpdQ1ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1ShiftSin%dWestVz%d",iShift,iVz);
      p_mEpdQ1ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

double StEpdEpManager::getPsi1ShiftEast()
{
  TVector2 Q1Vector = getQ1VecReCtrEast();
  const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi

  double deltaPsi1 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mEpdQ1ShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mEpdQ1ShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mEpdQ1ShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mEpdQ1ShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1ReCtr));
  }

  double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
  double Psi1Shift    = transPsi1(Psi1ShiftRaw);

  return Psi1Shift;
}

double StEpdEpManager::getPsi1ShiftWest()
{
  TVector2 Q1Vector = getQ1VecReCtrWest();
  const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi

  double deltaPsi1 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mEpdQ1ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mEpdQ1ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mEpdQ1ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mEpdQ1ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1ReCtr));
  }

  double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
  double Psi1Shift    = transPsi1(Psi1ShiftRaw);

  return Psi1Shift;
}

double StEpdEpManager::getPsi1ShiftFull()
{
  TVector2 Q1Vector = getQ1VecShiftFull();
  double Psi1ShiftRaw = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Shift    = transPsi1(Psi1ShiftRaw);

  return Psi1Shift;
}

double StEpdEpManager::transPsi1(double Psi1)
{
  double Psi1Corr = Psi1;
  if(Psi1 >  1.0*TMath::Pi()) Psi1Corr = Psi1 - TMath::TwoPi();
  if(Psi1 < -1.0*TMath::Pi()) Psi1Corr = Psi1 + TMath::TwoPi();

  return Psi1Corr;
}
//---------------------------------------------------------------------------------
// Shift Correction Full EP
void StEpdEpManager::initEpdShiftFull()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1ShiftCos%dFullVz%d",iShift,iVz);
      p_mEpdQ1ShiftCosFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1ShiftSin%dFullVz%d",iShift,iVz);
      p_mEpdQ1ShiftSinFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    }
  }
}

void StEpdEpManager::fillEpdShiftFull()
{
  TVector2 Q1Vector = getQ1VecShiftFull();
  const double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X());
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
    const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
    p_mEpdQ1ShiftCosFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
    p_mEpdQ1ShiftSinFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
  }
}

void StEpdEpManager::writeEpdShiftFull()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mEpdQ1ShiftCosFull[iVz][iShift]->Write();
      p_mEpdQ1ShiftSinFull[iVz][iShift]->Write();
    }
  }
}

void StEpdEpManager::readEpdShiftFull()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftPar/file_EpdShiftParFull_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(inputFile.c_str());

  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < 20; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1ShiftCos%dFullVz%d",iShift,iVz);
      p_mEpdQ1ShiftCosFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1ShiftSin%dFullVz%d",iShift,iVz);
      p_mEpdQ1ShiftSinFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

double StEpdEpManager::getPsi1ShiftFullCorr()
{
  TVector2 Q1Vector = getQ1VecShiftFull();
  const double Psi1Shift = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi

  double deltaPsi1 = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mEpdQ1ShiftCosFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mEpdQ1ShiftCosFull[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mEpdQ1ShiftSinFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mEpdQ1ShiftSinFull[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1Shift)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1Shift));
  }

  double Psi1ShiftCorrRaw = Psi1Shift + deltaPsi1;
  double Psi1ShiftCorr    = transPsi1(Psi1ShiftCorrRaw);

  return Psi1ShiftCorr;
}

TVector2 StEpdEpManager::getQ1VecShiftFullCorr()
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1FullCorr = getPsi1ShiftFullCorr();
  const double Q1x = TMath::Cos(1.0*Psi1FullCorr);
  const double Q1y = TMath::Sin(1.0*Psi1FullCorr);
  Q1Vector.Set(Q1x,Q1y);

  return Q1Vector;
}

//---------------------------------------------------------------------------------
// Sub Event Plane Resolution
void StEpdEpManager::initEpdResolution()
{
  p_mEpdSubEp1Res = new TProfile("p_mEpdSubEp1Res","p_mEpdSubEp1Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
}

void StEpdEpManager::fillEpdResolution(double Psi1East, double Psi1West)
{
  double res1Sub = TMath::Cos(1.0*(Psi1West-Psi1East));
  p_mEpdSubEp1Res->Fill((double)mCent9,res1Sub);
}

void StEpdEpManager::writeEpdResolution()
{
  p_mEpdSubEp1Res->Write();
}

void StEpdEpManager::readEpdResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_EpdEpResolution_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mEpdSubEp1Res = (TProfile*)file_mResolution->Get("p_mEpdSubEp1Res");
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    mEpdSubEp1ResVal[iCent]  = 0.0;
    mEpdSubEp1ResErr[iCent]  = 0.0;
    mEpdFullEp1ResVal[iCent] = 0.0;
    mEpdFullEp1ResErr[iCent] = 0.0;
  }

  TF1 *f_EpdEpResFull = new TF1("f_EpdEpResFull",funcEpdEpResFull,0,10,0);
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    double valRes1Sub  = -999.9;
    double errRes1Sub  = 1.0;
    double valRes1Full = -999.9;
    double errRes1Full = 1.0;
    double valRes1Raw  = p_mEpdSubEp1Res->GetBinContent(iCent+1);
    double errRes1Raw  = p_mEpdSubEp1Res->GetBinError(iCent+1);
    if(valRes1Raw > 0)
    {
      valRes1Sub = TMath::Sqrt(valRes1Raw);
      errRes1Sub = errRes1Raw/(2.0*valRes1Sub);

      // calculate full event plane resolution
      double chiSub  = f_EpdEpResFull->GetX(valRes1Sub);
      double chiFull = chiSub*TMath::Sqrt(2.0);
      valRes1Full    = f_EpdEpResFull->Eval(chiFull);
      // error propagation
      double errChiSub = errRes1Sub/f_EpdEpResFull->Derivative(chiSub);
      errRes1Full = f_EpdEpResFull->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
    }

    mEpdSubEp1ResVal[iCent]  = valRes1Sub;
    mEpdSubEp1ResErr[iCent]  = errRes1Sub;
    mEpdFullEp1ResVal[iCent] = valRes1Full;
    mEpdFullEp1ResErr[iCent] = errRes1Full;
  }
  file_mResolution->Close();
}

double StEpdEpManager::getEpdSubEp1ResVal(int cent9)
{
  return mEpdSubEp1ResVal[cent9];
}

double StEpdEpManager::getEpdSubEp1ResErr(int cent9)
{
  return mEpdSubEp1ResErr[cent9];
}

double StEpdEpManager::getEpdFullEp1ResVal(int cent9)
{
  return mEpdFullEp1ResVal[cent9];
}

double StEpdEpManager::getEpdFullEp1ResErr(int cent9)
{
  return mEpdFullEp1ResErr[cent9];
}
//---------------------------------------------------------------------------------
// Charged Hadron Directed Flow
void StEpdEpManager::initEpdSubEpFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    std::string proName = Form("p_mEpdSubEpDFlowCent%d",i_cent);
    p_mEpdSubEpDFlow[i_cent] = new TProfile(proName.c_str(),proName.c_str(),40,-6.0,6.0);
  }
}

void StEpdEpManager::fillEpdSubEpDFlow(double eta, double v1, double reweight)
{
  p_mEpdSubEpDFlow[mCent9]->Fill(eta, v1, reweight);
}

void StEpdEpManager::writeEpdSubEpFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    p_mEpdSubEpDFlow[i_cent]->Write();
  }
}
//---------------------------------------------------------------------------------
// Q1Vector
TVector2 StEpdEpManager::getQ1VecRawEast()
{
  TVector2 Q1Vector(0.0, 0.0);
  const double Q1x = -1.0*v_mQ1SideRawEast.X()/mQ1WgtSideRawEast; // flip the sign for Q1VecEast
  const double Q1y = -1.0*v_mQ1SideRawEast.Y()/mQ1WgtSideRawEast;
  Q1Vector.Set(Q1x, Q1y);

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecRawWest()
{
  TVector2 Q1Vector(0.0, 0.0);
  const double Q1x = v_mQ1SideRawWest.X()/mQ1WgtSideRawWest;
  const double Q1y = v_mQ1SideRawWest.Y()/mQ1WgtSideRawWest;
  Q1Vector.Set(Q1x, Q1y);

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecRawFull()
{
  TVector2 Q1VecEast = getQ1VecRawEast();
  TVector2 Q1VecWest = getQ1VecRawWest();
  TVector2 Q1VecFull = Q1VecWest + Q1VecEast;

  return Q1VecFull;
}

TVector2 StEpdEpManager::getQ1VecWgtEast()
{
  TVector2 Q1Vector(0.0, 0.0);
  const double Q1x = -1.0*v_mQ1SideWgtEast.X()/mQ1WgtSideWgtEast; // flip the sign for Q1VecEast
  const double Q1y = -1.0*v_mQ1SideWgtEast.Y()/mQ1WgtSideWgtEast;
  Q1Vector.Set(Q1x, Q1y);

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecWgtWest()
{
  TVector2 Q1Vector(0.0, 0.0);
  const double Q1x = v_mQ1SideWgtWest.X()/mQ1WgtSideWgtWest;
  const double Q1y = v_mQ1SideWgtWest.Y()/mQ1WgtSideWgtWest;
  Q1Vector.Set(Q1x, Q1y);

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecWgtFull()
{
  TVector2 Q1VecEast = getQ1VecWgtEast();
  TVector2 Q1VecWest = getQ1VecWgtWest();
  TVector2 Q1VecFull = Q1VecWest + Q1VecEast;

  return Q1VecFull;
}

TVector2 StEpdEpManager::getQ1VecReCtrEast()
{
  TVector2 Q1Vector(0.0, 0.0);
  const double Q1x = -1.0*v_mQ1SideReCtrEast.X()/mQ1WgtSideReCtrEast; // flip the sign for Q1VecEast
  const double Q1y = -1.0*v_mQ1SideReCtrEast.Y()/mQ1WgtSideReCtrEast;
  Q1Vector.Set(Q1x, Q1y);

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecReCtrWest()
{
  TVector2 Q1Vector(0.0, 0.0);
  const double Q1x = v_mQ1SideReCtrWest.X()/mQ1WgtSideReCtrWest;
  const double Q1y = v_mQ1SideReCtrWest.Y()/mQ1WgtSideReCtrWest;
  Q1Vector.Set(Q1x, Q1y);

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecReCtrFull()
{
  TVector2 Q1VecEast = getQ1VecReCtrEast();
  TVector2 Q1VecWest = getQ1VecReCtrWest();
  TVector2 Q1VecFull = Q1VecWest + Q1VecEast;

  return Q1VecFull;
}

TVector2 StEpdEpManager::getQ1VecShiftEast()
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1East = getPsi1ShiftEast();
  const double Q1x = TMath::Cos(1.0*Psi1East);
  const double Q1y = TMath::Sin(1.0*Psi1East);
  Q1Vector.Set(Q1x,Q1y);

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecShiftWest()
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1West = getPsi1ShiftWest();
  const double Q1x = TMath::Cos(1.0*Psi1West);
  const double Q1y = TMath::Sin(1.0*Psi1West);
  Q1Vector.Set(Q1x,Q1y);

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecShiftFull()
{
  TVector2 Q1VecEast = getQ1VecShiftEast();
  TVector2 Q1VecWest = getQ1VecShiftWest();
  TVector2 Q1VecFull = Q1VecWest + Q1VecEast;

  return Q1VecFull;
}

double StEpdEpManager::getPsi1RawEast()
{
  TVector2 Q1Vector = getQ1VecRawEast();
  double Psi1       = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Raw    = transPsi1(Psi1);

  return Psi1Raw;
}

double StEpdEpManager::getPsi1RawWest()
{
  TVector2 Q1Vector = getQ1VecRawWest();
  double Psi1       = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Raw    = transPsi1(Psi1);

  return Psi1Raw;
}

double StEpdEpManager::getPsi1RawFull()
{
  TVector2 Q1Vector = getQ1VecRawFull();
  double Psi1       = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Raw    = transPsi1(Psi1);

  return Psi1Raw;
}

double StEpdEpManager::getPsi1WgtEast()
{
  TVector2 Q1Vector = getQ1VecWgtEast();
  double Psi1       = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Wgt    = transPsi1(Psi1);

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1WgtWest()
{
  TVector2 Q1Vector = getQ1VecWgtWest();
  double Psi1       = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Wgt    = transPsi1(Psi1);

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1WgtFull()
{
  TVector2 Q1Vector = getQ1VecWgtFull();
  double Psi1       = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Wgt    = transPsi1(Psi1);

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1ReCtrEast()
{
  TVector2 Q1Vector = getQ1VecReCtrEast();
  double Psi1       = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Wgt    = transPsi1(Psi1);

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1ReCtrWest()
{
  TVector2 Q1Vector = getQ1VecReCtrWest();
  double Psi1       = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Wgt    = transPsi1(Psi1);

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1ReCtrFull()
{
  TVector2 Q1Vector = getQ1VecReCtrFull();
  double Psi1       = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
  double Psi1Wgt    = transPsi1(Psi1);

  return Psi1Wgt;
}
//---------------------------------------------------------------------------------
// raw EP
void StEpdEpManager::initEpdSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1RawEastCent%d",iCent); // 2nd EP
    h_mEpdEp1RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1RawWestCent%d",iCent);
    h_mEpdEp1RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1RawFullCent%d",iCent);
    h_mEpdEp1RawFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1RawCorrCent%d",iCent);
    h_mEpdEp1RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpRaw(double Psi1East, double Psi1West, double Psi1Full)
{
  h_mEpdEp1RawEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1RawWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1RawFull[mCent9]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1RawCorr[mCent9]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1RawEast[iCent]->Write();
    h_mEpdEp1RawWest[iCent]->Write();
    h_mEpdEp1RawFull[iCent]->Write();
    h_mEpdEp1RawCorr[iCent]->Write();
  }
}

// phi weighted EP
void StEpdEpManager::initEpdSubEpWgt()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1WgtEastCent%d",iCent); // 1st EP
    h_mEpdEp1WgtEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1WgtWestCent%d",iCent);
    h_mEpdEp1WgtWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1WgtFullCent%d",iCent);
    h_mEpdEp1WgtFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1WgtCorrCent%d",iCent);
    h_mEpdEp1WgtCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpWgt(double Psi1East, double Psi1West, double Psi1Full)
{
  h_mEpdEp1WgtEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1WgtWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1WgtFull[mCent9]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1WgtCorr[mCent9]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpWgt()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1WgtEast[iCent]->Write();
    h_mEpdEp1WgtWest[iCent]->Write();
    h_mEpdEp1WgtFull[iCent]->Write();
    h_mEpdEp1WgtCorr[iCent]->Write();
  }
}

// recenter EP
void StEpdEpManager::initEpdSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1ReCtrEastCent%d",iCent); // 2nd EP
    h_mEpdEp1ReCtrEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1ReCtrWestCent%d",iCent);
    h_mEpdEp1ReCtrWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1ReCtrFullCent%d",iCent);
    h_mEpdEp1ReCtrFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1ReCtrCorrCent%d",iCent);
    h_mEpdEp1ReCtrCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpReCtr(double Psi1East, double Psi1West, double Psi1Full)
{
  h_mEpdEp1ReCtrEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1ReCtrWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1ReCtrFull[mCent9]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1ReCtrCorr[mCent9]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1ReCtrEast[iCent]->Write();
    h_mEpdEp1ReCtrWest[iCent]->Write();
    h_mEpdEp1ReCtrFull[iCent]->Write();
    h_mEpdEp1ReCtrCorr[iCent]->Write();
  }
}

// shift EP
void StEpdEpManager::initEpdSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1ShiftEastCent%d",iCent); // 2nd EP
    h_mEpdEp1ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1ShiftWestCent%d",iCent);
    h_mEpdEp1ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1ShiftFullCent%d",iCent);
    h_mEpdEp1ShiftFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1ShiftCorrCent%d",iCent);
    h_mEpdEp1ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpShift(double Psi1East, double Psi1West, double Psi1Full)
{
  h_mEpdEp1ShiftEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1ShiftWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1ShiftFull[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1ShiftCorr[mCent9]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1ShiftEast[iCent]->Write();
    h_mEpdEp1ShiftWest[iCent]->Write();
    h_mEpdEp1ShiftFull[iCent]->Write();
    h_mEpdEp1ShiftCorr[iCent]->Write();
  }
}

// shift Full EP
void StEpdEpManager::initEpdFullEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1ShiftFullCorrCent%d",iCent); // 1st EP
    h_mEpdEp1ShiftFullCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdFullEpShift(double Psi1FullCorr)
{
  h_mEpdEp1ShiftFullCorr[mCent9]->Fill(mRunIndex,Psi1FullCorr);
}

void StEpdEpManager::writeEpdFullEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1ShiftFullCorr[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
