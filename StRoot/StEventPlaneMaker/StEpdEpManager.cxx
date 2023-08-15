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
  mCent9    = -1;
  mRunIndex = -1;
  mVzBin    = -1;

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

  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
  {
    mQ1WgtGrpRawEast[iGrp] = 0.0;
    mQ1WgtGrpRawWest[iGrp] = 0.0;
    v_mQ1GrpRawEast[iGrp].Set(0.0,0.0);
    v_mQ1GrpRawWest[iGrp].Set(0.0,0.0);

    mQ1WgtGrpWgtEast[iGrp] = 0.0;
    mQ1WgtGrpWgtWest[iGrp] = 0.0;
    v_mQ1GrpWgtEast[iGrp].Set(0.0,0.0);
    v_mQ1GrpWgtWest[iGrp].Set(0.0,0.0);

    mQ1WgtGrpReCtrEast[iGrp] = 0.0;
    mQ1WgtGrpReCtrWest[iGrp] = 0.0;
    v_mQ1GrpReCtrEast[iGrp].Set(0.0,0.0);
    v_mQ1GrpReCtrWest[iGrp].Set(0.0,0.0);
  }
}

void StEpdEpManager::initEpdEpManager(int cent9, int runIndex, int vzBin)
{
  mCent9    = cent9;
  mRunIndex = runIndex;
  mVzBin    = vzBin;
}
//---------------------------------------------------------------------------------
// Utilities
TVector3 StEpdEpManager::getEpdCtrVec(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  TVector3 EpdPoint  = mEpdGeom->TileCenter(picoEpdHit->id()); // get a tile center position
  TVector3 EpdVector = EpdPoint - primVtx;

  return EpdVector;
}

TVector3 StEpdEpManager::getEpdRanVec(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  TVector3 EpdPoint  = mEpdGeom->RandomPointOnTile(picoEpdHit->id()); // get a random position within the tile
  TVector3 EpdVector = EpdPoint - primVtx;

  return EpdVector;
}

TVector2 StEpdEpManager::calq1Vector(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  TVector2 q1Vector(0.0,0.0);
  // TVector3 EpdVector = getEpdCtrVec(picoEpdHit, primVtx);
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

int StEpdEpManager::getEpdEpGrp(StPicoEpdHit *picoEpdHit)
{
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  int grpId = (ringId < 8) ? 0 : 1; // Group 0: 0-7 rings | Group 1: 8-15 rings
  if(ringId < 0 || ringId > 15) grpId = -1;

  return grpId;
}

double StEpdEpManager::transPsi1(double Psi1)
{
  double Psi1Corr = Psi1;
  if(Psi1 >  TMath::Pi()) Psi1Corr = Psi1 - TMath::TwoPi();
  if(Psi1 < -TMath::Pi()) Psi1Corr = Psi1 + TMath::TwoPi();

  return Psi1Corr;
}

bool StEpdEpManager::isPsi1InRange(double Psi1)
{
  if(Psi1 < -TMath::Pi() || Psi1 > TMath::Pi()) 
  {
    return false;
  }

  return true;
}
//---------------------------------------------------------------------------------
// Calculate Q1Vector
void StEpdEpManager::addHitSideRawEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId >= 0 && ringId < mNumRingsUsed)
  {
    v_mQ1SideRawEast  += tileWgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtSideRawEast += tileWgt;
  }
}

void StEpdEpManager::addHitSideRawWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId >= 0 && ringId < mNumRingsUsed)
  {
    v_mQ1SideRawWest  += tileWgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtSideRawWest += tileWgt;
  }
}

void StEpdEpManager::addHitSideWgtEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId >= 0 && ringId < mNumRingsUsed)
  {
    v_mQ1SideWgtEast  += wgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtSideWgtEast += wgt;
  }
}

void StEpdEpManager::addHitSideWgtWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId >= 0 && ringId < mNumRingsUsed)
  {
    v_mQ1SideWgtWest  += wgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtSideWgtWest += wgt;
  }
}

void StEpdEpManager::addHitSideReCtrEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId >= 0 && ringId < mNumRingsUsed)
  {
    v_mQ1SideReCtrEast  += wgt*(calq1Vector(picoEpdHit,primVtx) - getq1VecSideCtrEast());
    mQ1WgtSideReCtrEast += wgt;
  }
}

void StEpdEpManager::addHitSideReCtrWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId >= 0 && ringId < mNumRingsUsed)
  {
    v_mQ1SideReCtrWest  += wgt*(calq1Vector(picoEpdHit,primVtx) - getq1VecSideCtrWest());
    mQ1WgtSideReCtrWest += wgt;
  }
}

void StEpdEpManager::addHitGrpRawEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const int grpId = getEpdEpGrp(picoEpdHit);
  if(grpId >= 0)
  {
    v_mQ1GrpRawEast[grpId]  += tileWgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtGrpRawEast[grpId] += tileWgt;
  }
}

void StEpdEpManager::addHitGrpRawWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const int grpId = getEpdEpGrp(picoEpdHit);
  if(grpId >= 0)
  {
    v_mQ1GrpRawWest[grpId]  += tileWgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtGrpRawWest[grpId] += tileWgt;
  }
}

void StEpdEpManager::addHitGrpWgtEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int grpId = getEpdEpGrp(picoEpdHit);
  if(grpId >= 0)
  {
    v_mQ1GrpWgtEast[grpId]  += wgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtGrpWgtEast[grpId] += wgt;
  }
}

void StEpdEpManager::addHitGrpWgtWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int grpId = getEpdEpGrp(picoEpdHit);
  if(grpId >= 0)
  {
    v_mQ1GrpWgtWest[grpId]  += wgt*calq1Vector(picoEpdHit,primVtx);
    mQ1WgtGrpWgtWest[grpId] += wgt;
  }
}

void StEpdEpManager::addHitGrpReCtrEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int grpId = getEpdEpGrp(picoEpdHit);
  if(grpId >= 0)
  {
    v_mQ1GrpReCtrEast[grpId]  += wgt*(calq1Vector(picoEpdHit,primVtx) - getq1VecGrpCtrEast(grpId));
    mQ1WgtGrpReCtrEast[grpId] += wgt;
  }
}

void StEpdEpManager::addHitGrpReCtrWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;
  const int grpId = getEpdEpGrp(picoEpdHit);
  if(grpId >= 0)
  {
    v_mQ1GrpReCtrWest[grpId]  += wgt*(calq1Vector(picoEpdHit,primVtx) - getq1VecGrpCtrWest(grpId));
    mQ1WgtGrpReCtrWest[grpId] += wgt;
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
// ReCtrPar Correction
void StEpdEpManager::initEpdSideReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mEpdQ1SideReCtrXEastVz%d",iVz); // 1st EP
    p_mEpdQ1SideReCtrXEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mEpdQ1SideReCtrYEastVz%d",iVz);
    p_mEpdQ1SideReCtrYEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mEpdQ1SideReCtrXWestVz%d",iVz);
    p_mEpdQ1SideReCtrXWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mEpdQ1SideReCtrYWestVz%d",iVz);
    p_mEpdQ1SideReCtrYWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  }
}

void StEpdEpManager::fillEpdSideReCtrEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;

  const TVector2 q1Vector = calq1Vector(picoEpdHit,primVtx);
  const double q1x = q1Vector.X();
  const double q1y = q1Vector.Y();

  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId >= 0 && ringId < mNumRingsUsed)
  {
    p_mEpdQ1SideReCtrXEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1x,wgt);
    p_mEpdQ1SideReCtrYEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1y,wgt);
  }
}

void StEpdEpManager::fillEpdSideReCtrWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;

  const TVector2 q1Vector = calq1Vector(picoEpdHit,primVtx);
  const double q1x = q1Vector.X();
  const double q1y = q1Vector.Y();

  const int ringId = picoEpdHit->row() - 1; // convert ring Id to 0-15
  if(ringId >= 0 && ringId < mNumRingsUsed)
  {
    p_mEpdQ1SideReCtrXWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1x,wgt);
    p_mEpdQ1SideReCtrYWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,q1y,wgt);
  }
}

void StEpdEpManager::writeEpdSideReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mEpdQ1SideReCtrXEast[iVz]->Write();
    p_mEpdQ1SideReCtrYEast[iVz]->Write();
    p_mEpdQ1SideReCtrXWest[iVz]->Write();
    p_mEpdQ1SideReCtrYWest[iVz]->Write();
  }
}

void StEpdEpManager::readEpdSideReCtr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_EpdReCenterPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mReCtrPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mEpdQ1SideReCtrXEastVz%d",iVz); // 1st EP
    p_mEpdQ1SideReCtrXEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mEpdQ1SideReCtrYEastVz%d",iVz);
    p_mEpdQ1SideReCtrYEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mEpdQ1SideReCtrXWestVz%d",iVz);
    p_mEpdQ1SideReCtrXWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mEpdQ1SideReCtrYWestVz%d",iVz);
    p_mEpdQ1SideReCtrYWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
  }
}

TVector2 StEpdEpManager::getq1VecSideCtrEast()
{
  const int binX   = p_mEpdQ1SideReCtrXEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1x = p_mEpdQ1SideReCtrXEast[mVzBin]->GetBinContent(binX);

  const int binY   = p_mEpdQ1SideReCtrYEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1y = p_mEpdQ1SideReCtrYEast[mVzBin]->GetBinContent(binY);

  TVector2 q1Vector(0.0,0.0);
  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

TVector2 StEpdEpManager::getq1VecSideCtrWest()
{
  const int binX   = p_mEpdQ1SideReCtrXWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1x = p_mEpdQ1SideReCtrXWest[mVzBin]->GetBinContent(binX);

  const int binY   = p_mEpdQ1SideReCtrYWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1y = p_mEpdQ1SideReCtrYWest[mVzBin]->GetBinContent(binY);

  TVector2 q1Vector(0.0,0.0);
  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

void StEpdEpManager::initEpdGrpReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      std::string proName = Form("p_mEpdQ1Grp%dReCtrXEastVz%d",iGrp,iVz); // 1st EP Trk Ave ReCtr
      p_mEpdQ1GrpReCtrXEast[iVz][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1Grp%dReCtrYEastVz%d",iGrp,iVz);
      p_mEpdQ1GrpReCtrYEast[iVz][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mEpdQ1Grp%dReCtrXWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpReCtrXWest[iVz][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1Grp%dReCtrYWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpReCtrYWest[iVz][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mEpdQ1GrpAve%dReCtrXEastVz%d",iGrp,iVz); // 1st EP Evt Ave ReCtr
      p_mEpdQ1GrpAveReCtrXEast[iVz][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1GrpAve%dReCtrYEastVz%d",iGrp,iVz);
      p_mEpdQ1GrpAveReCtrYEast[iVz][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mEpdQ1GrpAve%dReCtrXWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpAveReCtrXWest[iVz][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1GrpAve%dReCtrYWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpAveReCtrYWest[iVz][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    }
  }
}

void StEpdEpManager::fillEpdGrpReCtrEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;

  const TVector2 q1Vector = calq1Vector(picoEpdHit,primVtx);
  const double q1x = q1Vector.X();
  const double q1y = q1Vector.Y();

  const int grpId = getEpdEpGrp(picoEpdHit);
  if(grpId >= 0)
  {
    p_mEpdQ1GrpReCtrXEast[mVzBin][grpId]->Fill((double)mRunIndex,(double)mCent9,q1x,wgt);
    p_mEpdQ1GrpReCtrYEast[mVzBin][grpId]->Fill((double)mRunIndex,(double)mCent9,q1y,wgt);
  }
}

void StEpdEpManager::fillEpdGrpReCtrWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx)
{
  const double tileWgt = getTileWeight(picoEpdHit);
  const double phiWgt  = getPhiWeight(picoEpdHit); // phiWgt sets to 1.0 if no correction file
  const double etaWgt  = getEtaWeight(picoEpdHit); // etaWgt sets to 1.0 if no correction file
  const double wgt     = tileWgt * phiWgt * etaWgt;

  const TVector2 q1Vector = calq1Vector(picoEpdHit,primVtx);
  const double q1x = q1Vector.X();
  const double q1y = q1Vector.Y();

  const int grpId = getEpdEpGrp(picoEpdHit);
  if(grpId >= 0)
  {
    p_mEpdQ1GrpReCtrXWest[mVzBin][grpId]->Fill((double)mRunIndex,(double)mCent9,q1x,wgt);
    p_mEpdQ1GrpReCtrYWest[mVzBin][grpId]->Fill((double)mRunIndex,(double)mCent9,q1y,wgt);
  }
}

void StEpdEpManager::fillEpdGrpAveReCtrEast(TVector2 Q1VecGrp, int grpId)
{
  const double Q1x = Q1VecGrp.X();
  const double Q1y = Q1VecGrp.Y();
  if(grpId >= 0)
  {
    p_mEpdQ1GrpAveReCtrXEast[mVzBin][grpId]->Fill((double)mRunIndex,(double)mCent9,Q1x);
    p_mEpdQ1GrpAveReCtrYEast[mVzBin][grpId]->Fill((double)mRunIndex,(double)mCent9,Q1y);
  }
}

void StEpdEpManager::fillEpdGrpAveReCtrWest(TVector2 Q1VecGrp, int grpId)
{
  const double Q1x = Q1VecGrp.X();
  const double Q1y = Q1VecGrp.Y();
  if(grpId >= 0)
  {
    p_mEpdQ1GrpAveReCtrXWest[mVzBin][grpId]->Fill((double)mRunIndex,(double)mCent9,Q1x);
    p_mEpdQ1GrpAveReCtrYWest[mVzBin][grpId]->Fill((double)mRunIndex,(double)mCent9,Q1y);
  }
}

void StEpdEpManager::writeEpdGrpReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      p_mEpdQ1GrpReCtrXEast[iVz][iGrp]->Write();
      p_mEpdQ1GrpReCtrYEast[iVz][iGrp]->Write();
      p_mEpdQ1GrpReCtrXWest[iVz][iGrp]->Write();
      p_mEpdQ1GrpReCtrYWest[iVz][iGrp]->Write();

      p_mEpdQ1GrpAveReCtrXEast[iVz][iGrp]->Write();
      p_mEpdQ1GrpAveReCtrYEast[iVz][iGrp]->Write();
      p_mEpdQ1GrpAveReCtrXWest[iVz][iGrp]->Write();
      p_mEpdQ1GrpAveReCtrYWest[iVz][iGrp]->Write();
    }
  }
}

void StEpdEpManager::readEpdGrpReCtr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_EpdReCenterPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mReCtrPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      std::string proName = Form("p_mEpdQ1Grp%dReCtrXEastVz%d",iGrp,iVz); // 1st EP
      p_mEpdQ1GrpReCtrXEast[iVz][iGrp] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1Grp%dReCtrYEastVz%d",iGrp,iVz);
      p_mEpdQ1GrpReCtrYEast[iVz][iGrp] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

      proName = Form("p_mEpdQ1Grp%dReCtrXWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpReCtrXWest[iVz][iGrp] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1Grp%dReCtrYWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpReCtrYWest[iVz][iGrp] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

      proName = Form("p_mEpdQ1GrpAve%dReCtrXEastVz%d",iGrp,iVz); // 1st EP
      p_mEpdQ1GrpAveReCtrXEast[iVz][iGrp] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1GrpAve%dReCtrYEastVz%d",iGrp,iVz);
      p_mEpdQ1GrpAveReCtrYEast[iVz][iGrp] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

      proName = Form("p_mEpdQ1GrpAve%dReCtrXWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpAveReCtrXWest[iVz][iGrp] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1GrpAve%dReCtrYWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpAveReCtrYWest[iVz][iGrp] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    }
  }
}

TVector2 StEpdEpManager::getq1VecGrpCtrEast(int grpId)
{
  const int binX   = p_mEpdQ1GrpReCtrXEast[mVzBin][grpId]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1x = p_mEpdQ1GrpReCtrXEast[mVzBin][grpId]->GetBinContent(binX);

  const int binY   = p_mEpdQ1GrpReCtrYEast[mVzBin][grpId]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1y = p_mEpdQ1GrpReCtrYEast[mVzBin][grpId]->GetBinContent(binY);

  TVector2 q1Vector(0.0,0.0);
  q1Vector.Set(q1x,q1y);

  return q1Vector;
}

TVector2 StEpdEpManager::getq1VecGrpCtrWest(int grpId)
{
  const int binX   = p_mEpdQ1GrpReCtrXWest[mVzBin][grpId]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1x = p_mEpdQ1GrpReCtrXWest[mVzBin][grpId]->GetBinContent(binX);

  const int binY   = p_mEpdQ1GrpReCtrYWest[mVzBin][grpId]->FindBin((double)mRunIndex,(double)mCent9);
  const double q1y = p_mEpdQ1GrpReCtrYWest[mVzBin][grpId]->GetBinContent(binY);

  TVector2 q1Vector(0.0,0.0);
  q1Vector.Set(q1x,q1y);

  return q1Vector;
}
//---------------------------------------------------------------------------------
// Shift Correction
void StEpdEpManager::initEpdSideShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1SideShiftCos%dEastVz%d",iShift,iVz); // 1st EP
      p_mEpdQ1SideShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1SideShiftSin%dEastVz%d",iShift,iVz);
      p_mEpdQ1SideShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mEpdQ1SideShiftCos%dWestVz%d",iShift,iVz);
      p_mEpdQ1SideShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1SideShiftSin%dWestVz%d",iShift,iVz);
      p_mEpdQ1SideShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    }
  }
}

void StEpdEpManager::fillEpdSideShiftEast()
{
  TVector2 Q1VecSide = getQ1VecSideReCtrEast();
  if(Q1VecSide.Mod() > 0.0)
  {
    const double Psi1 = TMath::ATan2(Q1VecSide.Y(),Q1VecSide.X()); // -pi to pi
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
      const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
      p_mEpdQ1SideShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
      p_mEpdQ1SideShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
    }
  }
}

void StEpdEpManager::fillEpdSideShiftWest()
{
  TVector2 Q1Vector = getQ1VecSideReCtrWest();
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X());
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
      const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
      p_mEpdQ1SideShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
      p_mEpdQ1SideShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
    }
  }
}

void StEpdEpManager::writeEpdSideShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mEpdQ1SideShiftCosEast[iVz][iShift]->Write();
      p_mEpdQ1SideShiftSinEast[iVz][iShift]->Write();
      p_mEpdQ1SideShiftCosWest[iVz][iShift]->Write();
      p_mEpdQ1SideShiftSinWest[iVz][iShift]->Write();
    }
  }
}

void StEpdEpManager::readEpdSideShift()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftPar/file_EpdShiftPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1SideShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mEpdQ1SideShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1SideShiftSin%dEastVz%d",iShift,iVz);
      p_mEpdQ1SideShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mEpdQ1SideShiftCos%dWestVz%d",iShift,iVz);
      p_mEpdQ1SideShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1SideShiftSin%dWestVz%d",iShift,iVz);
      p_mEpdQ1SideShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

double StEpdEpManager::getPsi1SideShiftEast()
{
  double Psi1Shift = -999.9;
  TVector2 Q1Vector = getQ1VecSideReCtrEast();
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mEpdQ1SideShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mEpdQ1SideShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mEpdQ1SideShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mEpdQ1SideShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1ReCtr));
    }

    double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StEpdEpManager::getPsi1SideShiftWest()
{
  double Psi1Shift = -999.9;
  TVector2 Q1Vector = getQ1VecSideReCtrWest();
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mEpdQ1SideShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mEpdQ1SideShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mEpdQ1SideShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mEpdQ1SideShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1ReCtr));
    }

    double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StEpdEpManager::getPsi1SideShiftFull()
{
  double Psi1Shift = -999.9;
  TVector2 Q1Vector = getQ1VecSideShiftFull();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1ShiftRaw = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

TVector2 StEpdEpManager::getQ1VecSideShiftEast()
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1East = getPsi1SideShiftEast();
  if(isPsi1InRange(Psi1East))
  {
    const double Q1x = TMath::Cos(1.0*Psi1East);
    const double Q1y = TMath::Sin(1.0*Psi1East);
    Q1Vector.Set(Q1x,Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecSideShiftWest()
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1West = getPsi1SideShiftWest();
  if(isPsi1InRange(Psi1West))
  {
    const double Q1x = TMath::Cos(1.0*Psi1West);
    const double Q1y = TMath::Sin(1.0*Psi1West);
    Q1Vector.Set(Q1x,Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecSideShiftFull()
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecSideShiftEast();
  TVector2 Q1VecWest = getQ1VecSideShiftWest();
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}

void StEpdEpManager::initEpdGrpShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	std::string proName = Form("p_mEpdQ1Grp%dShiftCos%dEastVz%d",iGrp,iShift,iVz); // 1st EP
	p_mEpdQ1GrpShiftCosEast[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
	proName = Form("p_mEpdQ1Grp%dShiftSin%dEastVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftSinEast[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

	proName = Form("p_mEpdQ1Grp%dShiftCos%dWestVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftCosWest[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
	proName = Form("p_mEpdQ1Grp%dShiftSin%dWestVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftSinWest[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

	proName = Form("p_mEpdQ1GrpAve%dShiftCos%dEastVz%d",iGrp,iShift,iVz); // 1st EP
	p_mEpdQ1GrpAveShiftCosEast[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
	proName = Form("p_mEpdQ1GrpAve%dShiftSin%dEastVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpAveShiftSinEast[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

	proName = Form("p_mEpdQ1GrpAve%dShiftCos%dWestVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpAveShiftCosWest[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
	proName = Form("p_mEpdQ1GrpAve%dShiftSin%dWestVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpAveShiftSinWest[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      }
    }
  }
}

void StEpdEpManager::fillEpdGrpShiftEast(int grpId)
{
  TVector2 Q1VecGrp = getQ1VecGrpReCtrEast(grpId);
  const double Psi1 = TMath::ATan2(Q1VecGrp.Y(),Q1VecGrp.X()); // -pi to pi
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
    const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
    p_mEpdQ1GrpShiftCosEast[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
    p_mEpdQ1GrpShiftSinEast[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
  }
}

void StEpdEpManager::fillEpdGrpShiftWest(int grpId)
{
  TVector2 Q1VecGrp = getQ1VecGrpReCtrWest(grpId);
  const double Psi1 = TMath::ATan2(Q1VecGrp.Y(),Q1VecGrp.X()); // -pi to pi
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
    const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
    p_mEpdQ1GrpShiftCosWest[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
    p_mEpdQ1GrpShiftSinWest[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
  }
}

void StEpdEpManager::fillEpdGrpAveShiftEast(int grpId)
{
  TVector2 Q1VecGrp = getQ1VecGrpAveReCtrEast(grpId);
  const double Psi1 = TMath::ATan2(Q1VecGrp.Y(),Q1VecGrp.X()); // -pi to pi
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
    const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
    p_mEpdQ1GrpAveShiftCosEast[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
    p_mEpdQ1GrpAveShiftSinEast[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
  }
}

void StEpdEpManager::fillEpdGrpAveShiftWest(int grpId)
{
  TVector2 Q1VecGrp = getQ1VecGrpAveReCtrWest(grpId);
  const double Psi1 = TMath::ATan2(Q1VecGrp.Y(),Q1VecGrp.X()); // -pi to pi
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
    const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
    p_mEpdQ1GrpAveShiftCosWest[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
    p_mEpdQ1GrpAveShiftSinWest[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
  }
}

void StEpdEpManager::writeEpdGrpShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	p_mEpdQ1GrpShiftCosEast[iVz][iShift][iGrp]->Write();
	p_mEpdQ1GrpShiftSinEast[iVz][iShift][iGrp]->Write();
	p_mEpdQ1GrpShiftCosWest[iVz][iShift][iGrp]->Write();
	p_mEpdQ1GrpShiftSinWest[iVz][iShift][iGrp]->Write();

	p_mEpdQ1GrpAveShiftCosEast[iVz][iShift][iGrp]->Write();
	p_mEpdQ1GrpAveShiftSinEast[iVz][iShift][iGrp]->Write();
	p_mEpdQ1GrpAveShiftCosWest[iVz][iShift][iGrp]->Write();
	p_mEpdQ1GrpAveShiftSinWest[iVz][iShift][iGrp]->Write();
      }
    }
  }
}

void StEpdEpManager::readEpdGrpShift()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftPar/file_EpdShiftPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	std::string proName = Form("p_mEpdQ1Grp%dShiftCos%dEastVz%d",iGrp,iShift,iVz); // 1st EP
	p_mEpdQ1GrpShiftCosEast[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
	proName = Form("p_mEpdQ1Grp%dShiftSin%dEastVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftSinEast[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

	proName = Form("p_mEpdQ1Grp%dShiftCos%dWestVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftCosWest[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
	proName = Form("p_mEpdQ1Grp%dShiftSin%dWestVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftSinWest[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

	proName = Form("p_mEpdQ1GrpAve%dShiftCos%dEastVz%d",iGrp,iShift,iVz); // 1st EP
	p_mEpdQ1GrpAveShiftCosEast[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
	proName = Form("p_mEpdQ1GrpAve%dShiftSin%dEastVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpAveShiftSinEast[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

	proName = Form("p_mEpdQ1GrpAve%dShiftCos%dWestVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpAveShiftCosWest[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
	proName = Form("p_mEpdQ1GrpAve%dShiftSin%dWestVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpAveShiftSinWest[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      }
    }
  }
}

double StEpdEpManager::getPsi1GrpShiftEast(int grpId)
{
  double Psi1Shift = -999.9;
  TVector2 Q1Vector = getQ1VecGrpReCtrEast(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mEpdQ1GrpShiftCosEast[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mEpdQ1GrpShiftCosEast[mVzBin][iShift][grpId]->GetBinContent(binCos);

      const int binSin     = p_mEpdQ1GrpShiftSinEast[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mEpdQ1GrpShiftSinEast[mVzBin][iShift][grpId]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1ReCtr));
    }

    double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StEpdEpManager::getPsi1GrpShiftWest(int grpId)
{
  double Psi1Shift = -999.9;
  TVector2 Q1Vector = getQ1VecGrpReCtrWest(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mEpdQ1GrpShiftCosWest[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mEpdQ1GrpShiftCosWest[mVzBin][iShift][grpId]->GetBinContent(binCos);

      const int binSin     = p_mEpdQ1GrpShiftSinWest[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mEpdQ1GrpShiftSinWest[mVzBin][iShift][grpId]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1ReCtr));
    }

    double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StEpdEpManager::getPsi1GrpShiftFull(int grpId)
{
  double Psi1Shift = -999.9;
  TVector2 Q1Vector = getQ1VecGrpShiftFull(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1ShiftRaw = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StEpdEpManager::getPsi1GrpAveShiftEast(int grpId)
{
  double Psi1Shift = -999.9;
  TVector2 Q1Vector = getQ1VecGrpAveReCtrEast(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mEpdQ1GrpAveShiftCosEast[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mEpdQ1GrpAveShiftCosEast[mVzBin][iShift][grpId]->GetBinContent(binCos);

      const int binSin     = p_mEpdQ1GrpAveShiftSinEast[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mEpdQ1GrpAveShiftSinEast[mVzBin][iShift][grpId]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1ReCtr));
    }

    double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StEpdEpManager::getPsi1GrpAveShiftWest(int grpId)
{
  double Psi1Shift = -999.9;
  TVector2 Q1Vector = getQ1VecGrpAveReCtrWest(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1ReCtr = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mEpdQ1GrpAveShiftCosWest[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mEpdQ1GrpAveShiftCosWest[mVzBin][iShift][grpId]->GetBinContent(binCos);

      const int binSin     = p_mEpdQ1GrpAveShiftSinWest[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mEpdQ1GrpAveShiftSinWest[mVzBin][iShift][grpId]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1ReCtr)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1ReCtr));
    }

    double Psi1ShiftRaw = Psi1ReCtr + deltaPsi1;
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

double StEpdEpManager::getPsi1GrpAveShiftFull(int grpId)
{
  double Psi1Shift = -999.9;
  TVector2 Q1Vector = getQ1VecGrpAveShiftFull(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1ShiftRaw = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Shift = transPsi1(Psi1ShiftRaw);
  }

  return Psi1Shift;
}

TVector2 StEpdEpManager::getQ1VecGrpShiftEast(int grpId)
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1East = getPsi1GrpShiftEast(grpId);
  if(isPsi1InRange(Psi1East))
  {
    const double Q1x = TMath::Cos(1.0*Psi1East);
    const double Q1y = TMath::Sin(1.0*Psi1East);
    Q1Vector.Set(Q1x,Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpShiftWest(int grpId)
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1West = getPsi1GrpShiftWest(grpId);
  if(isPsi1InRange(Psi1West))
  {
    const double Q1x = TMath::Cos(1.0*Psi1West);
    const double Q1y = TMath::Sin(1.0*Psi1West);
    Q1Vector.Set(Q1x,Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpShiftFull(int grpId)
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecGrpShiftEast(grpId);
  TVector2 Q1VecWest = getQ1VecGrpShiftWest(grpId);
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}

TVector2 StEpdEpManager::getQ1VecGrpAveShiftEast(int grpId)
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1East = getPsi1GrpAveShiftEast(grpId);
  if(isPsi1InRange(Psi1East))
  {
    const double Q1x = TMath::Cos(1.0*Psi1East);
    const double Q1y = TMath::Sin(1.0*Psi1East);
    Q1Vector.Set(Q1x,Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpAveShiftWest(int grpId)
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1West = getPsi1GrpAveShiftWest(grpId);
  if(isPsi1InRange(Psi1West))
  {
    const double Q1x = TMath::Cos(1.0*Psi1West);
    const double Q1y = TMath::Sin(1.0*Psi1West);
    Q1Vector.Set(Q1x,Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpAveShiftFull(int grpId)
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecGrpAveShiftEast(grpId);
  TVector2 Q1VecWest = getQ1VecGrpAveShiftWest(grpId);
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}
//---------------------------------------------------------------------------------
// Shift Correction Full EP
void StEpdEpManager::initEpdSideShiftFull()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1SideShiftCos%dFullVz%d",iShift,iVz);
      p_mEpdQ1SideShiftCosFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mEpdQ1SideShiftSin%dFullVz%d",iShift,iVz);
      p_mEpdQ1SideShiftSinFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    }
  }
}

void StEpdEpManager::fillEpdSideShiftFull()
{
  TVector2 Q1VecSide = getQ1VecSideShiftFull();
  if(Q1VecSide.Mod() > 0.0)
  {
    const double Psi1 = TMath::ATan2(Q1VecSide.Y(),Q1VecSide.X());
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
      const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
      p_mEpdQ1SideShiftCosFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
      p_mEpdQ1SideShiftSinFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
    }
  }
}

void StEpdEpManager::writeEpdSideShiftFull()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mEpdQ1SideShiftCosFull[iVz][iShift]->Write();
      p_mEpdQ1SideShiftSinFull[iVz][iShift]->Write();
    }
  }
}

void StEpdEpManager::readEpdSideShiftFull()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftPar/file_EpdShiftParFull_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < 20; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1SideShiftCos%dFullVz%d",iShift,iVz);
      p_mEpdQ1SideShiftCosFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mEpdQ1SideShiftSin%dFullVz%d",iShift,iVz);
      p_mEpdQ1SideShiftSinFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

double StEpdEpManager::getPsi1SideShiftFullCorr()
{
  double Psi1ShiftCorr = -999.9;
  TVector2 Q1Vector = getQ1VecSideShiftFull();
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1Shift = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mEpdQ1SideShiftCosFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mEpdQ1SideShiftCosFull[mVzBin][iShift]->GetBinContent(binCos);

      const int binSin     = p_mEpdQ1SideShiftSinFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mEpdQ1SideShiftSinFull[mVzBin][iShift]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1Shift)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1Shift));
    }

    double Psi1ShiftCorrRaw = Psi1Shift + deltaPsi1;
    Psi1ShiftCorr = transPsi1(Psi1ShiftCorrRaw);
  }

  return Psi1ShiftCorr;
}

TVector2 StEpdEpManager::getQ1VecSideShiftFullCorr()
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1FullCorr = getPsi1SideShiftFullCorr();
  if(isPsi1InRange(Psi1FullCorr))
  {
    const double Q1x = TMath::Cos(1.0*Psi1FullCorr);
    const double Q1y = TMath::Sin(1.0*Psi1FullCorr);
    Q1Vector.Set(Q1x,Q1y);
  }

  return Q1Vector;
}

void StEpdEpManager::initEpdGrpShiftFull()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	std::string proName = Form("p_mEpdQ1Grp%dShiftCos%dFullVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftCosFull[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
	proName = Form("p_mEpdQ1Grp%dShiftSin%dFullVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftSinFull[iVz][iShift][iGrp] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      }
    }
  }
}

void StEpdEpManager::fillEpdGrpShiftFull(int grpId)
{
  TVector2 Q1VecGrp = getQ1VecGrpShiftFull(grpId);
  if(Q1VecGrp.Mod() > 0.0)
  {
    const double Psi1 = TMath::ATan2(Q1VecGrp.Y(),Q1VecGrp.X());
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
    {
      const double Psi1Cos = TMath::Cos(((double)iShift+1.0)*Psi1);
      const double Psi1Sin = TMath::Sin(((double)iShift+1.0)*Psi1);
      p_mEpdQ1GrpShiftCosFull[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Cos);
      p_mEpdQ1GrpShiftSinFull[mVzBin][iShift][grpId]->Fill((double)mRunIndex,(double)mCent9,Psi1Sin);
    }
  }
}

void StEpdEpManager::writeEpdGrpShiftFull()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	p_mEpdQ1GrpShiftCosFull[iVz][iShift][iGrp]->Write();
	p_mEpdQ1GrpShiftSinFull[iVz][iShift][iGrp]->Write();
      }
    }
  }
}

void StEpdEpManager::readEpdGrpShiftFull()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftPar/file_EpdShiftParFull_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mShiftPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < 20; ++iShift) // Shift Order
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	std::string proName = Form("p_mEpdQ1Grp%dShiftCos%dFullVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftCosFull[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
	proName = Form("p_mEpdQ1Grp%dShiftSin%dFullVz%d",iGrp,iShift,iVz);
	p_mEpdQ1GrpShiftSinFull[iVz][iShift][iGrp] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      }
    }
  }
}

double StEpdEpManager::getPsi1GrpShiftFullCorr(int grpId)
{
  double Psi1ShiftCorr = -999.9;
  TVector2 Q1Vector = getQ1VecGrpShiftFull(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    const double Psi1Shift = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    double deltaPsi1 = 0.0;
    for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
    {
      const int binCos     = p_mEpdQ1GrpShiftCosFull[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanCos = p_mEpdQ1GrpShiftCosFull[mVzBin][iShift][grpId]->GetBinContent(binCos);

      const int binSin     = p_mEpdQ1GrpShiftSinFull[mVzBin][iShift][grpId]->FindBin((double)mRunIndex,(double)mCent9);
      const double meanSin = p_mEpdQ1GrpShiftSinFull[mVzBin][iShift][grpId]->GetBinContent(binSin);

      deltaPsi1 += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*Psi1Shift)+meanCos*TMath::Sin(((double)iShift+1.0)*Psi1Shift));
    }

    double Psi1ShiftCorrRaw = Psi1Shift + deltaPsi1;
    Psi1ShiftCorr = transPsi1(Psi1ShiftCorrRaw);
  }

  return Psi1ShiftCorr;
}

TVector2 StEpdEpManager::getQ1VecGrpShiftFullCorr(int grpId)
{
  TVector2 Q1Vector(0.0,0.0);
  const double Psi1FullCorr = getPsi1GrpShiftFullCorr(grpId);
  if(isPsi1InRange(Psi1FullCorr))
  {
    const double Q1x = TMath::Cos(1.0*Psi1FullCorr);
    const double Q1y = TMath::Sin(1.0*Psi1FullCorr);
    Q1Vector.Set(Q1x,Q1y);
  }

  return Q1Vector;
}
//---------------------------------------------------------------------------------
// Sub Event Plane Resolution
void StEpdEpManager::initEpdSideResolution()
{
  p_mEpdSubEp1SideRes = new TProfile("p_mEpdSubEp1SideRes","p_mEpdSubEp1SideRes",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
}

void StEpdEpManager::fillEpdSideResolution(double Psi1East, double Psi1West)
{
  double res1Sub = TMath::Cos(1.0*(Psi1West-Psi1East));
  p_mEpdSubEp1SideRes->Fill((double)mCent9,res1Sub);
}

void StEpdEpManager::writeEpdSideResolution()
{
  p_mEpdSubEp1SideRes->Write();
}

void StEpdEpManager::readEpdSideResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_EpdEpResolution_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mEpdSubEp1SideRes = (TProfile*)file_mResolution->Get("p_mEpdSubEp1SideRes");

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    mEpdSubEp1SideResVal[iCent]  = 0.0;
    mEpdSubEp1SideResErr[iCent]  = 0.0;
    mEpdFullEp1SideResVal[iCent] = 0.0;
    mEpdFullEp1SideResErr[iCent] = 0.0;
  }

  TF1 *f_EpdEpResFull = new TF1("f_EpdEpResFull",funcEpdEpResFull,0,10,0);
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    double valRes1Sub  = -999.9;
    double errRes1Sub  = 1.0;
    double valRes1Full = -999.9;
    double errRes1Full = 1.0;
    double valRes1Raw  = p_mEpdSubEp1SideRes->GetBinContent(iCent+1);
    double errRes1Raw  = p_mEpdSubEp1SideRes->GetBinError(iCent+1);
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

    mEpdSubEp1SideResVal[iCent]  = valRes1Sub;
    mEpdSubEp1SideResErr[iCent]  = errRes1Sub;
    mEpdFullEp1SideResVal[iCent] = valRes1Full;
    mEpdFullEp1SideResErr[iCent] = errRes1Full;
  }
  file_mResolution->Close();
}

double StEpdEpManager::getEpdSubEp1SideResVal(int cent9)
{
  return mEpdSubEp1SideResVal[cent9];
}

double StEpdEpManager::getEpdSubEp1SideResErr(int cent9)
{
  return mEpdSubEp1SideResErr[cent9];
}

double StEpdEpManager::getEpdFullEp1SideResVal(int cent9)
{
  return mEpdFullEp1SideResVal[cent9];
}

double StEpdEpManager::getEpdFullEp1SideResErr(int cent9)
{
  return mEpdFullEp1SideResErr[cent9];
}

void StEpdEpManager::initEpdGrpResolution()
{
  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
  {
    std::string proName = Form("p_mEpdSubEp1Grp%dRes",iGrp);
    p_mEpdSubEp1GrpRes[iGrp] = new TProfile(proName.c_str(),proName.c_str(),mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  }
}

void StEpdEpManager::fillEpdGrpResolution(double Psi1East, double Psi1West, int grpId)
{
  double res1Sub = TMath::Cos(1.0*(Psi1West-Psi1East));
  p_mEpdSubEp1GrpRes[grpId]->Fill((double)mCent9,res1Sub);
}

void StEpdEpManager::writeEpdGrpResolution()
{
  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
  {
    p_mEpdSubEp1GrpRes[iGrp]->Write();
  }
}

void StEpdEpManager::readEpdGrpResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_EpdEpResolution_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
  {
    std::string proName = Form("p_mEpdSubEp1Grp%dRes",iGrp);
    p_mEpdSubEp1GrpRes[iGrp] = (TProfile*)file_mResolution->Get(proName.c_str());
  }

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      mEpdSubEp1GrpResVal[iCent][iGrp]  = 0.0;
      mEpdSubEp1GrpResErr[iCent][iGrp]  = 0.0;
      mEpdFullEp1GrpResVal[iCent][iGrp] = 0.0;
      mEpdFullEp1GrpResErr[iCent][iGrp] = 0.0;
    }
  }

  TF1 *f_EpdEpResFull = new TF1("f_EpdEpResFull",funcEpdEpResFull,0,10,0);
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      double valRes1Sub  = -999.9;
      double errRes1Sub  = 1.0;
      double valRes1Full = -999.9;
      double errRes1Full = 1.0;
      double valRes1Raw  = p_mEpdSubEp1GrpRes[iGrp]->GetBinContent(iCent+1);
      double errRes1Raw  = p_mEpdSubEp1GrpRes[iGrp]->GetBinError(iCent+1);
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

      mEpdSubEp1GrpResVal[iCent][iGrp]  = valRes1Sub;
      mEpdSubEp1GrpResErr[iCent][iGrp]  = errRes1Sub;
      mEpdFullEp1GrpResVal[iCent][iGrp] = valRes1Full;
      mEpdFullEp1GrpResErr[iCent][iGrp] = errRes1Full;
    }
  }
  file_mResolution->Close();
}

double StEpdEpManager::getEpdSubEp1GrpResVal(int cent9, int grpId)
{
  return mEpdSubEp1GrpResVal[cent9][grpId];
}

double StEpdEpManager::getEpdSubEp1GrpResErr(int cent9, int grpId)
{
  return mEpdSubEp1GrpResErr[cent9][grpId];
}

double StEpdEpManager::getEpdFullEp1GrpResVal(int cent9, int grpId)
{
  return mEpdFullEp1GrpResVal[cent9][grpId];
}

double StEpdEpManager::getEpdFullEp1GrpResErr(int cent9, int grpId)
{
  return mEpdFullEp1GrpResErr[cent9][grpId];
}
//---------------------------------------------------------------------------------
// Charged Hadron Directed Flow
void StEpdEpManager::initEpdSubEpSideFlow()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string proName = Form("p_mEpdSubEpSideV1Cent%d",iCent);
    p_mEpdSubEpSideV1[iCent] = new TProfile(proName.c_str(),proName.c_str(),100,-10.0,10.0);
  }
}

void StEpdEpManager::fillEpdSubEpSideV1(double eta, double v1, double reweight)
{
  p_mEpdSubEpSideV1[mCent9]->Fill(eta, v1, reweight);
}

void StEpdEpManager::writeEpdSubEpSideFlow()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    p_mEpdSubEpSideV1[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
// Q1Vector
TVector2 StEpdEpManager::getQ1VecSideRawEast()
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtSideRawEast > 0.0)
  {
    const double Q1x = -1.0*v_mQ1SideRawEast.X()/mQ1WgtSideRawEast; // flip the sign for Q1VecEast
    const double Q1y = -1.0*v_mQ1SideRawEast.Y()/mQ1WgtSideRawEast;
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecSideRawWest()
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtSideRawWest > 0.0)
  {
    const double Q1x = v_mQ1SideRawWest.X()/mQ1WgtSideRawWest;
    const double Q1y = v_mQ1SideRawWest.Y()/mQ1WgtSideRawWest;
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecSideRawFull()
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecSideRawEast();
  TVector2 Q1VecWest = getQ1VecSideRawWest();
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}

TVector2 StEpdEpManager::getQ1VecSideWgtEast()
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtSideWgtEast > 0.0)
  {
    const double Q1x = -1.0*v_mQ1SideWgtEast.X()/mQ1WgtSideWgtEast; // flip the sign for Q1VecEast
    const double Q1y = -1.0*v_mQ1SideWgtEast.Y()/mQ1WgtSideWgtEast;
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecSideWgtWest()
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtSideWgtWest > 0.0)
  {
    const double Q1x = v_mQ1SideWgtWest.X()/mQ1WgtSideWgtWest;
    const double Q1y = v_mQ1SideWgtWest.Y()/mQ1WgtSideWgtWest;
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecSideWgtFull()
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecSideWgtEast();
  TVector2 Q1VecWest = getQ1VecSideWgtWest();
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}

TVector2 StEpdEpManager::getQ1VecSideReCtrEast()
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtSideReCtrEast > 0.0)
  {
    const double Q1x = -1.0*v_mQ1SideReCtrEast.X()/mQ1WgtSideReCtrEast; // flip the sign for Q1VecEast
    const double Q1y = -1.0*v_mQ1SideReCtrEast.Y()/mQ1WgtSideReCtrEast;
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecSideReCtrWest()
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtSideReCtrWest > 0.0)
  {
    const double Q1x = v_mQ1SideReCtrWest.X()/mQ1WgtSideReCtrWest;
    const double Q1y = v_mQ1SideReCtrWest.Y()/mQ1WgtSideReCtrWest;
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecSideReCtrFull()
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecSideReCtrEast();
  TVector2 Q1VecWest = getQ1VecSideReCtrWest();
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}

double StEpdEpManager::getPsi1SideRawEast()
{
  double Psi1Raw = -999.9;
  TVector2 Q1Vector = getQ1VecSideRawEast();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Raw = transPsi1(Psi1);
  }

  return Psi1Raw;
}

double StEpdEpManager::getPsi1SideRawWest()
{
  double Psi1Raw = -999.9;
  TVector2 Q1Vector = getQ1VecSideRawWest();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Raw = transPsi1(Psi1);
  }

  return Psi1Raw;
}

double StEpdEpManager::getPsi1SideRawFull()
{
  double Psi1Raw = -999.9;
  TVector2 Q1Vector = getQ1VecSideRawFull();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Raw = transPsi1(Psi1);
  }

  return Psi1Raw;
}

double StEpdEpManager::getPsi1SideWgtEast()
{
  double Psi1Wgt = -999.9;
  TVector2 Q1Vector = getQ1VecSideWgtEast();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Wgt = transPsi1(Psi1);
  }

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1SideWgtWest()
{
  double Psi1Wgt = -999.9;
  TVector2 Q1Vector = getQ1VecSideWgtWest();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Wgt = transPsi1(Psi1);
  }

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1SideWgtFull()
{
  double Psi1Wgt = -999.9;
  TVector2 Q1Vector = getQ1VecSideWgtFull();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Wgt = transPsi1(Psi1);
  }

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1SideReCtrEast()
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector = getQ1VecSideReCtrEast();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

double StEpdEpManager::getPsi1SideReCtrWest()
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector = getQ1VecSideReCtrWest();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

double StEpdEpManager::getPsi1SideReCtrFull()
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector = getQ1VecSideReCtrFull();
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

TVector2 StEpdEpManager::getQ1VecGrpRawEast(int grpId)
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtGrpRawEast[grpId] > 0.0)
  {
    const double Q1x = -1.0*v_mQ1GrpRawEast[grpId].X()/mQ1WgtGrpRawEast[grpId]; // flip the sign for Q1VecEast
    const double Q1y = -1.0*v_mQ1GrpRawEast[grpId].Y()/mQ1WgtGrpRawEast[grpId];
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpRawWest(int grpId)
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtGrpRawWest[grpId] > 0.0)
  {
    const double Q1x = v_mQ1GrpRawWest[grpId].X()/mQ1WgtGrpRawWest[grpId];
    const double Q1y = v_mQ1GrpRawWest[grpId].Y()/mQ1WgtGrpRawWest[grpId];
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpRawFull(int grpId)
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecGrpRawEast(grpId);
  TVector2 Q1VecWest = getQ1VecGrpRawWest(grpId);
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}

TVector2 StEpdEpManager::getQ1VecGrpWgtEast(int grpId)
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtGrpWgtEast[grpId] > 0.0)
  {
    const double Q1x = -1.0*v_mQ1GrpWgtEast[grpId].X()/mQ1WgtGrpWgtEast[grpId]; // flip the sign for Q1VecEast
    const double Q1y = -1.0*v_mQ1GrpWgtEast[grpId].Y()/mQ1WgtGrpWgtEast[grpId];
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpWgtWest(int grpId)
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtGrpWgtWest[grpId] > 0.0)
  {
    const double Q1x = v_mQ1GrpWgtWest[grpId].X()/mQ1WgtGrpWgtWest[grpId];
    const double Q1y = v_mQ1GrpWgtWest[grpId].Y()/mQ1WgtGrpWgtWest[grpId];
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpWgtFull(int grpId)
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecGrpWgtEast(grpId);
  TVector2 Q1VecWest = getQ1VecGrpWgtWest(grpId);
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}

TVector2 StEpdEpManager::getQ1VecGrpReCtrEast(int grpId)
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtGrpReCtrEast[grpId] > 0.0)
  {
    const double Q1x = -1.0*v_mQ1GrpReCtrEast[grpId].X()/mQ1WgtGrpReCtrEast[grpId]; // flip the sign for Q1VecEast
    const double Q1y = -1.0*v_mQ1GrpReCtrEast[grpId].Y()/mQ1WgtGrpReCtrEast[grpId];
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpReCtrWest(int grpId)
{
  TVector2 Q1Vector(0.0, 0.0);
  if(mQ1WgtGrpReCtrWest[grpId] > 0.0)
  {
    const double Q1x = v_mQ1GrpReCtrWest[grpId].X()/mQ1WgtGrpReCtrWest[grpId];
    const double Q1y = v_mQ1GrpReCtrWest[grpId].Y()/mQ1WgtGrpReCtrWest[grpId];
    Q1Vector.Set(Q1x, Q1y);
  }

  return Q1Vector;
}

TVector2 StEpdEpManager::getQ1VecGrpReCtrFull(int grpId)
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecGrpReCtrEast(grpId);
  TVector2 Q1VecWest = getQ1VecGrpReCtrWest(grpId);
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}

TVector2 StEpdEpManager::getQ1VecGrpAveReCtrEast(int grpId)
{
  TVector2 Q1VecReCtr(0.0,0.0);
  TVector2 Q1VecRaw = getQ1VecGrpRawEast(grpId);
  if(Q1VecRaw.Mod() > 0.0)
  {
    const int binX      = p_mEpdQ1GrpAveReCtrXEast[mVzBin][grpId]->FindBin((double)mRunIndex,(double)mCent9);
    const double Q1xCtr = p_mEpdQ1GrpAveReCtrXEast[mVzBin][grpId]->GetBinContent(binX);

    const int binY      = p_mEpdQ1GrpAveReCtrYEast[mVzBin][grpId]->FindBin((double)mRunIndex,(double)mCent9);
    const double Q1yCtr = p_mEpdQ1GrpAveReCtrYEast[mVzBin][grpId]->GetBinContent(binY);

    const double Q1x = Q1VecRaw.X() - Q1xCtr;
    const double Q1y = Q1VecRaw.Y() - Q1yCtr;

    Q1VecReCtr.Set(Q1x, Q1y);
  }

  return Q1VecReCtr;
}

TVector2 StEpdEpManager::getQ1VecGrpAveReCtrWest(int grpId)
{
  TVector2 Q1VecReCtr(0.0,0.0);
  TVector2 Q1VecRaw = getQ1VecGrpRawWest(grpId);
  if(Q1VecRaw.Mod() > 0.0)
  {
    const int binX      = p_mEpdQ1GrpAveReCtrXWest[mVzBin][grpId]->FindBin((double)mRunIndex,(double)mCent9);
    const double Q1xCtr = p_mEpdQ1GrpAveReCtrXWest[mVzBin][grpId]->GetBinContent(binX);

    const int binY      = p_mEpdQ1GrpAveReCtrYWest[mVzBin][grpId]->FindBin((double)mRunIndex,(double)mCent9);
    const double Q1yCtr = p_mEpdQ1GrpAveReCtrYWest[mVzBin][grpId]->GetBinContent(binY);

    const double Q1x = Q1VecRaw.X() - Q1xCtr;
    const double Q1y = Q1VecRaw.Y() - Q1yCtr;

    Q1VecReCtr.Set(Q1x, Q1y);
  }

  return Q1VecReCtr;
}

TVector2 StEpdEpManager::getQ1VecGrpAveReCtrFull(int grpId)
{
  TVector2 Q1VecFull(0.0,0.0);
  TVector2 Q1VecEast = getQ1VecGrpAveReCtrEast(grpId);
  TVector2 Q1VecWest = getQ1VecGrpAveReCtrWest(grpId);
  if(Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() >0.0)
  {
    Q1VecFull = Q1VecWest + Q1VecEast;
  }

  return Q1VecFull;
}

double StEpdEpManager::getPsi1GrpRawEast(int grpId)
{
  double Psi1Raw = -999.9;
  TVector2 Q1Vector = getQ1VecGrpRawEast(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Raw = transPsi1(Psi1);
  }

  return Psi1Raw;
}

double StEpdEpManager::getPsi1GrpRawWest(int grpId)
{
  double Psi1Raw = -999.9;
  TVector2 Q1Vector = getQ1VecGrpRawWest(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Raw = transPsi1(Psi1);
  }

  return Psi1Raw;
}

double StEpdEpManager::getPsi1GrpRawFull(int grpId)
{
  double Psi1Raw = -999.9;
  TVector2 Q1Vector = getQ1VecGrpRawFull(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Raw = transPsi1(Psi1);
  }

  return Psi1Raw;
}

double StEpdEpManager::getPsi1GrpWgtEast(int grpId)
{
  double Psi1Wgt = -999.9;
  TVector2 Q1Vector = getQ1VecGrpWgtEast(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Wgt = transPsi1(Psi1);
  }

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1GrpWgtWest(int grpId)
{
  double Psi1Wgt = -999.9;
  TVector2 Q1Vector = getQ1VecGrpWgtWest(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Wgt = transPsi1(Psi1);
  }

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1GrpWgtFull(int grpId)
{
  double Psi1Wgt = -999.9;
  TVector2 Q1Vector = getQ1VecGrpWgtFull(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1Wgt = transPsi1(Psi1);
  }

  return Psi1Wgt;
}

double StEpdEpManager::getPsi1GrpReCtrEast(int grpId)
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector = getQ1VecGrpReCtrEast(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

double StEpdEpManager::getPsi1GrpReCtrWest(int grpId)
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector = getQ1VecGrpReCtrWest(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

double StEpdEpManager::getPsi1GrpReCtrFull(int grpId)
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector = getQ1VecGrpReCtrFull(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

double StEpdEpManager::getPsi1GrpAveReCtrEast(int grpId)
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector = getQ1VecGrpAveReCtrEast(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

double StEpdEpManager::getPsi1GrpAveReCtrWest(int grpId)
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector = getQ1VecGrpAveReCtrWest(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}

double StEpdEpManager::getPsi1GrpAveReCtrFull(int grpId)
{
  double Psi1ReCtr = -999.9;
  TVector2 Q1Vector = getQ1VecGrpAveReCtrFull(grpId);
  if(Q1Vector.Mod() > 0.0)
  {
    double Psi1 = TMath::ATan2(Q1Vector.Y(),Q1Vector.X()); // -pi to pi
    Psi1ReCtr = transPsi1(Psi1);
  }

  return Psi1ReCtr;
}
//---------------------------------------------------------------------------------
// raw EP
void StEpdEpManager::initEpdSubEpSideRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1SideRawEastCent%d",iCent);
    h_mEpdEp1SideRawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideRawWestCent%d",iCent);
    h_mEpdEp1SideRawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideRawFullCent%d",iCent);
    h_mEpdEp1SideRawFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideRawCorrCent%d",iCent);
    h_mEpdEp1SideRawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpSideRaw(double Psi1East, double Psi1West, double Psi1Full)
{
  h_mEpdEp1SideRawEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1SideRawWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1SideRawFull[mCent9]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1SideRawCorr[mCent9]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpSideRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1SideRawEast[iCent]->Write();
    h_mEpdEp1SideRawWest[iCent]->Write();
    h_mEpdEp1SideRawFull[iCent]->Write();
    h_mEpdEp1SideRawCorr[iCent]->Write();
  }
}

void StEpdEpManager::initEpdSubEpGrpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      std::string histName = Form("h_mEpdEp1Grp%dRawEastCent%d",iGrp,iCent);
      h_mEpdEp1GrpRawEast[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dRawWestCent%d",iGrp,iCent);
      h_mEpdEp1GrpRawWest[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dRawFullCent%d",iGrp,iCent);
      h_mEpdEp1GrpRawFull[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dRawCorrCent%d",iGrp,iCent);
      h_mEpdEp1GrpRawCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}
void StEpdEpManager::fillEpdSubEpGrpRaw(double Psi1East, double Psi1West, double Psi1Full, int grpId)
{
  h_mEpdEp1GrpRawEast[mCent9][grpId]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1GrpRawWest[mCent9][grpId]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1GrpRawFull[mCent9][grpId]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1GrpRawCorr[mCent9][grpId]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpGrpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      h_mEpdEp1GrpRawEast[iCent][iGrp]->Write();
      h_mEpdEp1GrpRawWest[iCent][iGrp]->Write();
      h_mEpdEp1GrpRawFull[iCent][iGrp]->Write();
      h_mEpdEp1GrpRawCorr[iCent][iGrp]->Write();
    }
  }
}

// phi weighted EP
void StEpdEpManager::initEpdSubEpSideWgt()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1SideWgtEastCent%d",iCent); // 1st EP
    h_mEpdEp1SideWgtEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideWgtWestCent%d",iCent);
    h_mEpdEp1SideWgtWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideWgtFullCent%d",iCent);
    h_mEpdEp1SideWgtFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideWgtCorrCent%d",iCent);
    h_mEpdEp1SideWgtCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpSideWgt(double Psi1East, double Psi1West, double Psi1Full)
{
  h_mEpdEp1SideWgtEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1SideWgtWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1SideWgtFull[mCent9]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1SideWgtCorr[mCent9]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpSideWgt()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1SideWgtEast[iCent]->Write();
    h_mEpdEp1SideWgtWest[iCent]->Write();
    h_mEpdEp1SideWgtFull[iCent]->Write();
    h_mEpdEp1SideWgtCorr[iCent]->Write();
  }
}

void StEpdEpManager::initEpdSubEpGrpWgt()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      std::string histName = Form("h_mEpdEp1Grp%dWgtEastCent%d",iGrp,iCent);
      h_mEpdEp1GrpWgtEast[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dWgtWestCent%d",iGrp,iCent);
      h_mEpdEp1GrpWgtWest[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dWgtFullCent%d",iGrp,iCent);
      h_mEpdEp1GrpWgtFull[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dWgtCorrCent%d",iGrp,iCent);
      h_mEpdEp1GrpWgtCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StEpdEpManager::fillEpdSubEpGrpWgt(double Psi1East, double Psi1West, double Psi1Full, int grpId)
{
  h_mEpdEp1GrpWgtEast[mCent9][grpId]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1GrpWgtWest[mCent9][grpId]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1GrpWgtFull[mCent9][grpId]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1GrpWgtCorr[mCent9][grpId]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpGrpWgt()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      h_mEpdEp1GrpWgtEast[iCent][iGrp]->Write();
      h_mEpdEp1GrpWgtWest[iCent][iGrp]->Write();
      h_mEpdEp1GrpWgtFull[iCent][iGrp]->Write();
      h_mEpdEp1GrpWgtCorr[iCent][iGrp]->Write();
    }
  }
}

// recenter EP
void StEpdEpManager::initEpdSubEpSideReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1SideReCtrEastCent%d",iCent); // 2nd EP
    h_mEpdEp1SideReCtrEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideReCtrWestCent%d",iCent);
    h_mEpdEp1SideReCtrWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideReCtrFullCent%d",iCent);
    h_mEpdEp1SideReCtrFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideReCtrCorrCent%d",iCent);
    h_mEpdEp1SideReCtrCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpSideReCtr(double Psi1East, double Psi1West, double Psi1Full)
{
  h_mEpdEp1SideReCtrEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1SideReCtrWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1SideReCtrFull[mCent9]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1SideReCtrCorr[mCent9]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpSideReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1SideReCtrEast[iCent]->Write();
    h_mEpdEp1SideReCtrWest[iCent]->Write();
    h_mEpdEp1SideReCtrFull[iCent]->Write();
    h_mEpdEp1SideReCtrCorr[iCent]->Write();
  }
}

void StEpdEpManager::initEpdSubEpGrpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      std::string histName = Form("h_mEpdEp1Grp%dReCtrEastCent%d",iGrp,iCent);
      h_mEpdEp1GrpReCtrEast[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dReCtrWestCent%d",iGrp,iCent);
      h_mEpdEp1GrpReCtrWest[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dReCtrFullCent%d",iGrp,iCent);
      h_mEpdEp1GrpReCtrFull[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dReCtrCorrCent%d",iGrp,iCent);
      h_mEpdEp1GrpReCtrCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StEpdEpManager::fillEpdSubEpGrpReCtr(double Psi1East, double Psi1West, double Psi1Full, int grpId)
{
  h_mEpdEp1GrpReCtrEast[mCent9][grpId]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1GrpReCtrWest[mCent9][grpId]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1GrpReCtrFull[mCent9][grpId]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1GrpReCtrCorr[mCent9][grpId]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpGrpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      h_mEpdEp1GrpReCtrEast[iCent][iGrp]->Write();
      h_mEpdEp1GrpReCtrWest[iCent][iGrp]->Write();
      h_mEpdEp1GrpReCtrFull[iCent][iGrp]->Write();
      h_mEpdEp1GrpReCtrCorr[iCent][iGrp]->Write();
    }
  }
}

// shift EP
void StEpdEpManager::initEpdSubEpSideShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1SideShiftEastCent%d",iCent); // 2nd EP
    h_mEpdEp1SideShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideShiftWestCent%d",iCent);
    h_mEpdEp1SideShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideShiftFullCent%d",iCent);
    h_mEpdEp1SideShiftFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1SideShiftCorrCent%d",iCent);
    h_mEpdEp1SideShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdSubEpSideShift(double Psi1East, double Psi1West, double Psi1Full)
{
  h_mEpdEp1SideShiftEast[mCent9]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1SideShiftWest[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1SideShiftFull[mCent9]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1SideShiftCorr[mCent9]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpSideShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1SideShiftEast[iCent]->Write();
    h_mEpdEp1SideShiftWest[iCent]->Write();
    h_mEpdEp1SideShiftFull[iCent]->Write();
    h_mEpdEp1SideShiftCorr[iCent]->Write();
  }
}

void StEpdEpManager::initEpdSubEpGrpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      std::string histName = Form("h_mEpdEp1Grp%dShiftEastCent%d",iGrp,iCent);
      h_mEpdEp1GrpShiftEast[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dShiftWestCent%d",iGrp,iCent);
      h_mEpdEp1GrpShiftWest[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dShiftFullCent%d",iGrp,iCent);
      h_mEpdEp1GrpShiftFull[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
      histName = Form("h_mEpdEp1Grp%dShiftCorrCent%d",iGrp,iCent);
      h_mEpdEp1GrpShiftCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StEpdEpManager::fillEpdSubEpGrpShift(double Psi1East, double Psi1West, double Psi1Full, int grpId)
{
  h_mEpdEp1GrpShiftEast[mCent9][grpId]->Fill(mRunIndex,Psi1East);
  h_mEpdEp1GrpShiftWest[mCent9][grpId]->Fill(mRunIndex,Psi1West);
  h_mEpdEp1GrpShiftFull[mCent9][grpId]->Fill(mRunIndex,Psi1Full);
  h_mEpdEp1GrpShiftCorr[mCent9][grpId]->Fill(Psi1East,Psi1West);
}

void StEpdEpManager::writeEpdSubEpGrpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      h_mEpdEp1GrpShiftEast[iCent][iGrp]->Write();
      h_mEpdEp1GrpShiftWest[iCent][iGrp]->Write();
      h_mEpdEp1GrpShiftFull[iCent][iGrp]->Write();
      h_mEpdEp1GrpShiftCorr[iCent][iGrp]->Write();
    }
  }
}

// shift Full EP
void StEpdEpManager::initEpdFullEpSideShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1SideShiftFullCorrCent%d",iCent); // 1st EP
    h_mEpdEp1SideShiftFullCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
  }
}

void StEpdEpManager::fillEpdFullEpSideShift(double Psi1FullCorr)
{
  h_mEpdEp1SideShiftFullCorr[mCent9]->Fill(mRunIndex,Psi1FullCorr);
}

void StEpdEpManager::writeEpdFullEpSideShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1SideShiftFullCorr[iCent]->Write();
  }
}

void StEpdEpManager::initEpdFullEpGrpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      std::string histName = Form("h_mEpdEp1Grp%dShiftFullCorrCent%d",iGrp,iCent);
      h_mEpdEp1GrpShiftFullCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StEpdEpManager::fillEpdFullEpGrpShift(double Psi1FullCorr, int grpId)
{
  h_mEpdEp1GrpShiftFullCorr[mCent9][grpId]->Fill(mRunIndex,Psi1FullCorr);
}

void StEpdEpManager::writeEpdFullEpGrpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      h_mEpdEp1GrpShiftFullCorr[iCent][iGrp]->Write();
    }
  }
}
//---------------------------------------------------------------------------------
