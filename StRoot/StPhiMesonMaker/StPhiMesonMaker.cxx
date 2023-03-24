#include <algorithm>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

#include "StMessMgr.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
#include "StRoot/StEventPlaneMaker/StEpdEpManager.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonTree.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonMaker.h"

ClassImp(StPhiMesonMaker)

//-----------------------------------------------------------------------------
StPhiMesonMaker::StPhiMesonMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int mode, const int beamType, const int flagME) : StMaker(name), mMode(mode), mType(beamType), mFlagME(flagME)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;

  if(mMode == 0)
  { // Test phi Meson Reconstruction
    str_mOutPutPhiTree = Form("./file_RecoPhi%s_%s_%s.root",str_mFlagME[mFlagME].c_str(),globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
}

//----------------------------------------------------------------------------- 
StPhiMesonMaker::~StPhiMesonMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StPhiMesonMaker::Init() 
{
  mZdcEpManager = new StZdcEpManager(mType); // initialize ZDC EP Manager
  mEpdEpManager = new StEpdEpManager(mType); // initialize EPD EP Manager
  mTpcEpManager = new StTpcEpManager(mType); // initialize TPC EP Manager
  mAnaCut       = new StAnalysisCut(mType);
  mAnaUtils     = new StAnalysisUtils(mType);
  mAnaUtils->initRunIndex(); // initialize std::map for run index

  if(!mRefMultCorr)
  {
    if( mAnaCut->isIsobar() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr_Isobar();
  }

  // initialize EP manager
  mZdcEpManager->readZdcGain(); // ZDC
  mZdcEpManager->readZdcReCtr();
  mZdcEpManager->readZdcShift();
  mZdcEpManager->readZdcShiftFull();
  mZdcEpManager->readZdcResolution();

  mEpdEpManager->readEpdPhiWgt(); // EPD
  mEpdEpManager->readEpdReCtr();
  mEpdEpManager->readEpdShift();
  mEpdEpManager->readEpdShiftFull();
  mEpdEpManager->readEpdResolution();

  mTpcEpManager->readTpcReCtr(); // TPC
  mTpcEpManager->readTpcShift();
  mTpcEpManager->readTpcResolution();

  if(mMode == 0)
  { // Test phi Meson Reconstruction
    file_mOutPutPhiTree = new TFile(str_mOutPutPhiTree.c_str(),"RECREATE");
    mPhiMesonTree = new StPhiMesonTree(mType);
    mPhiMesonTree->initPhiTree();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StPhiMesonMaker::Finish() 
{
  if(mMode == 0)
  {
    if(str_mOutPutPhiTree != "")
    {
      file_mOutPutPhiTree->cd();
      mPhiMesonTree->writePhiTree();
      file_mOutPutPhiTree->Close();
    }
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StPhiMesonMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
int StPhiMesonMaker::Make() 
{
  if(!mPicoDstMaker) 
  {
    LOG_ERROR << " No PicoDstMaker! Skip! " << endm;
    return kStErr;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) 
  {
    LOG_ERROR << " No PicoDst! Skip! " << endm;
    return kStErr;
  }

  mPicoEvent = (StPicoEvent*)mPicoDst->event();
  if(!mPicoEvent)
  {
    LOG_ERROR << "Error opening picoDst Event, skip!" << endm;
    return kStErr;
  }

  // MinBias trigger
  if( mAnaCut->isMinBias(mPicoEvent) )
  {
    // Event Information
    const int runId        = mPicoEvent->runId();
    const int refMult      = mPicoEvent->refMult();
    // const int grefMult     = mPicoEvent->grefMult();
    const TVector3 primVtx = mPicoEvent->primaryVertex();
    // const double vx        = mPicoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
    // const double vy        = mPicoEvent->primaryVertex().y();
    const double vz        = mPicoEvent->primaryVertex().z();
    // const double vzVpd     = mPicoEvent->vzVpd();
    const double zdcX      = mPicoEvent->ZDCx();
    // const unsigned short nBTofHits  = mPicoEvent->btofTrayMultiplicity();
    // const unsigned int nBTofHits    = mPicoDst->numberOfBTofHits(); // get number of tof hits
    const unsigned short nBTofMatch = mPicoEvent->nBTOFMatch();     // get number of tof match points
    const unsigned int nTracks      = mPicoDst->numberOfTracks();   // get number of tracks
    const unsigned int nEpdHits     = mPicoDst->numberOfEpdHits();  // get number of EPD hits

    // StRefMultCorr Cut & centrality
    if(!mRefMultCorr)
    {
      LOG_WARN << " No mRefMultCorr! Skip! " << endl;
      return kStErr;
    }

    mRefMultCorr->init(runId);
    mRefMultCorr->initEvent(refMult,vz,zdcX);
    if(mRefMultCorr->isBadRun(runId))
    {
      LOG_ERROR << "Bad Run from StRefMultCorr! Skip!" << endm;
      return kStErr;
    }
    if(mAnaUtils->isBadRun(runId))
    {
      LOG_ERROR << "Bad Run from StAnalysisUtils! Skip!" << endm;
      return kStErr;
    }

    const int cent9       = mRefMultCorr->getCentralityBin9(); // get Centrality9
    const int cent16      = mRefMultCorr->getCentralityBin16(); // get Centrality16
    const double reweight = mRefMultCorr->getWeight(); // get weight
    const int runIndex    = mAnaUtils->findRunIndex(runId); // find run index for a specific run
    // const int triggerBin  = mAnaUtils->getTriggerBin(mPicoEvent);
    const int vzBin       = mAnaUtils->getVzBin(vz); // 0 for -vz || 1 for +vz
    // cout << "runId = " << runId << ", runIndex = " << runIndex << endl;
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StAnalysisUtils! Skip!" << endm;
      return kStErr;
    }

    bool isPileUpEventStAnalysisCut = mAnaCut->isPileUpEvent((double)refMult, (double)nBTofMatch,vz); // alway return false for Isobar
    bool isPileUpEventStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut((double)refMult, (double)nBTofMatch,vz); // valid for Isobar
    bool isPileUpEvent = isPileUpEventStAnalysisCut || isPileUpEventStRefMultCorr;
    // cout << "isPileUpEvent = " << isPileUpEvent << ", isPileUpEventStAnalysisCut = " << isPileUpEventStAnalysisCut << ", isPileUpEventStRefMultCorr = " << isPileUpEventStRefMultCorr << endl;

    if(!isPileUpEvent && mAnaCut->isGoodCent9(cent9) && mAnaCut->passEventCut(mPicoEvent))
    { // apply Event Cuts for anlaysis 
      mZdcEpManager->initZdcEpManager(cent9,runIndex,vzBin); // initialize ZDC EP Manager
      mEpdEpManager->initEpdEpManager(cent9,runIndex,vzBin); // initialize EPD EP Manager
      mTpcEpManager->initTpcEpManager(cent9,runIndex,vzBin); // initialize TPC EP Manager

      for(int iSlat = 0; iSlat < 8; ++iSlat) // calculate Q1Vector from ZDC
      {
	mZdcEpManager->setZdcSmdGainCorr(0,0,iSlat,mPicoEvent->ZdcSmdEastVertical(iSlat));
	mZdcEpManager->setZdcSmdGainCorr(0,1,iSlat,mPicoEvent->ZdcSmdEastHorizontal(iSlat));
	mZdcEpManager->setZdcSmdGainCorr(1,0,iSlat,mPicoEvent->ZdcSmdWestVertical(iSlat));
	mZdcEpManager->setZdcSmdGainCorr(1,1,iSlat,mPicoEvent->ZdcSmdWestHorizontal(iSlat));
      }

      for(unsigned int iEpdHit = 0; iEpdHit < nEpdHits; ++iEpdHit) // calculate recentered Q1Vector from EPD
      {
	StPicoEpdHit *picoEpdHit = (StPicoEpdHit*)mPicoDst->epdHit(iEpdHit); 
	if(!picoEpdHit) continue;

	if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	{
	  mEpdEpManager->addHitReCtrEast(picoEpdHit, primVtx);
	}
	if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	{
	  mEpdEpManager->addHitReCtrWest(picoEpdHit, primVtx);
	}
      }

      for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack) // calculate recentered Q2Vector and Q3Vector from TPC
      { // track loop
	StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	if(!picoTrack) continue;

	if( mAnaCut->passTrkTpcEpEast(picoTrack, primVtx) ) // negative eta
	{
	  mTpcEpManager->addTrackReCtrEast(picoTrack);
	}
	if( mAnaCut->passTrkTpcEpWest(picoTrack, primVtx) ) // positive eta
	{
	  mTpcEpManager->addTrackReCtrWest(picoTrack);
	}
	if( mAnaCut->passTrkTpcEpFull(picoTrack, primVtx) )
	{
	  mTpcEpManager->addTrackReCtrFull(picoTrack);
	}
      }

      if(mMode == 0) // fill phi meson TTree
      {
	const TVector2 vQ1ZdcEast = mZdcEpManager->getQ1VecEast(5); // get Shift Corrected Q1Vector from ZDC
	const TVector2 vQ1ZdcWest = mZdcEpManager->getQ1VecWest(5);
	const TVector2 vQ1ZdcFull = mZdcEpManager->getQ1VecFull(vQ1ZdcEast,vQ1ZdcWest,5); // TVector2 vQ1ZdcFull = vQ1ZdcWest-vQ1ZdcEast;
	int flagZdcEp = 0;
	if( vQ1ZdcEast.Mod() > 0.0 && vQ1ZdcWest.Mod() > 0.0 && vQ1ZdcFull.Mod() > 0.0 ) // ZDC EP
	{
	  flagZdcEp = 1; // flag for ZDC EP
	}

	const TVector2 vQ1EpdEast = mEpdEpManager->getQ1VecSideShiftEast(); // get Shift Corrected Q1Vector from EPD
	const TVector2 vQ1EpdWest = mEpdEpManager->getQ1VecSideShiftWest();
	const TVector2 vQ1EpdFull = mEpdEpManager->getQ1VecSideShiftFull();
	const TVector2 vQ1EpdFullCorr = mEpdEpManager->getQ1VecSideShiftFullCorr();
	int flagEpdEp = 0;
	if( vQ1EpdEast.Mod() > 0.0 && vQ1EpdWest.Mod() > 0.0 && vQ1EpdFull.Mod() > 0.0 ) // EPD EP
	{
	  flagEpdEp = 1;
	}

	const int numTrkReCtrEast = mTpcEpManager->getNumTrkReCtrEast(); // TPC EP
	const int numTrkReCtrWest = mTpcEpManager->getNumTrkReCtrWest();
	int flagTpcEp = 0;
	if(mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest))
	{ // fill phi meson TTree only if TPC has avaliable EP
	  flagTpcEp = 1;
	  const TVector2 vQ2TpcEast = mTpcEpManager->getQ2VecReCtrEast();
	  const TVector2 vQ2TpcWest = mTpcEpManager->getQ2VecReCtrWest();
	  const TVector2 vQ3TpcEast = mTpcEpManager->getQ3VecReCtrEast();
	  const TVector2 vQ3TpcWest = mTpcEpManager->getQ3VecReCtrWest();

	  const TVector2 vQ2TpcFull = mTpcEpManager->getQ2VecReCtrFull();
	  const double Psi2ShiftFull  = mTpcEpManager->getPsi2ShiftFull(vQ2TpcFull);

	  mPhiMesonTree->clearEvtInfo();
	  mPhiMesonTree->setEvtInfo(cent9, cent16, reweight, vz, Psi2ShiftFull);
	  mPhiMesonTree->setZdcQ1Vec(flagZdcEp, vQ1ZdcEast, vQ1ZdcWest, vQ1ZdcFull); // Shift Corrected ZDC Q1Vector: East & West & Full
	  mPhiMesonTree->setEpdQ1Vec(flagEpdEp, vQ1EpdEast, vQ1EpdWest, vQ1EpdFullCorr); // Shift Corrected EPD Q1Vector: East & West & Full
	  mPhiMesonTree->setTpcQVec(flagTpcEp, vQ2TpcEast, vQ2TpcWest, vQ3TpcEast, vQ3TpcWest); // ReCenter Corrected TPC Q2Vector & Q3Vector: East & West
	  mPhiMesonTree->setNumTrks(numTrkReCtrEast, numTrkReCtrWest); // Number of Tracks used in ReCenter Correction: East & West
	  mPhiMesonTree->fillPhiTree(mPicoDst, mFlagME);
	}
      }
      mZdcEpManager->clearZdcEpManager();
      mEpdEpManager->clearEpdEpManager();
      mTpcEpManager->clearTpcEpManager();
    }
  }

  return kStOK;
}

