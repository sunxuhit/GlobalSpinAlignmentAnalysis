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
#include "StRoot/StPileupUtil/StPileupUtil.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
#include "StRoot/StEventPlaneMaker/StEpdEpManager.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"
// #include "StRoot/StEventPlaneMaker/StMixEpManager.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonTree.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonMaker.h"

ClassImp(StPhiMesonMaker)

//-----------------------------------------------------------------------------
StPhiMesonMaker::StPhiMesonMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int mode, const int beamType, const int flagME) : StMaker(name), mMode(mode), mType(beamType), mFlagME(flagME)
{
  mPicoDstMaker      = picoMaker;
  mPicoDst           = NULL;
  mRefMultCorr       = NULL;
  mPileupUtilFxt3p85 = NULL;

  if(mMode == 0) // phi Meson TTree Production
  {
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
  // mMixEpManager = new StMixEpManager(mType); // initialize Mix EP Manager
  mAnaCut       = new StAnalysisCut(mType);
  mAnaUtils     = new StAnalysisUtils(mType);
  mAnaUtils->initRunIndex(); // initialize std::map for run index

  if(!mRefMultCorr)
  {
    if( mAnaCut->isIsobar() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr_Isobar();
    if( mAnaCut->isFxt3p85GeV_2018() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  if(mAnaCut->isFxt3p85GeV_2018() && !mPileupUtilFxt3p85)
  {
    mPileupUtilFxt3p85 = new StPileupUtil();
    mPileupUtilFxt3p85->init();
  }

  if(mAnaCut->isIsobar())
  {
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->readZdcReCtr();
    mZdcEpManager->readZdcShift();
    mZdcEpManager->readZdcShiftFull();
    // mZdcEpManager->readZdcResolution();
    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->readEpdSideReCtr();
    mEpdEpManager->readEpdSideShift();
    mEpdEpManager->readEpdSideShiftFull();
    // mEpdEpManager->readEpdSideResolution();
  }
  if(mAnaCut->isFxt3p85GeV_2018())
  {
    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->readEpdGrpReCtr();
    mEpdEpManager->readEpdGrpShift();
    // mMixEpManager->readMixEpRes(); // Mix
  }
  mTpcEpManager->readTpcReCtr(); // TPC
  mTpcEpManager->readTpcShift();
  // mTpcEpManager->readTpcResolution();

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
    int refMult            = mPicoEvent->refMult();
    const TVector3 primVtx = mPicoEvent->primaryVertex();
    const double vz        = mPicoEvent->primaryVertex().z();
    const double zdcX      = mPicoEvent->ZDCx();
    const unsigned short nBTofMatch = mPicoEvent->nBTOFMatch();     // get number of tof match points
    const unsigned int nTracks      = mPicoDst->numberOfTracks();   // get number of tracks
    const unsigned int nEpdHits     = mPicoDst->numberOfEpdHits();  // get number of EPD hits

    const int vzBin  = mAnaUtils->getVzBin(vz); // 0 for -vz || 1 for +vz
    const int runIdx = mAnaUtils->findRunIndex(runId); // find run index for a specific run
    // cout << "runId = " << runId << ", runIdx = " << runIdx << endl;
    if(runIdx < 0)
    {
      LOG_ERROR << "Could not find this run Index from StAnalysisUtils! Skip!" << endm;
      return kStErr;
    }

    // StRefMultCorr Cut & centrality
    if(!mRefMultCorr)
    {
      LOG_WARN << " No mRefMultCorr! Skip! " << endl;
      return kStErr;
    }
    mRefMultCorr->init(runId);
    mRefMultCorr->initEvent(refMult,vz,zdcX);
    int cent9     = mRefMultCorr->getCentralityBin9(); // get Centrality9
    int cent16    = mRefMultCorr->getCentralityBin16(); // get Centrality16
    double refWgt = mRefMultCorr->getWeight(); // get weight

    if(mAnaCut->isFxt3p85GeV_2018()) // only for Fxt3p85GeV_2018
    { 
      mPileupUtilFxt3p85->initEvent(mPicoDst);
      refMult = mPileupUtilFxt3p85->get_refMultPrim();
      cent9   = mPileupUtilFxt3p85->get_centrality9();
      cent16  = mPileupUtilFxt3p85->get_centrality16();
      refWgt  = mPileupUtilFxt3p85->get_centralityWeight();
    }

    // reject bad runs
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

    bool isPileUpStAnalysisCut = mAnaCut->isPileUpEvent((double)refMult, (double)nBTofMatch,vz); // alway return false for Isobar
    bool isPileUpStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut((double)refMult, (double)nBTofMatch,vz); // valid for Isobar
    if(mAnaCut->isFxt3p85GeV_2018()) isPileUpStRefMultCorr = mPileupUtilFxt3p85->isPileupEPD(); // valid only for Fxt3p85GeV_2018
    bool isPileUpEvt = isPileUpStAnalysisCut || isPileUpStRefMultCorr;
    // cout << "isPileUpEvt = " << isPileUpEvt << ", isPileUpStAnalysisCut = " << isPileUpStAnalysisCut << ", isPileUpStRefMultCorr = " << isPileUpStRefMultCorr << endl;

    if(!isPileUpEvt && mAnaCut->isGoodCent9(cent9) && mAnaCut->passEventCut(mPicoEvent))
    { // apply Event Cuts for anlaysis 
      mZdcEpManager->initZdcEpManager(cent9,runIdx,vzBin); // initialize ZDC EP Manager
      mEpdEpManager->initEpdEpManager(cent9,runIdx,vzBin); // initialize EPD EP Manager
      mTpcEpManager->initTpcEpManager(cent9,runIdx,vzBin); // initialize TPC EP Manager
      // mMixEpManager->initMixEpManager(cent9,runIdx,vzBin); // initialize Mix EP Manager

      // Calculate QVector
      if( mAnaCut->isIsobar() )
      {
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
	    mEpdEpManager->addHitSideReCtrEast(picoEpdHit, primVtx);
	  }
	  if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	  {
	    mEpdEpManager->addHitSideReCtrWest(picoEpdHit, primVtx);
	  }
	}
      }

      for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack) // calculate recentered Q1Vector, Q2Vector and Q3Vector from TPC
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

      if( mAnaCut->isFxt3p85GeV_2018() )
      {
	for(unsigned int iEpdHit = 0; iEpdHit < nEpdHits; ++iEpdHit) // EPD QVector
	{ // calculate phi & eta weighted (if avaliable) and recentered Q1Vector from EPD
	  StPicoEpdHit *picoEpdHit = (StPicoEpdHit*)mPicoDst->epdHit(iEpdHit); 
	  if(!picoEpdHit) continue;

	  if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	  {
	    mEpdEpManager->addHitGrpReCtrEast(picoEpdHit, primVtx);
	  }
	  if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	  {
	    mEpdEpManager->addHitGrpReCtrWest(picoEpdHit, primVtx);
	  }
	}
      }

      if(mMode == 0) // fill phi meson TTree
      {
	if( mAnaCut->isIsobar() )
	{
	  int flagZdcEp = 0; // ZDC EP
	  const TVector2 vQ1ZdcEast = mZdcEpManager->getQ1VecEast(5); // get Shift Corrected Q1Vector from ZDC
	  const TVector2 vQ1ZdcWest = mZdcEpManager->getQ1VecWest(5);
	  const TVector2 vQ1ZdcFull = mZdcEpManager->getQ1VecFull(vQ1ZdcEast,vQ1ZdcWest,5); // TVector2 vQ1ZdcFull = vQ1ZdcWest-vQ1ZdcEast;
	  if( mAnaCut->passQVecZdc(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull) ) // ZDC EP
	  {
	    flagZdcEp = 1; // flag for ZDC EP
	  }

	  int flagEpdSideEp = 0; // EPD Side EP
	  const TVector2 vQ1EpdSideEast = mEpdEpManager->getQ1VecSideShiftEast(); // get Shift Corrected Q1Vector from EPD Side
	  const TVector2 vQ1EpdSideWest = mEpdEpManager->getQ1VecSideShiftWest();
	  const TVector2 vQ1EpdSideFull = mEpdEpManager->getQ1VecSideShiftFull();
	  const TVector2 vQ1EpdSideFullCorr = mEpdEpManager->getQ1VecSideShiftFullCorr();
	  if( mAnaCut->passQVecEpdSide(vQ1EpdSideEast,vQ1EpdSideWest,vQ1EpdSideFull) ) // EPD EP
	  {
	    flagEpdSideEp = 1;
	  }

	  const int flagEpdGrp0Ep = -1; // EPD Grp0 EP & NOT used in Isobar
	  const TVector2 vQ1EpdGrp0East(0.0,0.0); // get Shift Corrected Q1Vector from EPD Grp0
	  const TVector2 vQ1EpdGrp0West(0.0,0.0);
	  const TVector2 vQ1EpdGrp0Full(0.0,0.0);
	  const TVector2 vQ1EpdGrp0FullCorr(0.0,0.0);
	  const int flagEpdGrp1Ep = -1; // EPD Grp0 EP & NOT used in Isobar
	  const TVector2 vQ1EpdGrp1East(0.0,0.0); // get Shift Corrected Q1Vector from EPD Grp1
	  const TVector2 vQ1EpdGrp1West(0.0,0.0);
	  const TVector2 vQ1EpdGrp1Full(0.0,0.0);
	  const TVector2 vQ1EpdGrp1FullCorr(0.0,0.0);

	  int flagTpcEp = 0; // TPC EP
	  const TVector2 vQ1TpcEast = mTpcEpManager->getQ1VecReCtrEast();
	  const TVector2 vQ1TpcWest = mTpcEpManager->getQ1VecReCtrWest();
	  const TVector2 vQ2TpcEast = mTpcEpManager->getQ2VecReCtrEast();
	  const TVector2 vQ2TpcWest = mTpcEpManager->getQ2VecReCtrWest();
	  const TVector2 vQ3TpcEast = mTpcEpManager->getQ3VecReCtrEast();
	  const TVector2 vQ3TpcWest = mTpcEpManager->getQ3VecReCtrWest();
	  const int numTrkReCtrEast = mTpcEpManager->getNumTrkReCtrEast();
	  const int numTrkReCtrWest = mTpcEpManager->getNumTrkReCtrWest();
	  if(mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest))
	  {
	    flagTpcEp = 1;
	  }

	  const TVector2 vQ2TpcFull = mTpcEpManager->getQ2VecReCtrFull();
	  const double Psi2ShiftFull  = mTpcEpManager->getPsi2ShiftFull(vQ2TpcFull); // Psi2 from Full TPC EP & Used for Event Mixing

	  mPhiMesonTree->clearEvtInfo();
	  mPhiMesonTree->setEvtInfo(runIdx, cent9, cent16, refWgt, vz, Psi2ShiftFull); // Use Psi2ShiftFull for Event Mixing
	  mPhiMesonTree->setZdcQ1Flag(flagZdcEp); // ZDC EP
	  mPhiMesonTree->setZdcQ1Vec(vQ1ZdcEast, vQ1ZdcWest, vQ1ZdcFull); // Shift Corrected ZDC Q1Vector: East & West & Full
	  mPhiMesonTree->setEpdQ1SideFlag(flagEpdSideEp); // EPD Side EP
	  mPhiMesonTree->setEpdQ1SideVec(vQ1EpdSideEast, vQ1EpdSideWest, vQ1EpdSideFullCorr); // Shift Corrected EPD Side Q1Vector: East & West & Full
	  mPhiMesonTree->setEpdQ1Grp0Flag(flagEpdGrp0Ep); // EPD Grp0 EP & NOT Used in Isobar
	  mPhiMesonTree->setEpdQ1Grp0Vec(vQ1EpdGrp0East, vQ1EpdGrp0West, vQ1EpdGrp0FullCorr); // Shift Corrected EPD Grp0 Q1Vector: East & West & Full
	  mPhiMesonTree->setEpdQ1Grp1Flag(flagEpdGrp1Ep); // EPD Grp1 EP & NOT Used in Isobar
	  mPhiMesonTree->setEpdQ1Grp1Vec(vQ1EpdGrp1East, vQ1EpdGrp1West, vQ1EpdGrp1FullCorr); // Shift Corrected EPD Grp1 Q1Vector: East & West & Full
	  mPhiMesonTree->setTpcQFlag(flagTpcEp); // TPC EP
	  mPhiMesonTree->setTpcQ1Vec(vQ1TpcEast, vQ1TpcWest); // ReCenter Corrected TPC Q1Vector: East & West
	  mPhiMesonTree->setTpcQ2Vec(vQ2TpcEast, vQ2TpcWest); // ReCenter Corrected TPC Q2Vector: East & West
	  mPhiMesonTree->setTpcQ3Vec(vQ3TpcEast, vQ3TpcWest); // ReCenter Corrected TPC Q3Vector: East & West
	  mPhiMesonTree->setNumTrks(numTrkReCtrEast, numTrkReCtrWest); // Number of Tracks used in ReCenter Correction: East & West
	  mPhiMesonTree->fillPhiTree(mPicoDst, mFlagME);
	}
	if( mAnaCut->isFxt3p85GeV_2018() )
	{
	  const int flagZdcEp = -1; // ZDC EP & NOT Used in FXT
	  const TVector2 vQ1ZdcEast(0.0,0.0);
	  const TVector2 vQ1ZdcWest(0.0,0.0);
	  const TVector2 vQ1ZdcFull(0.0,0.0);

	  const int flagEpdSideEp = -1; // EPD Side EP & NOT Used in FXT
	  const TVector2 vQ1EpdSideEast(0.0,0.0);
	  const TVector2 vQ1EpdSideWest(0.0,0.0);
	  const TVector2 vQ1EpdSideFull(0.0,0.0);
	  const TVector2 vQ1EpdSideFullCorr(0.0,0.0);

	  int flagEpdGrp0Ep = 0; // EPD Grp0 EP
	  const TVector2 vQ1EpdGrp0East = mEpdEpManager->getQ1VecGrpShiftEast(0); // get Shift Corrected Q1Vector from EPD Grp0
	  const TVector2 vQ1EpdGrp0West = mEpdEpManager->getQ1VecGrpShiftWest(0);
	  const TVector2 vQ1EpdGrp0Full = mEpdEpManager->getQ1VecGrpShiftFull(0);
	  // const TVector2 vQ1EpdGrp0FullCorr = mEpdEpManager->getQ1VecGrpShiftFullCorr(0);
	  if( mAnaCut->passQVecEpdGrp(vQ1EpdGrp0East,vQ1EpdGrp0West,vQ1EpdGrp0Full,0) ) // EPD EP Grp0
	  {
	    flagEpdGrp0Ep = 1;
	  }

	  int flagEpdGrp1Ep = 0; // EPD Grp1 EP
	  const TVector2 vQ1EpdGrp1East = mEpdEpManager->getQ1VecGrpShiftEast(1); // get Shift Corrected Q1Vector from EPD Grp0
	  const TVector2 vQ1EpdGrp1West = mEpdEpManager->getQ1VecGrpShiftWest(1);
	  const TVector2 vQ1EpdGrp1Full = mEpdEpManager->getQ1VecGrpShiftFull(1);
	  // const TVector2 vQ1EpdGrp1FullCorr = mEpdEpManager->getQ1VecGrpShiftFullCorr(1);
	  if( mAnaCut->passQVecEpdGrp(vQ1EpdGrp1East,vQ1EpdGrp1West,vQ1EpdGrp1Full,1) ) // EPD EP Grp1
	  {
	    flagEpdGrp1Ep = 1;
	  }

	  int flagTpcEp = 0; // TPC EP
	  const TVector2 vQ1TpcEast = mTpcEpManager->getQ1VecReCtrEast();
	  const TVector2 vQ1TpcWest = mTpcEpManager->getQ1VecReCtrWest();
	  const TVector2 vQ2TpcEast = mTpcEpManager->getQ2VecReCtrEast();
	  const TVector2 vQ2TpcWest = mTpcEpManager->getQ2VecReCtrWest();
	  const TVector2 vQ3TpcEast = mTpcEpManager->getQ3VecReCtrEast();
	  const TVector2 vQ3TpcWest = mTpcEpManager->getQ3VecReCtrWest();
	  const int numTrkReCtrEast = mTpcEpManager->getNumTrkReCtrEast();
	  const int numTrkReCtrWest = mTpcEpManager->getNumTrkReCtrWest();
	  if(mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest))
	  {
	    flagTpcEp = 1;
	  }

	  const double Psi1EpdGrp0 = mEpdEpManager->getPsi1GrpShiftEast(0); // Psi1 from EPD Grp 0 East & Used for Event Mixing

	  mPhiMesonTree->clearEvtInfo();
	  mPhiMesonTree->setEvtInfo(runIdx, refMult, cent9, cent16, refWgt, vz, Psi1EpdGrp0); // Use Psi1EpdGrp0 for Event Mixing
	  mPhiMesonTree->setZdcQ1Flag(flagZdcEp); // ZDC EP
	  mPhiMesonTree->setZdcQ1Vec(vQ1ZdcEast, vQ1ZdcWest, vQ1ZdcFull); // Shift Corrected ZDC Q1Vector: East & West & Full
	  mPhiMesonTree->setEpdQ1SideFlag(flagEpdSideEp); // EPD Side EP
	  mPhiMesonTree->setEpdQ1SideVec(vQ1EpdSideEast, vQ1EpdSideWest, vQ1EpdSideFullCorr); // Shift Corrected EPD Side Q1Vector: East & West & Full
	  mPhiMesonTree->setEpdQ1Grp0Flag(flagEpdGrp0Ep); // EPD Grp0 EP & NOT Used in Isobar
	  mPhiMesonTree->setEpdQ1Grp0Vec(vQ1EpdGrp0East, vQ1EpdGrp0West, vQ1EpdGrp0Full); // Shift Corrected EPD Grp0 Q1Vector: East & West & Full
	  mPhiMesonTree->setEpdQ1Grp1Flag(flagEpdGrp1Ep); // EPD Grp1 EP & NOT Used in Isobar
	  mPhiMesonTree->setEpdQ1Grp1Vec(vQ1EpdGrp1East, vQ1EpdGrp1West, vQ1EpdGrp1Full); // Shift Corrected EPD Grp1 Q1Vector: East & West & Full
	  mPhiMesonTree->setTpcQFlag(flagTpcEp); // TPC EP
	  mPhiMesonTree->setTpcQ1Vec(vQ1TpcEast, vQ1TpcWest); // ReCenter Corrected TPC Q1Vector: East & West
	  mPhiMesonTree->setTpcQ2Vec(vQ2TpcEast, vQ2TpcWest); // ReCenter Corrected TPC Q2Vector: East & West
	  mPhiMesonTree->setTpcQ3Vec(vQ3TpcEast, vQ3TpcWest); // ReCenter Corrected TPC Q3Vector: East & West
	  mPhiMesonTree->setNumTrks(numTrkReCtrEast, numTrkReCtrWest); // Number of Tracks used in ReCenter Correction: East & West
	  mPhiMesonTree->fillPhiTree(mPicoDst, mFlagME);
	}
      }

      mZdcEpManager->clearZdcEpManager();
      mEpdEpManager->clearEpdEpManager();
      mTpcEpManager->clearTpcEpManager();
      // mMixEpManager->clearMixEpManager();
    }
  }

  return kStOK;
}

