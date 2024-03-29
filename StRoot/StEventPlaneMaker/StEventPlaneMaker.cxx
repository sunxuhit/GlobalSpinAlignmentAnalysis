#include <algorithm>

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
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
#include "StRoot/StEventPlaneMaker/StEpdEpManager.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"
#include "StRoot/StEventPlaneMaker/StMixEpManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneMaker.h"

ClassImp(StEventPlaneMaker)

//-----------------------------------------------------------------------------
StEventPlaneMaker::StEventPlaneMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int mode, const int beamType) : StMaker(name), mMode(mode), mType(beamType)
{
  mPicoDstMaker      = picoMaker;
  mPicoDst           = NULL;
  mRefMultCorr       = NULL;
  mPileupUtilFxt3p85 = NULL;

  if(mMode == 0) // fill Gain Correction for ZDC & phi Weight for EPD
  {
    str_mOutPutGainCorr = Form("./file_GainCorr_%s_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 1) // fill Re-Center Correction for ZDC & EPD & TPC Sub EP
  {
    str_mOutPutReCenterPar = Form("./file_ReCenterPar_%s_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 2) // fill Shift Correction for ZDC & EPD & TPC Sub EP
  {
    str_mOutPutShiftPar = Form("./file_ShiftPar_%s_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 3) // fill Shift Correction for ZDC & EPD Full EP
  {
    str_mOutPutShiftPar = Form("./file_ShiftParFull_%s_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 4) // fill Event Plane Resolution for ZDC & EPD & TPC Sub EP
  {
    str_mOutPutResolution = Form("./file_EpResolution_%s_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 5) // fill Charged Hadron v1 for ZDC and EPD & v2 and v3 for TPC
  {
    str_mOutPutFlow = Form("./file_ChargedFlow_%s_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
}

StEventPlaneMaker::~StEventPlaneMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Init() 
{
  mZdcEpManager = new StZdcEpManager(mType); // initialize ZDC EP Manager
  mEpdEpManager = new StEpdEpManager(mType); // initialize EPD EP Manager
  mTpcEpManager = new StTpcEpManager(mType); // initialize TPC EP Manager
  mMixEpManager = new StMixEpManager(mType); // initialize Mix EP Manager
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

  if(mMode == 0)
  { // fill Gain Correction Factors for ZDC and phi Weight for EPD
    file_mOutPutGainCorr = new TFile(str_mOutPutGainCorr.c_str(),"RECREATE");
    if(mAnaCut->isIsobar()) 
    {
      mZdcEpManager->initZdcGain(); // ZDC
      mEpdEpManager->initEpdPhiWgt(); // EPD
      mEpdEpManager->initEpdSubEpSideRaw();
    }
    if(mAnaCut->isFxt3p85GeV_2018()) 
    {
      mEpdEpManager->initEpdPhiWgt(); // EPD
      mEpdEpManager->initEpdSubEpGrpRaw();
    }
  }
  if(mMode == 1)
  { // fill ReCenter Correction Parameters for ZDC & EPD & TPC Sub EP
    file_mOutPutReCenterPar = new TFile(str_mOutPutReCenterPar.c_str(),"RECREATE");
    if(mAnaCut->isIsobar()) 
    {
      mZdcEpManager->readZdcGain(); // ZDC
      mZdcEpManager->initZdcReCtr();
      mZdcEpManager->initZdcSubEpRaw();
      mEpdEpManager->readEpdPhiWgt(); // EPD
      mEpdEpManager->initEpdSideReCtr();
      mEpdEpManager->initEpdSubEpSideWgt();
    }
    if(mAnaCut->isFxt3p85GeV_2018()) 
    {
      // mEpdEpManager->readEpdPhiWgt(); // EPD
      mEpdEpManager->initEpdGrpReCtr();
      mEpdEpManager->initEpdSubEpGrpWgt();
      mMixEpManager->initMixSubEpRaw(); // Mix
    }
    mTpcEpManager->initTpcReCtr(); // TPC
    mTpcEpManager->initTpcSubEpRaw();
  }
  if(mMode == 2)
  { // fill Shift Correction Parameters for ZDC & EPD & TPC Sub EP
    file_mOutPutShiftPar = new TFile(str_mOutPutShiftPar.c_str(),"RECREATE");
    if(mAnaCut->isIsobar()) 
    {
      mZdcEpManager->readZdcGain(); // ZDC
      mZdcEpManager->readZdcReCtr();
      mZdcEpManager->initZdcShift();
      mZdcEpManager->initZdcSubEpReCtr();
      mEpdEpManager->readEpdPhiWgt(); // EPD
      mEpdEpManager->readEpdSideReCtr();
      mEpdEpManager->initEpdSideShift();
      mEpdEpManager->initEpdSubEpSideReCtr();
    }
    if(mAnaCut->isFxt3p85GeV_2018()) 
    {
      // mEpdEpManager->readEpdPhiWgt(); // EPD
      mEpdEpManager->readEpdGrpReCtr();
      mEpdEpManager->initEpdGrpShift();
      mEpdEpManager->initEpdSubEpGrpReCtr();
      mMixEpManager->initMixSubEpReCtr(); // Mix
    }
    mTpcEpManager->readTpcReCtr(); // TPC
    mTpcEpManager->initTpcShift();
    mTpcEpManager->initTpcSubEpReCtr();
  }
  if(mMode == 3)
  { // fill Shift Correction Parameters for ZDC & EPD Full EP
    file_mOutPutShiftPar = new TFile(str_mOutPutShiftPar.c_str(),"RECREATE");
    if(mAnaCut->isIsobar()) 
    {
      mZdcEpManager->readZdcGain(); // ZDC
      mZdcEpManager->readZdcReCtr();
      mZdcEpManager->readZdcShift();
      mZdcEpManager->initZdcShiftFull();
      mZdcEpManager->initZdcSubEpShift();
      mEpdEpManager->readEpdPhiWgt(); // EPD
      mEpdEpManager->readEpdSideReCtr();
      mEpdEpManager->readEpdSideShift();
      mEpdEpManager->initEpdSideShiftFull();
      mEpdEpManager->initEpdSubEpSideShift();
    }
  }
  if(mMode == 4)
  { // fill Event Plane Resolution for ZDC & EPD & TPC Sub EP
    file_mOutPutResolution = new TFile(str_mOutPutResolution.c_str(),"RECREATE");
    if(mAnaCut->isIsobar()) 
    {
      mZdcEpManager->readZdcGain(); // ZDC
      mZdcEpManager->readZdcReCtr();
      mZdcEpManager->readZdcShift();
      mZdcEpManager->readZdcShiftFull();
      mZdcEpManager->initZdcResolution();
      mZdcEpManager->initZdcSubEpShift();
      mZdcEpManager->initZdcFullEpShift();
      mEpdEpManager->readEpdPhiWgt(); // EPD
      mEpdEpManager->readEpdSideReCtr();
      mEpdEpManager->readEpdSideShift();
      mEpdEpManager->readEpdSideShiftFull();
      mEpdEpManager->initEpdSideResolution();
      mEpdEpManager->initEpdSubEpSideShift();
      mEpdEpManager->initEpdFullEpSideShift();
    }
    if(mAnaCut->isFxt3p85GeV_2018()) 
    {
      // mEpdEpManager->readEpdPhiWgt(); // EPD
      mEpdEpManager->readEpdGrpReCtr();
      mEpdEpManager->readEpdGrpShift();
      mEpdEpManager->initEpdGrpResolution();
      mEpdEpManager->initEpdSubEpGrpShift();
      mMixEpManager->initMixEpRes(); // Mix
      mMixEpManager->initMixSubEpShift();
    }
    mTpcEpManager->readTpcReCtr(); // TPC
    mTpcEpManager->readTpcShift();
    mTpcEpManager->initTpcResolution();
    mTpcEpManager->initTpcSubEpShift();
  }
  if(mMode == 5)
  { // calculate charged hadron v1 from ZDC & EPD and charged hadron v2 and v3 from TPC
    file_mOutPutFlow = new TFile(str_mOutPutFlow.c_str(),"RECREATE");
    if(mAnaCut->isIsobar()) 
    {
      mZdcEpManager->readZdcGain(); // ZDC
      mZdcEpManager->readZdcReCtr();
      mZdcEpManager->readZdcShift();
      mZdcEpManager->readZdcShiftFull();
      mZdcEpManager->readZdcResolution();
      mZdcEpManager->initZdcFullEpFlow();
      mEpdEpManager->readEpdPhiWgt(); // EPD
      mEpdEpManager->readEpdSideReCtr();
      mEpdEpManager->readEpdSideShift();
      mEpdEpManager->readEpdSideShiftFull();
      mEpdEpManager->readEpdSideResolution();
      mEpdEpManager->initEpdSubEpSideFlow();
    }
    if(mAnaCut->isFxt3p85GeV_2018()) 
    {
      // mEpdEpManager->readEpdPhiWgt(); // EPD
      mEpdEpManager->readEpdGrpReCtr();
      mEpdEpManager->readEpdGrpShift();
      // mEpdEpManager->readEpdGrpShiftFull();
      // mEpdEpManager->readEpdGrpResolution();
      mMixEpManager->readMixEpRes(); // Mix
      mMixEpManager->readDeuEfficiency(); // Mix
      mMixEpManager->initMixSubEpFlow();
    }
    mTpcEpManager->readTpcReCtr(); // TPC
    mTpcEpManager->readTpcShift();
    mTpcEpManager->readTpcResolution();
    mTpcEpManager->initTpcSubEpFlow();
  }

  return kStOK;
}

int StEventPlaneMaker::Finish() 
{
  if(mMode == 0)
  {
    if(str_mOutPutGainCorr != "")
    {
      file_mOutPutGainCorr->cd();
      if(mAnaCut->isIsobar()) 
      {
	mZdcEpManager->writeZdcGain(); // ZDC
	mEpdEpManager->writeEpdPhiWgt(); // EPD
	mEpdEpManager->writeEpdSubEpSideRaw();
      }
      if(mAnaCut->isFxt3p85GeV_2018()) 
      {
	mEpdEpManager->writeEpdPhiWgt(); // EPD
	mEpdEpManager->writeEpdSubEpGrpRaw();
      }
      file_mOutPutGainCorr->Close();
    }
  }
  if(mMode == 1)
  {
    if(str_mOutPutReCenterPar != "")
    {
      file_mOutPutReCenterPar->cd();
      if(mAnaCut->isIsobar()) 
      {
	mZdcEpManager->writeZdcReCtr(); // ZDC
	mZdcEpManager->writeZdcSubEpRaw();
	mEpdEpManager->writeEpdSideReCtr(); // EPD
	mEpdEpManager->writeEpdSubEpSideWgt();
      }
      if(mAnaCut->isFxt3p85GeV_2018()) 
      {
	mEpdEpManager->writeEpdGrpReCtr(); // EPD
	mEpdEpManager->writeEpdSubEpGrpWgt();
	mMixEpManager->writeMixSubEpRaw(); // Mix
      }
      mTpcEpManager->writeTpcReCtr(); // TPC
      mTpcEpManager->writeTpcSubEpRaw();
      file_mOutPutReCenterPar->Close();
    }
  }
  if(mMode == 2)
  {
    if(str_mOutPutShiftPar != "")
    {
      file_mOutPutShiftPar->cd();
      if(mAnaCut->isIsobar()) 
      {
	mZdcEpManager->writeZdcShift(); // ZDC
	mZdcEpManager->writeZdcSubEpReCtr();
	mEpdEpManager->writeEpdSideShift(); // EPD
	mEpdEpManager->writeEpdSubEpSideReCtr();
      }
      if(mAnaCut->isFxt3p85GeV_2018()) 
      {
	mEpdEpManager->writeEpdGrpShift(); // EPD
	mEpdEpManager->writeEpdSubEpGrpReCtr();
	mMixEpManager->writeMixSubEpReCtr(); // Mix
      }
      mTpcEpManager->writeTpcShift(); // TPC
      mTpcEpManager->writeTpcSubEpReCtr();
      file_mOutPutShiftPar->Close();
    }
  }
  if(mMode == 3)
  {
    if(str_mOutPutShiftPar != "")
    {
      file_mOutPutShiftPar->cd();
      if(mAnaCut->isIsobar()) 
      {
	mZdcEpManager->writeZdcShiftFull(); // ZDC
	mZdcEpManager->writeZdcSubEpShift();
	mEpdEpManager->writeEpdSideShiftFull(); // EPD
	mEpdEpManager->writeEpdSubEpSideShift();
      }
      file_mOutPutShiftPar->Close();
    }
  }
  if(mMode == 4)
  {
    if(str_mOutPutResolution != "")
    {
      file_mOutPutResolution->cd();
      if(mAnaCut->isIsobar()) 
      {
	mZdcEpManager->writeZdcResolution(); // ZDC
	mZdcEpManager->writeZdcSubEpShift();
	mZdcEpManager->writeZdcFullEpShift();
	mEpdEpManager->writeEpdSideResolution(); // EPD
	mEpdEpManager->writeEpdSubEpSideShift();
	mEpdEpManager->writeEpdFullEpSideShift();
      }
      if(mAnaCut->isFxt3p85GeV_2018()) 
      {
	mEpdEpManager->writeEpdGrpResolution(); // EPD
	mEpdEpManager->writeEpdSubEpGrpShift();
	mMixEpManager->writeMixEpRes(); // Mix
	mMixEpManager->writeMixSubEpShift();
      }
      mTpcEpManager->writeTpcResolution(); // TPC
      mTpcEpManager->writeTpcSubEpShift();
      file_mOutPutResolution->Close();
    }
  }
  if(mMode == 5)
  {
    if(str_mOutPutFlow != "")
    {
      file_mOutPutFlow->cd();
      if(mAnaCut->isIsobar()) 
      {
	mZdcEpManager->writeZdcFullEpFlow(); // ZDC
	mEpdEpManager->writeEpdSubEpSideFlow(); // EPD
      }
      if(mAnaCut->isFxt3p85GeV_2018()) 
      {
	mMixEpManager->writeMixSubEpFlow(); // Mix
      }
      mTpcEpManager->writeTpcSubEpFlow(); // TPC
      file_mOutPutFlow->Close();
    }
  }

  return kStOK;
}

void StEventPlaneMaker::Clear(Option_t *opt) {
}
//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Make() 
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
    const int evtId        = mPicoEvent->eventId();
    int refMult            = mPicoEvent->refMult(); // initialize refMult
    const TVector3 primVtx = mPicoEvent->primaryVertex();
    const double primVz    = mPicoEvent->primaryVertex().z();
    const double zdcX      = mPicoEvent->ZDCx();
    const unsigned short nBTofMatch = mPicoEvent->nBTOFMatch();     // get number of tof match points
    const unsigned int nTracks      = mPicoDst->numberOfTracks();   // get number of tracks
    const unsigned int nEpdHits     = mPicoDst->numberOfEpdHits();  // get number of EPD hits

    const int vzBin  = mAnaUtils->getVzBin(primVz); // 0 for -vz || 1 for +vz
    const int runIdx = mAnaUtils->findRunIndex(runId); // find run index for a specific run
    // cout << "runId = " << runId << ", runIdx = " << runIdx << endl;
    if(runIdx < 0)
    {
      LOG_ERROR << "Could not find this run Index from StAnalysisUtils! Skip!" << endm;
      return kStErr;
    }

    if(!mRefMultCorr) // StRefMultCorr Cut & centrality
    {
      LOG_WARN << " No mRefMultCorr! Skip! " << endl;
      return kStErr;
    }
    mRefMultCorr->init(runId); // get centrality9 and reweight
    mRefMultCorr->initEvent(refMult,primVz,zdcX);
    int cent9     = mRefMultCorr->getCentralityBin9(); // get Centrality9
    double refWgt = mRefMultCorr->getWeight(); // get weight

    if(mAnaCut->isFxt3p85GeV_2018()) // only for Fxt3p85GeV_2018
    { 
      mPileupUtilFxt3p85->initEvent(mPicoDst);
      refMult = mPileupUtilFxt3p85->get_refMultPrim();
      cent9   = mPileupUtilFxt3p85->get_centrality9();
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

    bool isPileUpStAnalysisCut = mAnaCut->isPileUpEvent((double)refMult, (double)nBTofMatch,primVz); // alway return false for Isobar & Fxt3p85GeV_2018
    bool isPileUpStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut((double)refMult, (double)nBTofMatch,primVz); // valid for Isobar
    if(mAnaCut->isFxt3p85GeV_2018()) isPileUpStRefMultCorr = mPileupUtilFxt3p85->isPileupEPD(); // valid only for Fxt3p85GeV_2018
    bool isPileUpEvt = isPileUpStAnalysisCut || isPileUpStRefMultCorr;
    // cout << "isPileUpEvt = " << isPileUpEvt << ", isPileUpStAnalysisCut = " << isPileUpStAnalysisCut << ", isPileUpStRefMultCorr = " << isPileUpStRefMultCorr << endl;

    // cout << "runId = " << runId << ", runIdx = " << runIdx << ", evtId = " << evtId << ", cent9 = " << cent9 << ", isPileUpEvt = " << isPileUpEvt << ", passEventCut = " << mAnaCut->passEventCut(mPicoEvent) << ", isGoodCent9 = " << mAnaCut->isGoodCent9(cent9) << endl;

    if( !isPileUpEvt && mAnaCut->isGoodCent9(cent9) && mAnaCut->passEventCut(mPicoEvent) )
    { // apply Event Cuts for anlaysis 
      mZdcEpManager->initZdcEpManager(cent9,runIdx,vzBin); // initialize ZDC EP Manager
      mEpdEpManager->initEpdEpManager(cent9,runIdx,vzBin); // initialize EPD EP Manager
      mTpcEpManager->initTpcEpManager(cent9,runIdx,vzBin); // initialize TPC EP Manager
      mMixEpManager->initMixEpManager(cent9,runIdx,vzBin); // initialize Mix EP Manager

      // Calculate QVector
      if( mAnaCut->isIsobar() )
      {
	for(int iSlat = 0; iSlat < 8; ++iSlat) // ZDC QVector
	{ 
	  if(mMode == 0)
	  { // fill Gain Correction Factors for ZDC
	    mZdcEpManager->setZdcSmd(0,0,iSlat,mPicoEvent->ZdcSmdEastVertical(iSlat)); // read in raw ADC value from ZDC-SMD
	    mZdcEpManager->setZdcSmd(0,1,iSlat,mPicoEvent->ZdcSmdEastHorizontal(iSlat));
	    mZdcEpManager->setZdcSmd(1,0,iSlat,mPicoEvent->ZdcSmdWestVertical(iSlat));
	    mZdcEpManager->setZdcSmd(1,1,iSlat,mPicoEvent->ZdcSmdWestHorizontal(iSlat));
	  }
	  if(mMode > 0)
	  { // read and apply gain correction to raw ADC value for ZDC
	    mZdcEpManager->setZdcSmdGainCorr(0,0,iSlat,mPicoEvent->ZdcSmdEastVertical(iSlat));
	    mZdcEpManager->setZdcSmdGainCorr(0,1,iSlat,mPicoEvent->ZdcSmdEastHorizontal(iSlat));
	    mZdcEpManager->setZdcSmdGainCorr(1,0,iSlat,mPicoEvent->ZdcSmdWestVertical(iSlat));
	    mZdcEpManager->setZdcSmdGainCorr(1,1,iSlat,mPicoEvent->ZdcSmdWestHorizontal(iSlat));
	  }
	}

	for(unsigned int iEpdHit = 0; iEpdHit < nEpdHits; ++iEpdHit) // EPD QVector
	{
	  StPicoEpdHit *picoEpdHit = (StPicoEpdHit*)mPicoDst->epdHit(iEpdHit); 
	  if(!picoEpdHit) continue;

	  if(mMode == 0)
	  { // fill phi weight for EPD
	    if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	    {
	      mEpdEpManager->addHitSideRawEast(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdPhiWgtEast(picoEpdHit);
	    }
	    if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	    {
	      mEpdEpManager->addHitSideRawWest(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdPhiWgtWest(picoEpdHit);
	    }
	  }
	  if(mMode == 1) // calculate phi and eta weighted (if avaliable) Q1Vector from EPD
	  { // fill recenter correction for EPD
	    if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	    {
	      mEpdEpManager->addHitSideWgtEast(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdSideReCtrEast(picoEpdHit, primVtx); // fill EPD ReCenter Parameters East
	    }
	    if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	    {
	      mEpdEpManager->addHitSideWgtWest(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdSideReCtrWest(picoEpdHit, primVtx); // fill EPD ReCenter Parameters West
	    }
	  }
	  if(mMode > 1) // calculate phi & eta weighted (if avaliable) and recentered Q1Vector from EPD
	  {
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
      }

      for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack) // calculate Q1Vector & Q2Vector & Q3Vector from TPC
      { // track loop
	StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	if(!picoTrack) continue;

	if(mMode == 1) // calculate raw Q1Vector, Q2Vector and Q3Vector from TPC
	{
	  if( mAnaCut->passTrkTpcEpEast(picoTrack, primVtx) ) // negative eta
	  {
	    mTpcEpManager->addTrackRawEast(picoTrack);
	    mTpcEpManager->fillTpcReCtrEast(picoTrack); // fill TPC ReCenter Parameters East
	  }
	  if( mAnaCut->passTrkTpcEpWest(picoTrack, primVtx) ) // positive eta
	  {
	    mTpcEpManager->addTrackRawWest(picoTrack);
	    mTpcEpManager->fillTpcReCtrWest(picoTrack); // fill TPC ReCenter Parameters West
	  }
	  if( mAnaCut->passTrkTpcEpFull(picoTrack, primVtx) )
	  {
	    mTpcEpManager->addTrackRawFull(picoTrack);
	    mTpcEpManager->fillTpcReCtrFull(picoTrack); // fill TPC ReCenter Parameters Full
	  }
	}
	if(mMode > 1 && mMode != 3) // calculate recentered Q1Vector, Q2Vector and Q3Vector from TPC
	{
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
      }

      if( mAnaCut->isFxt3p85GeV_2018() )
      {
	for(unsigned int iEpdHit = 0; iEpdHit < nEpdHits; ++iEpdHit) // EPD QVector
	{
	  StPicoEpdHit *picoEpdHit = (StPicoEpdHit*)mPicoDst->epdHit(iEpdHit); 
	  if(!picoEpdHit) continue;

	  if(mMode == 0)
	  { // fill phi weight for EPD
	    if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	    {
	      mEpdEpManager->addHitGrpRawEast(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdPhiWgtEast(picoEpdHit);
	    }
	    if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	    {
	      mEpdEpManager->addHitGrpRawWest(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdPhiWgtWest(picoEpdHit);
	    }
	  }
	  if(mMode == 1) // calculate phi and eta weighted (if avaliable) Q1Vector from EPD
	  { // fill recentered correction for EPD
	    if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	    {
	      mEpdEpManager->addHitGrpRawEast(picoEpdHit, primVtx); // used by Evt Ave ReCenter Mehtod
	      mEpdEpManager->addHitGrpWgtEast(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdGrpReCtrTrkAveEast(picoEpdHit, primVtx);
	    }
	    if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	    {
	      mEpdEpManager->addHitGrpRawWest(picoEpdHit, primVtx); // used by Evt Ave ReCenter Mehtod
	      mEpdEpManager->addHitGrpWgtWest(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdGrpReCtrTrkAveWest(picoEpdHit, primVtx);
	    }
	  }
	  if(mMode > 1) // calculate phi & eta weighted (if avaliable) and recentered Q1Vector from EPD
	  {
	    if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	    {
	      mEpdEpManager->addHitGrpRawEast(picoEpdHit, primVtx); // used by Evt Ave ReCenter Mehtod
	      mEpdEpManager->addHitGrpReCtrTrkAveEast(picoEpdHit, primVtx);
	    }
	    if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	    {
	      mEpdEpManager->addHitGrpRawWest(picoEpdHit, primVtx); // used by Evt Ave ReCenter Mehtod
	      mEpdEpManager->addHitGrpReCtrTrkAveWest(picoEpdHit, primVtx);
	    }
	  }
	}
      }

      if(mMode == 0) // fill gain correction parameters for ZDC & EPD
      {
	if( mAnaCut->isIsobar() )
	{
	  for(int iEastWest = 0; iEastWest < 2; ++iEastWest) // fill ZDC Gain Correction Histograms
	  {
	    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
	    {
	      for(int iSlat = 0; iSlat < 8; ++iSlat)
	      {
		mZdcEpManager->fillZdcGain(iEastWest,iVertHori,iSlat,mZdcEpManager->getZdcSmd(iEastWest,iVertHori,iSlat));
		// std::cout << "iEastWest = " << iEastWest << ", iVertHori = " << iVertHori << ", iSlat = " << iSlat << ", zdc = " << mZdcEpManager->getZdcSmd(iEastWest,iVertHori,iSlat) << std::endl;
	      }
	    }
	  }

	  TVector2 vQ1EpdSideEast = mEpdEpManager->getQ1VecSideRawEast(); // get Q1Vector from EPD
	  TVector2 vQ1EpdSideWest = mEpdEpManager->getQ1VecSideRawWest();
	  TVector2 vQ1EpdSideFull = mEpdEpManager->getQ1VecSideRawFull();
	  if( mAnaCut->passQVecEpdSide(vQ1EpdSideEast,vQ1EpdSideWest,vQ1EpdSideFull) ) // EPD EP
	  {
	    const double Psi1SideRawEast = mEpdEpManager->getPsi1SideRawEast();
	    const double Psi1SideRawWest = mEpdEpManager->getPsi1SideRawWest();
	    const double Psi1SideRawFull = mEpdEpManager->getPsi1SideRawFull();
	    mEpdEpManager->fillEpdSubEpSideRaw(Psi1SideRawEast,Psi1SideRawWest,Psi1SideRawFull);
	  }
	}

	if( mAnaCut->isFxt3p85GeV_2018() )
	{
	  TVector2 vQ1EpdGrpEast[mNumRingsGrps], vQ1EpdGrpWest[mNumRingsGrps], vQ1EpdGrpFull[mNumRingsGrps];
	  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	  {
	    vQ1EpdGrpEast[iGrp] = mEpdEpManager->getQ1VecGrpRawEast(iGrp); // get Q1Vector from EPD groups
	    vQ1EpdGrpWest[iGrp] = mEpdEpManager->getQ1VecGrpRawWest(iGrp);
	    vQ1EpdGrpFull[iGrp] = mEpdEpManager->getQ1VecGrpRawFull(iGrp);
	    if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[iGrp],vQ1EpdGrpWest[iGrp],vQ1EpdGrpFull[iGrp],iGrp) ) // EPD EP
	    {
	      const double Psi1GrpRawEast = mEpdEpManager->getPsi1GrpRawEast(iGrp);
	      const double Psi1GrpRawWest = mEpdEpManager->getPsi1GrpRawWest(iGrp);
	      const double Psi1GrpRawFull = mEpdEpManager->getPsi1GrpRawFull(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpRaw(Psi1GrpRawEast,Psi1GrpRawWest,Psi1GrpRawFull,iGrp);

	      // cout << "iGrp = " << iGrp << ", Psi1GrpRawEast = " << Psi1GrpRawEast << ", Psi1GrpRawWest = " << Psi1GrpRawWest << ", Psi1GrpRawFull = " << Psi1GrpRawFull << endl;
	    }
	  }
	}
      }
      if(mMode == 1) // fill recenter correction parameter for ZDC & EPD & TPC Sub EP
      {
	if( mAnaCut->isIsobar() )
	{
	  TVector2 vQ1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 vQ1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 vQ1ZdcFull = mZdcEpManager->getQ1VecFull(vQ1ZdcEast,vQ1ZdcWest,mMode); // TVector2 vQ1ZdcFull = vQ1ZdcWest-vQ1ZdcEast;
	  if( mAnaCut->passQVecZdc(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull) ) // ZDC EP
	  {
	    mZdcEpManager->fillZdcSubEpRaw(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull);
	    mZdcEpManager->fillZdcReCtrEast(vQ1ZdcEast);
	    mZdcEpManager->fillZdcReCtrWest(vQ1ZdcWest);
	  }

	  TVector2 vQ1EpdSideEast = mEpdEpManager->getQ1VecSideWgtEast(); // get Q1Vector from EPD
	  TVector2 vQ1EpdSideWest = mEpdEpManager->getQ1VecSideWgtWest();
	  TVector2 vQ1EpdSideFull = mEpdEpManager->getQ1VecSideWgtFull();
	  if( mAnaCut->passQVecEpdSide(vQ1EpdSideEast,vQ1EpdSideWest,vQ1EpdSideFull) ) // EPD EP
	  {
	    const double Psi1SideWgtEast = mEpdEpManager->getPsi1SideWgtEast();
	    const double Psi1SideWgtWest = mEpdEpManager->getPsi1SideWgtWest();
	    const double Psi1SideWgtFull = mEpdEpManager->getPsi1SideWgtFull();
	    mEpdEpManager->fillEpdSubEpSideWgt(Psi1SideWgtEast, Psi1SideWgtWest, Psi1SideWgtFull);
	  }
	}

	const int numTrkRawEast = mTpcEpManager->getNumTrkRawEast(); // TPC EP
	const int numTrkRawWest = mTpcEpManager->getNumTrkRawWest();
	if(mAnaCut->passNumTrkTpcSubEpRaw(numTrkRawEast, numTrkRawWest))
	{
	  const double Psi1RawEast = mTpcEpManager->getPsi1RawEast(); // 1st EP
	  const double Psi1RawWest = mTpcEpManager->getPsi1RawWest();
	  const double Psi1RawFull = mTpcEpManager->getPsi1RawFull();
	  const double Psi2RawEast = mTpcEpManager->getPsi2RawEast(); // 2nd EP
	  const double Psi2RawWest = mTpcEpManager->getPsi2RawWest();
	  const double Psi2RawFull = mTpcEpManager->getPsi2RawFull();
	  const double Psi3RawEast = mTpcEpManager->getPsi3RawEast(); // 3rd EP
	  const double Psi3RawWest = mTpcEpManager->getPsi3RawWest();
	  const double Psi3RawFull = mTpcEpManager->getPsi3RawFull();
	  mTpcEpManager->fillTpcSubEpRaw(Psi1RawEast, Psi1RawWest, Psi1RawFull, Psi2RawEast, Psi2RawWest, Psi2RawFull, Psi3RawEast, Psi3RawWest, Psi3RawFull);
	}

	if( mAnaCut->isFxt3p85GeV_2018() )
	{
	  TVector2 vQ1EpdGrpEast[mNumRingsGrps], vQ1EpdGrpWest[mNumRingsGrps], vQ1EpdGrpFull[mNumRingsGrps];
	  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	  {
	    vQ1EpdGrpEast[iGrp] = mEpdEpManager->getQ1VecGrpWgtEast(iGrp); // get Q1Vector from EPD groups
	    vQ1EpdGrpWest[iGrp] = mEpdEpManager->getQ1VecGrpWgtWest(iGrp);
	    vQ1EpdGrpFull[iGrp] = mEpdEpManager->getQ1VecGrpWgtFull(iGrp);
	    if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[iGrp],vQ1EpdGrpWest[iGrp],vQ1EpdGrpFull[iGrp],iGrp) ) // EPD EP
	    {
	      const double Psi1GrpWgtEast = mEpdEpManager->getPsi1GrpWgtEast(iGrp);
	      const double Psi1GrpWgtWest = mEpdEpManager->getPsi1GrpWgtWest(iGrp);
	      const double Psi1GrpWgtFull = mEpdEpManager->getPsi1GrpWgtFull(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpWgt(Psi1GrpWgtEast, Psi1GrpWgtWest, Psi1GrpWgtFull, iGrp);
	      mEpdEpManager->fillEpdGrpReCtrEvtAveEast(iGrp); // fill Evt Ave ReCenter Parameter
	      mEpdEpManager->fillEpdGrpReCtrEvtAveWest(iGrp);
	      // TVector2 Q1VecEastDisplay = -1.0*vQ1EpdGrpEast[iGrp];
	      // cout << "runId = " << runId << ", runIdx = " << runIdx << ", evtId = " << evtId << ", cent9 = " << cent9 << ", iGrp = " << iGrp << ", Psi1GrpWgtEast = " << Q1VecEastDisplay.Phi() << endl;
	    }
	  }
	  // cout << endl;

	  if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[0],vQ1EpdGrpWest[0],vQ1EpdGrpFull[0],0) &&
	      mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[1],vQ1EpdGrpWest[1],vQ1EpdGrpFull[1],1) &&
	      mAnaCut->passNumTrkTpcSubEpRaw(numTrkRawEast, numTrkRawWest) )
	  { // 3-sub events method
	    const double Psi1EpdGrp0 = mEpdEpManager->getPsi1GrpWgtEast(0); // Psi1 from EPD Grp 0 East
	    const double Psi1EpdGrp1 = mEpdEpManager->getPsi1GrpWgtEast(1); // Psi1 from EPD Grp 1 East
	    const double Psi1TpcEast = mTpcEpManager->getPsi1RawEast(); // Psi1 from TPC East
	    const double Psi1TpcWest = mTpcEpManager->getPsi1RawWest(); // Psi1 from TPC West

	    mMixEpManager->fillMixSubEpRaw(Psi1EpdGrp0,Psi1EpdGrp1,Psi1TpcEast,Psi1TpcWest);
	  }
	}
      }
      if(mMode == 2) // fill shift correction parameter for ZDC & EPD & TPC Sub EP
      {
	if( mAnaCut->isIsobar() )
	{
	  TVector2 vQ1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 vQ1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 vQ1ZdcFull = mZdcEpManager->getQ1VecFull(vQ1ZdcEast,vQ1ZdcWest,mMode); // TVector2 vQ1ZdcFull = vQ1ZdcWest-vQ1ZdcEast;
	  if( mAnaCut->passQVecZdc(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull) ) // ZDC EP
	  {
	    mZdcEpManager->fillZdcSubEpReCtr(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull);
	    mZdcEpManager->fillZdcShiftEast(vQ1ZdcEast);
	    mZdcEpManager->fillZdcShiftWest(vQ1ZdcWest);
	  }

	  TVector2 vQ1EpdSideEast = mEpdEpManager->getQ1VecSideReCtrEast(); // get Q1Vector from EPD
	  TVector2 vQ1EpdSideWest = mEpdEpManager->getQ1VecSideReCtrWest();
	  TVector2 vQ1EpdSideFull = mEpdEpManager->getQ1VecSideReCtrFull();
	  if( mAnaCut->passQVecEpdSide(vQ1EpdSideEast,vQ1EpdSideWest,vQ1EpdSideFull) ) // EPD EP
	  {
	    const double Psi1SideReCtrEast = mEpdEpManager->getPsi1SideReCtrEast(); // EPD EP
	    const double Psi1SideReCtrWest = mEpdEpManager->getPsi1SideReCtrWest();
	    const double Psi1SideReCtrFull = mEpdEpManager->getPsi1SideReCtrFull();
	    mEpdEpManager->fillEpdSubEpSideReCtr(Psi1SideReCtrEast, Psi1SideReCtrWest, Psi1SideReCtrFull);
	    mEpdEpManager->fillEpdSideShiftEast();
	    mEpdEpManager->fillEpdSideShiftWest();
	  }
	}

	const int numTrkReCtrEast = mTpcEpManager->getNumTrkReCtrEast(); // TPC EP
	const int numTrkReCtrWest = mTpcEpManager->getNumTrkReCtrWest();
	if(mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest))
	{
	  const double Psi1ReCtrEast = mTpcEpManager->getPsi1ReCtrEast(); // 1st EP
	  const double Psi1ReCtrWest = mTpcEpManager->getPsi1ReCtrWest();
	  const double Psi1ReCtrFull = mTpcEpManager->getPsi1ReCtrFull();
	  const double Psi2ReCtrEast = mTpcEpManager->getPsi2ReCtrEast(); // 2nd EP
	  const double Psi2ReCtrWest = mTpcEpManager->getPsi2ReCtrWest();
	  const double Psi2ReCtrFull = mTpcEpManager->getPsi2ReCtrFull();
	  const double Psi3ReCtrEast = mTpcEpManager->getPsi3ReCtrEast(); // 3rd EP
	  const double Psi3ReCtrWest = mTpcEpManager->getPsi3ReCtrWest();
	  const double Psi3ReCtrFull = mTpcEpManager->getPsi3ReCtrFull();
	  mTpcEpManager->fillTpcSubEpReCtr(Psi1ReCtrEast, Psi1ReCtrWest, Psi1ReCtrFull, Psi2ReCtrEast, Psi2ReCtrWest, Psi2ReCtrFull, Psi3ReCtrEast, Psi3ReCtrWest, Psi3ReCtrFull);
	  mTpcEpManager->fillTpcShiftEast();
	  mTpcEpManager->fillTpcShiftWest();
	  mTpcEpManager->fillTpcShiftFull();
	}

	if( mAnaCut->isFxt3p85GeV_2018() )
	{
	  TVector2 vQ1EpdGrpTrkAveEast[mNumRingsGrps], vQ1EpdGrpTrkAveWest[mNumRingsGrps], vQ1EpdGrpTrkAveFull[mNumRingsGrps];
	  TVector2 vQ1EpdGrpEvtAveEast[mNumRingsGrps], vQ1EpdGrpEvtAveWest[mNumRingsGrps], vQ1EpdGrpEvtAveFull[mNumRingsGrps];
	  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	  {
	    vQ1EpdGrpTrkAveEast[iGrp] = mEpdEpManager->getQ1VecGrpReCtrTrkAveEast(iGrp); // get Trk Ave ReCtr Q1Vector from EPD groups
	    vQ1EpdGrpTrkAveWest[iGrp] = mEpdEpManager->getQ1VecGrpReCtrTrkAveWest(iGrp);
	    vQ1EpdGrpTrkAveFull[iGrp] = mEpdEpManager->getQ1VecGrpReCtrTrkAveFull(iGrp);
	    vQ1EpdGrpEvtAveEast[iGrp] = mEpdEpManager->getQ1VecGrpReCtrEvtAveEast(iGrp); // get Evt Ave ReCtr Q1Vector from EPD groups
	    vQ1EpdGrpEvtAveWest[iGrp] = mEpdEpManager->getQ1VecGrpReCtrEvtAveWest(iGrp);
	    vQ1EpdGrpEvtAveFull[iGrp] = mEpdEpManager->getQ1VecGrpReCtrEvtAveFull(iGrp);
	    if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpTrkAveEast[iGrp],vQ1EpdGrpTrkAveWest[iGrp],vQ1EpdGrpTrkAveFull[iGrp],iGrp) ) // EPD EP
	    {
	      const double Psi1GrpReCtrTrkAveEast = mEpdEpManager->getPsi1GrpReCtrTrkAveEast(iGrp);
	      const double Psi1GrpReCtrTrkAveWest = mEpdEpManager->getPsi1GrpReCtrTrkAveWest(iGrp);
	      const double Psi1GrpReCtrTrkAveFull = mEpdEpManager->getPsi1GrpReCtrTrkAveFull(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpReCtrTrkAve(Psi1GrpReCtrTrkAveEast, Psi1GrpReCtrTrkAveWest, Psi1GrpReCtrTrkAveFull, iGrp);
	      mEpdEpManager->fillEpdGrpShiftTrkAveEast(iGrp); // fill Trk Ave Shift Parameter
	      mEpdEpManager->fillEpdGrpShiftTrkAveWest(iGrp);

	      const double Psi1GrpReCtrEvtAveEast = mEpdEpManager->getPsi1GrpReCtrEvtAveEast(iGrp);
	      const double Psi1GrpReCtrEvtAveWest = mEpdEpManager->getPsi1GrpReCtrEvtAveWest(iGrp);
	      const double Psi1GrpReCtrEvtAveFull = mEpdEpManager->getPsi1GrpReCtrEvtAveFull(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpReCtrEvtAve(Psi1GrpReCtrEvtAveEast, Psi1GrpReCtrEvtAveWest, Psi1GrpReCtrEvtAveFull, iGrp);
	      mEpdEpManager->fillEpdGrpShiftEvtAveEast(iGrp); // fill Evt Ave Shift Parameter
	      mEpdEpManager->fillEpdGrpShiftEvtAveWest(iGrp);

	      // cout << "runId = " << runId << ", evtId = " << evtId << ", iGrp = " << iGrp << ", Psi1GrpReCtrEvtAveEast = " << Psi1GrpReCtrEvtAveEast << endl;
	      // TVector2 Q1VecEastDisplay = -1.0*vQ1EpdGrpEvtAveEast[iGrp];
	      // cout << "runId = " << runId << ", runIdx = " << runIdx << ", evtId = " << evtId << ", cent9 = " << cent9 << ", iGrp = " << iGrp << ", Psi1GrpReCtrEvtAveEast = " << Q1VecEastDisplay.Phi() << endl;
	    }
	  }
	  // cout << endl;

	  if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpTrkAveEast[0],vQ1EpdGrpTrkAveWest[0],vQ1EpdGrpTrkAveFull[0],0) &&
	      mAnaCut->passQVecEpdGrp(vQ1EpdGrpTrkAveEast[1],vQ1EpdGrpTrkAveWest[1],vQ1EpdGrpTrkAveFull[1],1) &&
	      mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest) )
	  { // 3-sub events method
	    const double Psi1EpdGrp0 = mEpdEpManager->getPsi1GrpReCtrTrkAveEast(0); // Psi1 from EPD Grp 0 East
	    const double Psi1EpdGrp1 = mEpdEpManager->getPsi1GrpReCtrTrkAveEast(1); // Psi1 from EPD Grp 1 East
	    const double Psi1TpcEast = mTpcEpManager->getPsi1ReCtrEast(); // Psi1 from TPC East
	    const double Psi1TpcWest = mTpcEpManager->getPsi1ReCtrWest(); // Psi1 from TPC West

	    mMixEpManager->fillMixSubEpReCtr(Psi1EpdGrp0,Psi1EpdGrp1,Psi1TpcEast,Psi1TpcWest);
	  }
	}
      }
      if(mMode == 3) // fill shift correction parameter for ZDC & EPD Full EP
      {
	if( mAnaCut->isIsobar() )
	{
	  TVector2 vQ1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 vQ1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 vQ1ZdcFull = mZdcEpManager->getQ1VecFull(vQ1ZdcEast,vQ1ZdcWest,mMode); // TVector2 vQ1ZdcFull = vQ1ZdcWest-vQ1ZdcEast;
	  if( mAnaCut->passQVecZdc(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull) ) // ZDC EP
	  {
	    mZdcEpManager->fillZdcSubEpShift(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull);
	    mZdcEpManager->fillZdcShiftFull(vQ1ZdcFull);
	  }

	  TVector2 vQ1EpdSideEast = mEpdEpManager->getQ1VecSideShiftEast(); // get Q1Vector from EPD
	  TVector2 vQ1EpdSideWest = mEpdEpManager->getQ1VecSideShiftWest();
	  TVector2 vQ1EpdSideFull = mEpdEpManager->getQ1VecSideShiftFull();
	  if( mAnaCut->passQVecEpdSide(vQ1EpdSideEast,vQ1EpdSideWest,vQ1EpdSideFull) ) // EPD EP
	  {
	    const double Psi1SideShiftEast = mEpdEpManager->getPsi1SideShiftEast();
	    const double Psi1SideShiftWest = mEpdEpManager->getPsi1SideShiftWest();
	    const double Psi1SideShiftFull = mEpdEpManager->getPsi1SideShiftFull();
	    mEpdEpManager->fillEpdSubEpSideShift(Psi1SideShiftEast, Psi1SideShiftWest, Psi1SideShiftFull);
	    mEpdEpManager->fillEpdSideShiftFull();
	  }
	}
      }
      if(mMode == 4) // fill event plane resolution for ZDC-SMD & EPD & TPC Sub EP
      {
	if( mAnaCut->isIsobar() )
	{
	  TVector2 vQ1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 vQ1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 vQ1ZdcFull = mZdcEpManager->getQ1VecFull(vQ1ZdcEast,vQ1ZdcWest,mMode); // TVector2 vQ1ZdcFull = vQ1ZdcWest-vQ1ZdcEast;
	  if( mAnaCut->passQVecZdc(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull) ) // ZDC EP
	  {
	    mZdcEpManager->fillZdcSubEpShift(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull);
	    mZdcEpManager->fillZdcFullEpShift(vQ1ZdcFull);
	    mZdcEpManager->fillZdcResolution(vQ1ZdcEast,vQ1ZdcWest);
	  }

	  TVector2 vQ1EpdSideEast = mEpdEpManager->getQ1VecSideShiftEast(); // get Q1Vector from EPD
	  TVector2 vQ1EpdSideWest = mEpdEpManager->getQ1VecSideShiftWest();
	  TVector2 vQ1EpdSideFull = mEpdEpManager->getQ1VecSideShiftFull();
	  // TVector2 vQ1EpdSideFullCorr = mEpdEpManager->getQ1VecSideShiftFullCorr();
	  if( mAnaCut->passQVecEpdSide(vQ1EpdSideEast,vQ1EpdSideWest,vQ1EpdSideFull) ) // EPD EP
	  {
	    const double Psi1SideShiftEast     = mEpdEpManager->getPsi1SideShiftEast();
	    const double Psi1SideShiftWest     = mEpdEpManager->getPsi1SideShiftWest();
	    const double Psi1SideShiftFull     = mEpdEpManager->getPsi1SideShiftFull();
	    const double Psi1SideShiftFullCorr = mEpdEpManager->getPsi1SideShiftFullCorr();
	    mEpdEpManager->fillEpdSubEpSideShift(Psi1SideShiftEast, Psi1SideShiftWest, Psi1SideShiftFull);
	    mEpdEpManager->fillEpdFullEpSideShift(Psi1SideShiftFullCorr);
	    mEpdEpManager->fillEpdSideResolution(Psi1SideShiftEast, Psi1SideShiftWest);
	  }
	}

	const int numTrkReCtrEast = mTpcEpManager->getNumTrkReCtrEast(); // TPC EP
	const int numTrkReCtrWest = mTpcEpManager->getNumTrkReCtrWest();
	if(mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest))
	{
	  const TVector2 vQ1TpcEast = mTpcEpManager->getQ1VecReCtrEast(); // ReCentered Q1Vector
	  const TVector2 vQ1TpcWest = mTpcEpManager->getQ1VecReCtrWest();
	  const TVector2 vQ1TpcFull = mTpcEpManager->getQ1VecReCtrFull();
	  const TVector2 vQ2TpcEast = mTpcEpManager->getQ2VecReCtrEast(); // ReCentered Q2Vector
	  const TVector2 vQ2TpcWest = mTpcEpManager->getQ2VecReCtrWest();
	  const TVector2 vQ2TpcFull = mTpcEpManager->getQ2VecReCtrFull();
	  const TVector2 vQ3TpcEast = mTpcEpManager->getQ3VecReCtrEast(); // ReCentered Q3Vector
	  const TVector2 vQ3TpcWest = mTpcEpManager->getQ3VecReCtrWest();
	  const TVector2 vQ3TpcFull = mTpcEpManager->getQ3VecReCtrFull();

	  const double Psi1ShiftEast = mTpcEpManager->getPsi1ShiftEast(vQ1TpcEast); // 1st EP
	  const double Psi1ShiftWest = mTpcEpManager->getPsi1ShiftWest(vQ1TpcWest);
	  const double Psi1ShiftFull = mTpcEpManager->getPsi1ShiftFull(vQ1TpcFull);
	  const double Psi2ShiftEast = mTpcEpManager->getPsi2ShiftEast(vQ2TpcEast); // 2nd EP
	  const double Psi2ShiftWest = mTpcEpManager->getPsi2ShiftWest(vQ2TpcWest);
	  const double Psi2ShiftFull = mTpcEpManager->getPsi2ShiftFull(vQ2TpcFull);
	  const double Psi3ShiftEast = mTpcEpManager->getPsi3ShiftEast(vQ3TpcEast); // 3rd EP
	  const double Psi3ShiftWest = mTpcEpManager->getPsi3ShiftWest(vQ3TpcWest);
	  const double Psi3ShiftFull = mTpcEpManager->getPsi3ShiftFull(vQ3TpcFull);

	  mTpcEpManager->fillTpcSubEpShift(Psi1ShiftEast, Psi1ShiftWest, Psi1ShiftFull, Psi2ShiftEast, Psi2ShiftWest, Psi2ShiftFull, Psi3ShiftEast, Psi3ShiftWest, Psi3ShiftFull);
	  mTpcEpManager->fillTpcResolution(Psi1ShiftEast, Psi1ShiftWest, Psi2ShiftEast, Psi2ShiftWest, Psi3ShiftEast, Psi3ShiftWest);
	}

	if( mAnaCut->isFxt3p85GeV_2018() )
	{
	  TVector2 vQ1EpdGrpTrkAveEast[mNumRingsGrps], vQ1EpdGrpTrkAveWest[mNumRingsGrps], vQ1EpdGrpTrkAveFull[mNumRingsGrps]; 
	  TVector2 vQ1EpdGrpEvtAveEast[mNumRingsGrps], vQ1EpdGrpEvtAveWest[mNumRingsGrps], vQ1EpdGrpEvtAveFull[mNumRingsGrps]; 
	  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	  {
	    vQ1EpdGrpTrkAveEast[iGrp] = mEpdEpManager->getQ1VecGrpShiftTrkAveEast(iGrp); // get Trk Ave ReCtr Q1Vector from EPD
	    vQ1EpdGrpTrkAveWest[iGrp] = mEpdEpManager->getQ1VecGrpShiftTrkAveWest(iGrp);
	    vQ1EpdGrpTrkAveFull[iGrp] = mEpdEpManager->getQ1VecGrpShiftTrkAveFull(iGrp);
	    vQ1EpdGrpEvtAveEast[iGrp] = mEpdEpManager->getQ1VecGrpShiftEvtAveEast(iGrp); // get Evt Ave ReCtr Q1Vector from EPD
	    vQ1EpdGrpEvtAveWest[iGrp] = mEpdEpManager->getQ1VecGrpShiftEvtAveWest(iGrp);
	    vQ1EpdGrpEvtAveFull[iGrp] = mEpdEpManager->getQ1VecGrpShiftEvtAveFull(iGrp);
	    if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpTrkAveEast[iGrp],vQ1EpdGrpTrkAveWest[iGrp],vQ1EpdGrpTrkAveFull[iGrp],iGrp) ) // EPD EP
	    {
	      const double Psi1GrpShiftTrkAveEast = mEpdEpManager->getPsi1GrpShiftTrkAveEast(iGrp);
	      const double Psi1GrpShiftTrkAveWest = mEpdEpManager->getPsi1GrpShiftTrkAveWest(iGrp);
	      const double Psi1GrpShiftTrkAveFull = mEpdEpManager->getPsi1GrpShiftTrkAveFull(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpShiftTrkAve(Psi1GrpShiftTrkAveEast, Psi1GrpShiftTrkAveWest, Psi1GrpShiftTrkAveFull, iGrp);
	      mEpdEpManager->fillEpdGrpResolution(Psi1GrpShiftTrkAveEast, Psi1GrpShiftTrkAveWest, iGrp);

	      const double Psi1GrpShiftEvtAveEast = mEpdEpManager->getPsi1GrpShiftEvtAveEast(iGrp);
	      const double Psi1GrpShiftEvtAveWest = mEpdEpManager->getPsi1GrpShiftEvtAveWest(iGrp);
	      const double Psi1GrpShiftEvtAveFull = mEpdEpManager->getPsi1GrpShiftEvtAveFull(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpShiftEvtAve(Psi1GrpShiftEvtAveEast, Psi1GrpShiftEvtAveWest, Psi1GrpShiftEvtAveFull, iGrp);

	      // cout << "runId = " << runId << ", evtId = " << evtId << ", iGrp = " << iGrp << ", Psi1GrpShiftEvtAveEast = " << Psi1GrpShiftEvtAveEast << endl;
	      // TVector2 Q1VecEastDisplay = -1.0*vQ1EpdGrpEvtAveEast[iGrp];
	      // cout << "runId = " << runId << ", runIdx = " << runIdx << ", evtId = " << evtId << ", cent9 = " << cent9 << ", iGrp = " << iGrp << ", Psi1GrpShiftEvtAveEast = " << Q1VecEastDisplay.Phi() << endl;
	    }
	  }
	  // cout << endl;

	  if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpTrkAveEast[0],vQ1EpdGrpTrkAveWest[0],vQ1EpdGrpTrkAveFull[0],0) &&
	      mAnaCut->passQVecEpdGrp(vQ1EpdGrpTrkAveEast[1],vQ1EpdGrpTrkAveWest[1],vQ1EpdGrpTrkAveFull[1],1) )
	  // if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpTrkAveEast[0],vQ1EpdGrpTrkAveWest[0],vQ1EpdGrpTrkAveFull[0],0) &&
	  //     mAnaCut->passQVecEpdGrp(vQ1EpdGrpTrkAveEast[1],vQ1EpdGrpTrkAveWest[1],vQ1EpdGrpTrkAveFull[1],1) &&
	  //     mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest) )
	  { // 3-sub events method
	    const double Psi1EpdGrp0 = mEpdEpManager->getPsi1GrpShiftTrkAveEast(0); // Psi1 from EPD Grp 0 East
	    const double Psi1EpdGrp1 = mEpdEpManager->getPsi1GrpShiftTrkAveEast(1); // Psi1 from EPD Grp 1 East

	    const TVector2 vQ1TpcEast = mTpcEpManager->getQ1VecReCtrEast();
	    const TVector2 vQ1TpcWest = mTpcEpManager->getQ1VecReCtrWest();
	    const double Psi1TpcEast  = mTpcEpManager->getPsi1ShiftEast(vQ1TpcEast); // Psi1 from TPC East
	    const double Psi1TpcWest  = mTpcEpManager->getPsi1ShiftWest(vQ1TpcWest); // Psi1 from TPC West

	    mMixEpManager->fillMixEpRes(Psi1EpdGrp0,Psi1EpdGrp1,Psi1TpcEast,Psi1TpcWest);
	    mMixEpManager->fillMixSubEpShift(Psi1EpdGrp0,Psi1EpdGrp1,Psi1TpcEast,Psi1TpcWest);
	  }
	}
      }
      if(mMode == 5) // calculate charged hadron v1 from ZDC & EPD and charged hadron v2 and v3 from TPC
      {
	if( mAnaCut->isIsobar() )
	{
	  TVector2 vQ1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 vQ1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 vQ1ZdcFull = mZdcEpManager->getQ1VecFull(vQ1ZdcEast,vQ1ZdcWest,mMode); // TVector2 vQ1ZdcFull = vQ1ZdcWest-vQ1ZdcEast;
	  if( mAnaCut->passQVecZdc(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull)  && mZdcEpManager->getZdcFullEpResVal(cent9) > 0.0 ) // ZDC EP
	  { // charged hadron v1 from ZDC
	    const double Psi1ZdcFull = TMath::ATan2(vQ1ZdcFull.Y(),vQ1ZdcFull.X()); // -pi to pi
	    for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack)	  
	    { // calculate charged hadron v1 in TPC w.r.t. ZDC
	      StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	      if(!picoTrack) continue;

	      const TVector3 primMom = picoTrack->pMom(); // primary Momentum
	      const double pt  = primMom.Pt();
	      const double eta = primMom.PseudoRapidity();
	      const double phi = primMom.Phi(); // -pi to pi

	      if(mAnaCut->passTrkTpcFlowFull(picoTrack, primVtx))
	      {
		const double v1Zdc = TMath::Cos(1.0*(phi-Psi1ZdcFull))/mZdcEpManager->getZdcFullEpResVal(cent9);
		mZdcEpManager->fillZdcFullEpV1(eta, pt, v1Zdc, refWgt);
	      }
	    }
	  }

	  TVector2 vQ1EpdSideEast = mEpdEpManager->getQ1VecSideShiftEast(); // get Q1Vector from EPD
	  TVector2 vQ1EpdSideWest = mEpdEpManager->getQ1VecSideShiftWest();
	  TVector2 vQ1EpdSideFull = mEpdEpManager->getQ1VecSideShiftFull();
	  if( mAnaCut->passQVecEpdSide(vQ1EpdSideEast,vQ1EpdSideWest,vQ1EpdSideFull) && mEpdEpManager->getEpdSubEp1SideResVal(cent9) > 0.0 ) // EPD EP
	  { // charged hadron v1 from EPD
	    const double Psi1SideEpdEast = mEpdEpManager->getPsi1SideShiftEast();
	    const double Psi1SideEpdWest = mEpdEpManager->getPsi1SideShiftWest();
	    const double Psi1SideEpdFull = mEpdEpManager->getPsi1SideShiftFullCorr();
	    for(unsigned int iEpdHit = 0; iEpdHit < nEpdHits; ++iEpdHit)
	    { // calculate charged hadron v1 in EPD w.r.t. EPD
	      StPicoEpdHit *picoEpdHit = (StPicoEpdHit*)mPicoDst->epdHit(iEpdHit); 
	      if(!picoEpdHit) continue;

	      TVector3 EpdVector = mEpdEpManager->getEpdRanVec(picoEpdHit,primVtx);
	      const double eta = EpdVector.PseudoRapidity();
	      const double phi = EpdVector.Phi(); // -pi to pi

	      if( mAnaCut->passHitEpdFlowEast(picoEpdHit) ) // negative eta
	      {
		const double v1Epd = TMath::Cos(1.0*(phi-Psi1SideEpdWest))/mEpdEpManager->getEpdSubEp1SideResVal(cent9);
		mEpdEpManager->fillEpdSubEpSideV1(eta, v1Epd, refWgt);
	      }
	      if( mAnaCut->passHitEpdFlowWest(picoEpdHit) ) // positive eta
	      {
		const double v1Epd = TMath::Cos(1.0*(phi-Psi1SideEpdEast))/mEpdEpManager->getEpdSubEp1SideResVal(cent9);
		mEpdEpManager->fillEpdSubEpSideV1(eta, v1Epd, refWgt);
	      }
	    }
	    for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack)	  
	    { // calculate charged hadron v1 in TPC w.r.t. EPD
	      StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	      if(!picoTrack) continue;

	      const TVector3 primMom = picoTrack->pMom(); // primary Momentum
	      const double pt  = primMom.Pt();
	      const double eta = primMom.PseudoRapidity();
	      const double phi = primMom.Phi(); // -pi to pi

	      // keep the same pt cuts as charged v1 from ZDC
	      if(mAnaCut->passTrkTpcFlowFull(picoTrack, primVtx) && pt > 0.15 && pt < 2.0 ) // positive eta
	      { // use full EPD EP
		const double v1Epd = TMath::Cos(1.0*(phi-Psi1SideEpdFull))/mEpdEpManager->getEpdFullEp1SideResVal(cent9);
		mEpdEpManager->fillEpdSubEpSideV1(eta, v1Epd, refWgt);
	      }
	      // if(mAnaCut->passTrkTpcFlowEast(picoTrack, primVtx) && pt > 0.15 && pt < 2.0) // negative eta
	      // { // correlate track from East to EP from West EPD
		// const double v1Epd = TMath::Cos(1.0*(phi-Psi1SideEpdWest))/mEpdEpManager->getEpdSubEp1SideResVal(cent9);
		// mEpdEpManager->fillEpdSubEpSideV1(eta, v1Epd, refWgt);
	      // }
	      // if(mAnaCut->passTrkTpcFlowWest(picoTrack, primVtx) && pt > 0.15 && pt < 2.0 ) // positive eta
	      // { // correlate track from West to EP from East EPD
		// const double v1Epd = TMath::Cos(1.0*(phi-Psi1SideEpdEast))/mEpdEpManager->getEpdSubEp1SideResVal(cent9);
		// mEpdEpManager->fillEpdSubEpSideV1(eta, v1Epd, refWgt);
	      // }
	    }
	  }
	}

	const int numTrkReCtrEast = mTpcEpManager->getNumTrkReCtrEast(); // TPC EP
	const int numTrkReCtrWest = mTpcEpManager->getNumTrkReCtrWest();
	if(mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest))
	{ // charged hadron v2 and v3 from TPC
	  const TVector2 vQ1TpcEast = mTpcEpManager->getQ1VecReCtrEast();
	  const TVector2 vQ1TpcWest = mTpcEpManager->getQ1VecReCtrWest();
	  const TVector2 vQ2TpcEast = mTpcEpManager->getQ2VecReCtrEast();
	  const TVector2 vQ2TpcWest = mTpcEpManager->getQ2VecReCtrWest();
	  const TVector2 vQ3TpcEast = mTpcEpManager->getQ3VecReCtrEast();
	  const TVector2 vQ3TpcWest = mTpcEpManager->getQ3VecReCtrWest();

	  const double Psi1TpcEast = mTpcEpManager->getPsi1ShiftEast(vQ1TpcEast); // -pi to pi
	  const double Psi1TpcWest = mTpcEpManager->getPsi1ShiftWest(vQ1TpcWest);
	  const double Psi2TpcEast = mTpcEpManager->getPsi2ShiftEast(vQ2TpcEast); // -pi/2 to pi/2
	  const double Psi2TpcWest = mTpcEpManager->getPsi2ShiftWest(vQ2TpcWest);
	  const double Psi3TpcEast = mTpcEpManager->getPsi3ShiftEast(vQ3TpcEast); // -pi/3 to pi/3
	  const double Psi3TpcWest = mTpcEpManager->getPsi3ShiftWest(vQ3TpcWest);

	  for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack)	  
	  { // calculate charged charged hadron v2 & v3 int TPC w.r.t TPC
	    StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	    if(!picoTrack) continue;

	    const TVector3 primMom = picoTrack->pMom(); // primary Momentum
	    const double pt  = primMom.Pt();
	    // const double eta = primMom.PseudoRapidity();
	    const double phi = primMom.Phi(); // -pi to pi

	    if(mAnaCut->passTrkTpcFlowEast(picoTrack, primVtx)) // negative eta
	    { // correlate track from East to EP from West
	      if(mTpcEpManager->getTpcSubEp1ResVal(cent9) > 0.0)
	      {
		const double v1Tpc = TMath::Cos(1.0*(phi-Psi1TpcWest))/mTpcEpManager->getTpcSubEp1ResVal(cent9);
		mTpcEpManager->fillTpcSubEpV1(pt, v1Tpc, refWgt);
	      }
	      if(mTpcEpManager->getTpcSubEp2ResVal(cent9) > 0.0)
	      {
		const double v2Tpc = TMath::Cos(2.0*(phi-Psi2TpcWest))/mTpcEpManager->getTpcSubEp2ResVal(cent9);
		mTpcEpManager->fillTpcSubEpV2(pt, v2Tpc, refWgt);
	      }
	      if(mTpcEpManager->getTpcSubEp3ResVal(cent9) > 0.0)
	      {
		const double v3Tpc = TMath::Cos(3.0*(phi-Psi3TpcWest))/mTpcEpManager->getTpcSubEp3ResVal(cent9);
		mTpcEpManager->fillTpcSubEpV3(pt, v3Tpc, refWgt);
	      }
	    }
	    if(mAnaCut->passTrkTpcFlowWest(picoTrack, primVtx)) // positive eta
	    { // correlate track from West to EP from East
	      if(mTpcEpManager->getTpcSubEp1ResVal(cent9) > 0.0)
	      {
		const double v1Tpc = TMath::Cos(1.0*(phi-Psi1TpcEast))/mTpcEpManager->getTpcSubEp1ResVal(cent9);
		mTpcEpManager->fillTpcSubEpV1(pt, v1Tpc, refWgt);
	      }
	      if(mTpcEpManager->getTpcSubEp2ResVal(cent9) > 0.0)
	      {
		const double v2Tpc = TMath::Cos(2.0*(phi-Psi2TpcEast))/mTpcEpManager->getTpcSubEp2ResVal(cent9);
		mTpcEpManager->fillTpcSubEpV2(pt, v2Tpc, refWgt);
	      }
	      if(mTpcEpManager->getTpcSubEp3ResVal(cent9) > 0.0)
	      {
		const double v3Tpc = TMath::Cos(3.0*(phi-Psi3TpcEast))/mTpcEpManager->getTpcSubEp3ResVal(cent9);
		mTpcEpManager->fillTpcSubEpV3(pt, v3Tpc, refWgt);
	      }
	    }
	  }
	}

	if(mAnaCut->isFxt3p85GeV_2018())
	{
	  TVector2 vQ1EpdGrpEast[mNumRingsGrps], vQ1EpdGrpWest[mNumRingsGrps], vQ1EpdGrpFull[mNumRingsGrps]; 
	  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	  {
	    vQ1EpdGrpEast[iGrp] = mEpdEpManager->getQ1VecGrpShiftTrkAveEast(iGrp); // get Q1Vector from EPD
	    vQ1EpdGrpWest[iGrp] = mEpdEpManager->getQ1VecGrpShiftTrkAveWest(iGrp);
	    vQ1EpdGrpFull[iGrp] = mEpdEpManager->getQ1VecGrpShiftTrkAveFull(iGrp);
	  }

	  if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[0],vQ1EpdGrpWest[0],vQ1EpdGrpFull[0],0) &&
	      mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[1],vQ1EpdGrpWest[1],vQ1EpdGrpFull[1],1) )
	  // if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[0],vQ1EpdGrpWest[0],vQ1EpdGrpFull[0],0) &&
	  //     mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[1],vQ1EpdGrpWest[1],vQ1EpdGrpFull[1],1) &&
	  //     mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest) )
	  { // 3-sub events method
	    const double Psi1EpdGrp0 = mEpdEpManager->getPsi1GrpShiftTrkAveEast(0); // Psi1 from EPD Grp 0 East 
	    for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack)	  
	    { // calculate charged charged hadron v2 & v3 int TPC w.r.t TPC
	      StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	      if(!picoTrack) continue;

	      const TVector3 primMom = picoTrack->pMom(); // primary Momentum
	      const double phi       = primMom.Phi(); // -pi to pi
	      const double pMag      = primMom.Mag();
	      const double pT        = primMom.Pt();
	      const double deuteronZ = mAnaUtils->calcNSigmaZ(1.0,anaUtils::mMassDeuteron,pMag,picoTrack->dEdx()); // assume every track is deutron
	      const double mass2     = mAnaUtils->getPrimMass2(mPicoDst,iTrack);

	      if(mAnaCut->passTrkTpcFlow(picoTrack, primVtx) && mMixEpManager->getMixSubEp1Res1Val(cent9,0) > 0.0)
	      {
		const double etaLab = primMom.PseudoRapidity();
		const double v1Epd  = TMath::Cos(1.0*(phi-Psi1EpdGrp0))/mMixEpManager->getMixSubEp1Res1Val(cent9,0);
		if(mAnaCut->passTrkProFlow(picoTrack))
		{
		  if( pT > 0.4 && pT < 1.0 && cent9 >=4 && cent9 <=6 )
		  { // 10-40% with the same cuts as PLB827,136941
		    const double yLab   = mAnaUtils->getRapidityLab(picoTrack, 2212); // y_proton in lab frame
		    const double yCms   = mAnaUtils->getRapidityCMS(yLab);
		    const double p_eff  = 1.0;
		    if(p_eff > 0.0 && yCms > 0.0) mMixEpManager->fillMixSubEpProV1(yCms, v1Epd, refWgt, p_eff);
		  }
		}
		if(mAnaCut->passTrkDeuFlow(pMag,deuteronZ,mass2))
		{
		  if( pT > 0.8 && pT < 2.0 && cent9 >=4 && cent9 <=6 )
		  { // 10-40% with the same cuts as PLB827,136941
		    const double yLab   = mAnaUtils->getRapidityLab(picoTrack, 1234); // y_deuteron in lab frame
		    const double yCms   = mAnaUtils->getRapidityCMS(yLab);
		    const double d_eff  = mMixEpManager->calcDeuEfficiency(pT,pMag,etaLab,yCms);
		    // const double d_eff  = 1.0;
		    if(d_eff > 0.0 && yCms > 0.0) mMixEpManager->fillMixSubEpDeuV1(yCms, v1Epd, refWgt, d_eff);
		  }
		}
	      }
	    }
	  }
	}
      }

      mZdcEpManager->clearZdcEpManager();
      mEpdEpManager->clearEpdEpManager();
      mTpcEpManager->clearTpcEpManager();
      mMixEpManager->clearMixEpManager();
    }
  }

  return kStOK;
}
