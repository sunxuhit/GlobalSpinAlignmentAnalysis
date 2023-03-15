#include <algorithm>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

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
#include "StRoot/StEventPlaneMaker/StEventPlaneMaker.h"

ClassImp(StEventPlaneMaker)

//-----------------------------------------------------------------------------
StEventPlaneMaker::StEventPlaneMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int mode, const int beamType) : StMaker(name), mMode(mode), mType(beamType)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;
  // mMode = Mode;
  // mType = beamType;

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
  mAnaCut       = new StAnalysisCut(mType);
  mAnaUtils     = new StAnalysisUtils(mType);
  mAnaUtils->initRunIndex(); // initialize std::map for run index

  if(!mRefMultCorr)
  {
    if( mAnaCut->isIsobar() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr_Isobar();
    if( mAnaCut->isFxt() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }

  if(mMode == 0)
  { // fill Gain Correction Factors for ZDC and phi Weight for EPD
    file_mOutPutGainCorr = new TFile(str_mOutPutGainCorr.c_str(),"RECREATE");
    mZdcEpManager->initZdcGain(); // ZDC

    mEpdEpManager->initEpdPhiWgt(); // EPD
    mEpdEpManager->initEpdSubEpRaw();
  }
  if(mMode == 1)
  { // fill ReCenter Correction Parameters for ZDC & EPD & TPC Sub EP
    file_mOutPutReCenterPar = new TFile(str_mOutPutReCenterPar.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->initZdcReCtr();
    mZdcEpManager->initZdcSubEpRaw();

    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->initEpdReCtr();
    mEpdEpManager->initEpdSubEpWgt();

    mTpcEpManager->initTpcReCtr(); // TPC
    mTpcEpManager->initTpcSubEpRaw();
  }
  if(mMode == 2)
  { // fill Shift Correction Parameters for ZDC & EPD & TPC Sub EP
    file_mOutPutShiftPar = new TFile(str_mOutPutShiftPar.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->readZdcReCtr();
    mZdcEpManager->initZdcShift();
    mZdcEpManager->initZdcSubEpReCtr();

    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->readEpdReCtr();
    mEpdEpManager->initEpdShift();
    mEpdEpManager->initEpdSubEpReCtr();

    mTpcEpManager->readTpcReCtr(); // TPC
    mTpcEpManager->initTpcShift();
    mTpcEpManager->initTpcSubEpReCtr();
  }
  if(mMode == 3)
  { // fill Shift Correction Parameters for ZDC & EPD Full EP
    file_mOutPutShiftPar = new TFile(str_mOutPutShiftPar.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->readZdcReCtr();
    mZdcEpManager->readZdcShift();
    mZdcEpManager->initZdcShiftFull();
    mZdcEpManager->initZdcSubEpShift();

    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->readEpdReCtr();
    mEpdEpManager->readEpdShift();
    mEpdEpManager->initEpdShiftFull();
    mEpdEpManager->initEpdSubEpShift();
  }
  if(mMode == 4)
  { // fill Event Plane Resolution for ZDC & EPD & TPC Sub EP
    file_mOutPutResolution = new TFile(str_mOutPutResolution.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->readZdcReCtr();
    mZdcEpManager->readZdcShift();
    mZdcEpManager->readZdcShiftFull();
    mZdcEpManager->initZdcResolution();
    mZdcEpManager->initZdcSubEpShift();
    mZdcEpManager->initZdcFullEpShift();

    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->readEpdReCtr();
    mEpdEpManager->readEpdShift();
    mEpdEpManager->readEpdShiftFull();
    mEpdEpManager->initEpdResolution();
    mEpdEpManager->initEpdSubEpShift();
    mEpdEpManager->initEpdFullEpShift();

    mTpcEpManager->readTpcReCtr(); // TPC
    mTpcEpManager->readTpcShift();
    mTpcEpManager->initTpcResolution();
    mTpcEpManager->initTpcSubEpShift();
  }
  if(mMode == 5)
  { // calculate charged hadron v1 from ZDC & EPD and charged hadron v2 and v3 from TPC
    file_mOutPutFlow = new TFile(str_mOutPutFlow.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->readZdcReCtr();
    mZdcEpManager->readZdcShift();
    mZdcEpManager->readZdcShiftFull();
    mZdcEpManager->readZdcResolution();
    mZdcEpManager->initZdcFullEpFlow();

    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->readEpdReCtr();
    mEpdEpManager->readEpdShift();
    mEpdEpManager->readEpdShiftFull();
    mEpdEpManager->readEpdResolution();
    mEpdEpManager->initEpdSubEpFlow();

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
      mZdcEpManager->writeZdcGain(); // ZDC

      mEpdEpManager->writeEpdPhiWgt(); // EPD
      mEpdEpManager->writeEpdSubEpRaw();
      file_mOutPutGainCorr->Close();
    }
  }
  if(mMode == 1)
  {
    if(str_mOutPutReCenterPar != "")
    {
      file_mOutPutReCenterPar->cd();
      mZdcEpManager->writeZdcReCtr(); // ZDC
      mZdcEpManager->writeZdcSubEpRaw();

      mEpdEpManager->writeEpdReCtr(); // EPD
      mEpdEpManager->writeEpdSubEpWgt();

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
      mZdcEpManager->writeZdcShift(); // ZDC
      mZdcEpManager->writeZdcSubEpReCtr();

      mEpdEpManager->writeEpdShift(); // EPD
      mEpdEpManager->writeEpdSubEpReCtr();

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
      mZdcEpManager->writeZdcShiftFull(); // ZDC
      mZdcEpManager->writeZdcSubEpShift();

      mEpdEpManager->writeEpdShiftFull(); // EPD
      mEpdEpManager->writeEpdSubEpShift(); // EPD
      file_mOutPutShiftPar->Close();
    }
  }
  if(mMode == 4)
  {
    if(str_mOutPutResolution != "")
    {
      file_mOutPutResolution->cd();
      mZdcEpManager->writeZdcResolution(); // ZDC
      mZdcEpManager->writeZdcSubEpShift();
      mZdcEpManager->writeZdcFullEpShift();

      mEpdEpManager->writeEpdResolution(); // EPD
      mEpdEpManager->writeEpdSubEpShift();
      mEpdEpManager->writeEpdFullEpShift();

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
      mZdcEpManager->writeZdcFullEpFlow(); // ZDC
      mEpdEpManager->writeEpdSubEpFlow(); // EPD
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

    bool isPileUpEventStAnalysisCut = mAnaCut->isPileUpEvent((double)refMult, (double)nBTofMatch,vz); // alway return false for Isobar & Fxt (for now)
    bool isPileUpEventStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut((double)refMult, (double)nBTofMatch,vz); // valid for Isobar || NOT for Fxt
    bool isPileUpEvent = isPileUpEventStAnalysisCut || isPileUpEventStRefMultCorr;
    // cout << "isPileUpEvent = " << isPileUpEvent << ", isPileUpEventStAnalysisCut = " << isPileUpEventStAnalysisCut << ", isPileUpEventStRefMultCorr = " << isPileUpEventStRefMultCorr << endl;

    if(!isPileUpEvent && mAnaCut->isGoodCent9(cent9) && mAnaCut->passEventCut(mPicoEvent))
    { // apply Event Cuts for anlaysis 
      mZdcEpManager->initZdcEpManager(cent9,runIndex,vzBin); // initialize ZDC EP Manager
      mEpdEpManager->initEpdEpManager(cent9,runIndex,vzBin); // initialize EPD EP Manager
      mTpcEpManager->initTpcEpManager(cent9,runIndex,vzBin); // initialize TPC EP Manager
      if(mMode == 0)
      { // fill Gain Correction Factors for ZDC
	for(int iSlat = 0; iSlat < 8; ++iSlat) // read in raw ADC value from ZDC-SMD
	{
	  mZdcEpManager->setZdcSmd(0,0,iSlat,mPicoEvent->ZdcSmdEastVertical(iSlat));
	  mZdcEpManager->setZdcSmd(0,1,iSlat,mPicoEvent->ZdcSmdEastHorizontal(iSlat));
	  mZdcEpManager->setZdcSmd(1,0,iSlat,mPicoEvent->ZdcSmdWestVertical(iSlat));
	  mZdcEpManager->setZdcSmd(1,1,iSlat,mPicoEvent->ZdcSmdWestHorizontal(iSlat));
	}
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

	for(unsigned int iEpdHit = 0; iEpdHit < nEpdHits; ++iEpdHit)
	{ // fill phi weight for EPD
	  StPicoEpdHit *picoEpdHit = (StPicoEpdHit*)mPicoDst->epdHit(iEpdHit); 
	  if(!picoEpdHit) continue;

	  if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	  {
	    mEpdEpManager->addHitRawEast(picoEpdHit, primVtx);
	    mEpdEpManager->fillEpdPhiWgtEast(picoEpdHit);
	  }
	  if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	  {
	    mEpdEpManager->addHitRawWest(picoEpdHit, primVtx);
	    mEpdEpManager->fillEpdPhiWgtWest(picoEpdHit);
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

	TVector2 vQ1EpdGrpEast[mNumGroups], vQ1EpdGrpWest[mNumGroups], vQ1EpdGrpFull[mNumGroups];
	for(int iGrp = 0; iGrp < mNumGroups; ++iGrp)
	{
	  vQ1EpdGrpEast[iGrp] = mEpdEpManager->getQ1VecGrpRawEast(iGrp); // get Q1Vector from EPD groups
	  vQ1EpdGrpWest[iGrp] = mEpdEpManager->getQ1VecGrpRawWest(iGrp);
	  vQ1EpdGrpFull[iGrp] = mEpdEpManager->getQ1VecGrpRawFull(iGrp);
	  if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[iGrp],vQ1EpdGrpWest[iGrp],vQ1EpdGrpFull[iGrp],iGrp)) // EPD EP
	  {
	    const double Psi1GrpRawEast = mEpdEpManager->getPsi1GrpRawEast(iGrp);
	    const double Psi1GrpRawWest = mEpdEpManager->getPsi1GrpRawWest(iGrp);
	    const double Psi1GrpRawFull = mEpdEpManager->getPsi1GrpRawFull(iGrp);
	    mEpdEpManager->fillEpdSubEpGrpRaw(Psi1GrpRawEast,Psi1GrpRawWest,Psi1GrpRawFull,iGrp);
	  }
	}
      }
      if(mMode > 0)
      {
	for(int iSlat = 0; iSlat < 8; ++iSlat) // read and apply gain correction to raw ADC value for ZDC
	{
	  mZdcEpManager->setZdcSmdGainCorr(0,0,iSlat,mPicoEvent->ZdcSmdEastVertical(iSlat));
	  mZdcEpManager->setZdcSmdGainCorr(0,1,iSlat,mPicoEvent->ZdcSmdEastHorizontal(iSlat));
	  mZdcEpManager->setZdcSmdGainCorr(1,0,iSlat,mPicoEvent->ZdcSmdWestVertical(iSlat));
	  mZdcEpManager->setZdcSmdGainCorr(1,1,iSlat,mPicoEvent->ZdcSmdWestHorizontal(iSlat));
	}

	for(unsigned int iEpdHit = 0; iEpdHit < nEpdHits; ++iEpdHit) // read and apply phi weight & eta weight (if avaliable) for EPD
	{ // fill phi weight for EPD
	  StPicoEpdHit *picoEpdHit = (StPicoEpdHit*)mPicoDst->epdHit(iEpdHit); 
	  if(!picoEpdHit) continue;

	  if(mMode == 1) // calculate phi and eta weighted (if avaliable) Q1Vector from EPD
	  {
	    if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	    {
	      mEpdEpManager->addHitWgtEast(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdSideReCtrEast(picoEpdHit, primVtx); // fill EPD ReCenter Parameters East
	      mEpdEpManager->fillEpdGrpReCtrEast(picoEpdHit, primVtx);
	    }
	    if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	    {
	      mEpdEpManager->addHitWgtWest(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdSideReCtrWest(picoEpdHit, primVtx); // fill EPD ReCenter Parameters West
	      mEpdEpManager->fillEpdGrpReCtrWest(picoEpdHit, primVtx);
	    }
	  }
	  else // calculate phi & eta weighted (if avaliable) and recentered Q1Vector from EPD
	  {
	    if( mAnaCut->passHitEpdEpEast(picoEpdHit) ) // negative eta
	    {
	      mEpdEpManager->addHitReCtrEast(picoEpdHit, primVtx);
	    }
	    if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	    {
	      mEpdEpManager->addHitReCtrWest(picoEpdHit, primVtx);
	    }
	  }
	}

	for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack) // calculate Q2Vector and Q3Vector from TPC
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
	  if(mMode == 2 || mMode == 4) // calculate recentered Q1Vector, Q2Vector and Q3Vector from TPC
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

	if(mMode == 1) // fill recenter correction parameter for ZDC & EPD & TPC Sub EP
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
	    const double Psi1SideWgtEast = mEpdEpManager->getPsi1SideWgtEast(); // EPD EP
	    const double Psi1SideWgtWest = mEpdEpManager->getPsi1SideWgtWest();
	    const double Psi1SideWgtFull = mEpdEpManager->getPsi1SideWgtFull();
	    mEpdEpManager->fillEpdSubEpSideWgt(Psi1SideWgtEast, Psi1SideWgtWest, Psi1SideWgtFull);
	  }

	  TVector2 vQ1EpdGrpEast[mNumGroups], vQ1EpdGrpWest[mNumGroups], vQ1EpdGrpFull[mNumGroups];
	  for(int iGrp = 0; iGrp < mNumGroups; ++iGrp)
	  {
	    vQ1EpdGrpEast[iGrp] = mEpdEpManager->getQ1VecGrpWgtEast(iGrp); // get Q1Vector from EPD groups
	    vQ1EpdGrpWest[iGrp] = mEpdEpManager->getQ1VecGrpWgtWest(iGrp);
	    vQ1EpdGrpFull[iGrp] = mEpdEpManager->getQ1VecGrpWgtFull(iGrp);
	    if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[iGrp],vQ1EpdGrpWest[iGrp],vQ1EpdGrpFull[iGrp],iGrp) ) // EPD EP
	    {
	      const double Psi1GrpWgtEast = mEpdEpManager->getPsi1GrpWgtEast(iGrp); // EPD EP
	      const double Psi1GrpWgtWest = mEpdEpManager->getPsi1GrpWgtWest(iGrp);
	      const double Psi1GrpWgtFull = mEpdEpManager->getPsi1GrpWgtFull(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpWgt(Psi1GrpWgtEast, Psi1GrpWgtWest, Psi1GrpWgtFull, iGrp);
	    }
	  }

	  const int numTrackRawEast = mTpcEpManager->getNumTrkRawEast(); // TPC EP
	  const int numTrackRawWest = mTpcEpManager->getNumTrkRawWest();
	  if(mAnaCut->passNumTrkTpcSubEpRaw(numTrackRawEast, numTrackRawWest))
	  {
	    const double Psi1RawEast = mTpcEpManager->getPsi1RawEast();
	    const double Psi1RawWest = mTpcEpManager->getPsi1RawWest();
	    const double Psi1RawFull = mTpcEpManager->getPsi1RawFull();
	    const double Psi2RawEast = mTpcEpManager->getPsi2RawEast();
	    const double Psi2RawWest = mTpcEpManager->getPsi2RawWest();
	    const double Psi2RawFull = mTpcEpManager->getPsi2RawFull();
	    const double Psi3RawEast = mTpcEpManager->getPsi3RawEast();
	    const double Psi3RawWest = mTpcEpManager->getPsi3RawWest();
	    const double Psi3RawFull = mTpcEpManager->getPsi3RawFull();
	    mTpcEpManager->fillTpcSubEpRaw(Psi1RawEast, Psi1RawWest, Psi1RawFull, Psi2RawEast, Psi2RawWest, Psi2RawFull, Psi3RawEast, Psi3RawWest, Psi3RawFull);
	  }
	}
	if(mMode == 2) // fill shift correction parameter for ZDC & EPD & TPC Sub EP
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

	  TVector2 vQ1EpdGrpEast[mNumGroups], vQ1EpdGrpWest[mNumGroups], vQ1EpdGrpFull[mNumGroups];
	  for(int iGrp = 0; iGrp < mNumGroups; ++iGrp)
	  {
	    vQ1EpdGrpEast[iGrp] = mEpdEpManager->getQ1VecGrpReCtrEast(iGrp); // get Q1Vector from EPD groups
	    vQ1EpdGrpWest[iGrp] = mEpdEpManager->getQ1VecGrpReCtrWest(iGrp);
	    vQ1EpdGrpFull[iGrp] = mEpdEpManager->getQ1VecGrpReCtrFull(iGrp);
	    if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[iGrp],vQ1EpdGrpWest[iGrp],vQ1EpdGrpFull[iGrp],iGrp) ) // EPD EP
	    {
	      const double Psi1GrpReCtrEast = mEpdEpManager->getPsi1GrpReCtrEast(iGrp); // EPD EP
	      const double Psi1GrpReCtrWest = mEpdEpManager->getPsi1GrpReCtrWest(iGrp);
	      const double Psi1GrpReCtrFull = mEpdEpManager->getPsi1GrpReCtrFull(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpReCtr(Psi1Grp0ReCtrEast, Psi1Grp0ReCtrWest, Psi1Grp0ReCtrFull, Psi1Grp1ReCtrEast, Psi1Grp1ReCtrWest, Psi1Grp1ReCtrFull);
	      mEpdEpManager->fillEpdGrpShiftEast(iGrp);
	      mEpdEpManager->fillEpdGrpShiftWest(iGrp);
	    }
	  }

	  const int numTrkReCtrEast = mTpcEpManager->getNumTrkReCtrEast(); // TPC EP
	  const int numTrkReCtrWest = mTpcEpManager->getNumTrkReCtrWest();
	  if(mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest))
	  {
	    const double Psi1ReCtrEast = mTpcEpManager->getPsi1ReCtrEast();
	    const double Psi1ReCtrWest = mTpcEpManager->getPsi1ReCtrWest();
	    const double Psi1ReCtrFull = mTpcEpManager->getPsi1ReCtrFull();
	    const double Psi2ReCtrEast = mTpcEpManager->getPsi2ReCtrEast();
	    const double Psi2ReCtrWest = mTpcEpManager->getPsi2ReCtrWest();
	    const double Psi2ReCtrFull = mTpcEpManager->getPsi2ReCtrFull();
	    const double Psi3ReCtrEast = mTpcEpManager->getPsi3ReCtrEast();
	    const double Psi3ReCtrWest = mTpcEpManager->getPsi3ReCtrWest();
	    const double Psi3ReCtrFull = mTpcEpManager->getPsi3ReCtrFull();
	    mTpcEpManager->fillTpcSubEpReCtr(Psi1ReCtrEast, Psi1ReCtrWest, Psi1ReCtrFull, Psi2ReCtrEast, Psi2ReCtrWest, Psi2ReCtrFull, Psi3ReCtrEast, Psi3ReCtrWest, Psi3ReCtrFull);
	    mTpcEpManager->fillTpcShiftEast();
	    mTpcEpManager->fillTpcShiftWest();
	    mTpcEpManager->fillTpcShiftFull();
	  }
	}
	if(mMode == 3) // fill shift correction parameter for ZDC & EPD Full EP
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
	    const double Psi1SideShiftEast = mEpdEpManager->getPsi1SideShiftEast(); // EPD EP
	    const double Psi1SideShiftWest = mEpdEpManager->getPsi1SideShiftWest();
	    const double Psi1SideShiftFull = mEpdEpManager->getPsi1SideShiftFull();
	    mEpdEpManager->fillEpdSubEpSideShift(Psi1SideShiftEast, Psi1SideShiftWest, Psi1SideShiftFull);
	    mEpdEpManager->fillEpdSideShiftFull();
	  }

	  TVector2 vQ1EpdGrpEast[mNumGroups], vQ1EpdGrpWest[mNumGroups], vQ1EpdGrpFull[mNumGroups];
	  for(int iGrp = 0; iGrp < mNumGroups; ++iGrp)
	  {
	    vQ1EpdGrpEast[iGrp] = mEpdEpManager->getQ1VecGrpShiftEast(iGrp); // get Q1Vector from EPD groups
	    vQ1EpdGrpWest[iGrp] = mEpdEpManager->getQ1VecGrpShiftWest(iGrp);
	    vQ1EpdGrpFull[iGrp] = mEpdEpManager->getQ1VecGrpShiftFull(iGrp);
	    if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[iGrp],vQ1EpdGrpWest[iGrp],vQ1EpdGrpFull[iGrp],iGrp) ) // EPD EP
	    {
	      const double Psi1GrpShiftEast = mEpdEpManager->getPsi1GrpShiftEast(iGrp); // EPD EP
	      const double Psi1GrpShiftWest = mEpdEpManager->getPsi1GrpShiftWest(iGrp);
	      const double Psi1GrpShiftFull = mEpdEpManager->getPsi1GrpShiftFull(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpShift(Psi1GrpShiftEast, Psi1GrpShiftWest, Psi1GrpShiftFull, iGrp);
	      mEpdEpManager->fillEpdGrpShiftFull(iGrp);
	    }
	  }
	}
	if(mMode == 4) // fill event plane resolution for ZDC-SMD & EPD & TPC Sub EP
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

	  TVector2 vQ1EpdSideEast     = mEpdEpManager->getQ1VecSideShiftEast(); // get Q1Vector from EPD
	  TVector2 vQ1EpdSideWest     = mEpdEpManager->getQ1VecSideShiftWest();
	  TVector2 vQ1EpdSideFull     = mEpdEpManager->getQ1VecSideShiftFull();
	  // TVector2 vQ1EpdSideFullCorr = mEpdEpManager->getQ1VecSideShiftFullCorr();
	  if( mAnaCut->passQVecEpdSide(vQ1EpdSideEast,vQ1EpdSideWest,vQ1EpdSideFull) ) // EPD EP
	  {
	    const double Psi1SideShiftEast     = mEpdEpManager->getPsi1SideShiftEast(); // EPD EP
	    const double Psi1SideShiftWest     = mEpdEpManager->getPsi1SideShiftWest();
	    const double Psi1SideShiftFull     = mEpdEpManager->getPsi1SideShiftFull();
	    const double Psi1SideShiftFullCorr = mEpdEpManager->getPsi1SideShiftFullCorr();
	    mEpdEpManager->fillEpdSubEpSideShift(Psi1SideShiftEast, Psi1SideShiftWest, Psi1SideShiftFull);
	    mEpdEpManager->fillEpdFullEpSideShift(Psi1SideShiftFullCorr);
	    mEpdEpManager->fillEpdSideResolution(Psi1SideShiftEast, Psi1SideShiftWest);
	  }

	  TVector2 vQ1EpdGrpEast[mNumGroups], vQ1EpdGrpWest[mNumGroups], vQ1EpdGrpFull[mNumGroups];
	  for(int iGrp = 0; iGrp < mNumGroups; ++iGrp)
	  {
	    TVector2 vQ1EpdGrpEast[iGrp]     = mEpdEpManager->getQ1VecGrpShiftEast(iGrp); // get Q1Vector from EPD
	    TVector2 vQ1EpdGrpWest[iGrp]     = mEpdEpManager->getQ1VecGrpShiftWest(iGrp);
	    TVector2 vQ1EpdGrpFull[iGrp]     = mEpdEpManager->getQ1VecGrpShiftFull(iGrp);
	    TVector2 vQ1EpdGrpFullCorr[iGrp] = mEpdEpManager->getQ1VecGrpShiftFullCorr(iGrp);
	    if( mAnaCut->passQVecEpdGrp(vQ1EpdGrpEast[iGrp],vQ1EpdGrpWest[iGrp],vQ1EpdGrpFull[iGrp],iGrp) ) // EPD EP
	    {
	      const double Psi1GrpShiftEast     = mEpdEpManager->getPsi1GrpShiftEast(iGrp); // EPD EP
	      const double Psi1GrpShiftWest     = mEpdEpManager->getPsi1GrpShiftWest(iGrp);
	      const double Psi1GrpShiftFull     = mEpdEpManager->getPsi1GrpShiftFull(iGrp);
	      const double Psi1GrpShiftFullCorr = mEpdEpManager->getPsi1GrpShiftFullCorr(iGrp);
	      mEpdEpManager->fillEpdSubEpGrpShift(Psi1GrpShiftEast, Psi1GrpShiftWest, Psi1GrpShiftFull, iGrp);
	      mEpdEpManager->fillEpdFullEpGrpShift(Psi1GrpShiftFullCorr,iGrp);
	      mEpdEpManager->fillEpdGrpResolution(Psi1GrpShiftEast, Psi1GrpShiftWest, iGrp);
	    }
	  }

	  const int numTrkReCtrEast = mTpcEpManager->getNumTrkReCtrEast(); // TPC EP
	  const int numTrkReCtrWest = mTpcEpManager->getNumTrkReCtrWest();
	  if(mAnaCut->passNumTrkTpcSubEpReCtr(numTrkReCtrEast, numTrkReCtrWest))
	  {
	    const TVector2 vQ1TpcEast = mTpcEpManager->getQ1VecReCtrEast();
	    const TVector2 vQ1TpcWest = mTpcEpManager->getQ1VecReCtrWest();
	    const TVector2 vQ1TpcFull = mTpcEpManager->getQ1VecReCtrFull();
	    const TVector2 vQ2TpcEast = mTpcEpManager->getQ2VecReCtrEast();
	    const TVector2 vQ2TpcWest = mTpcEpManager->getQ2VecReCtrWest();
	    const TVector2 vQ2TpcFull = mTpcEpManager->getQ2VecReCtrFull();
	    const TVector2 vQ3TpcEast = mTpcEpManager->getQ3VecReCtrEast();
	    const TVector2 vQ3TpcWest = mTpcEpManager->getQ3VecReCtrWest();
	    const TVector2 vQ3TpcFull = mTpcEpManager->getQ3VecReCtrFull();

	    const double Psi1ShiftEast = mTpcEpManager->getPsi1ShiftEast(vQ1TpcEast);
	    const double Psi1ShiftWest = mTpcEpManager->getPsi1ShiftWest(vQ1TpcWest);
	    const double Psi1ShiftFull = mTpcEpManager->getPsi1ShiftFull(vQ1TpcFull);
	    const double Psi2ShiftEast = mTpcEpManager->getPsi2ShiftEast(vQ2TpcEast);
	    const double Psi2ShiftWest = mTpcEpManager->getPsi2ShiftWest(vQ2TpcWest);
	    const double Psi2ShiftFull = mTpcEpManager->getPsi2ShiftFull(vQ2TpcFull);
	    const double Psi3ShiftEast = mTpcEpManager->getPsi3ShiftEast(vQ3TpcEast);
	    const double Psi3ShiftWest = mTpcEpManager->getPsi3ShiftWest(vQ3TpcWest);
	    const double Psi3ShiftFull = mTpcEpManager->getPsi3ShiftFull(vQ3TpcFull);

	    mTpcEpManager->fillTpcSubEpShift(Psi1ShiftEast, Psi1ShiftWest, Psi1ShiftFull, Psi2ShiftEast, Psi2ShiftWest, Psi2ShiftFull, Psi3ShiftEast, Psi3ShiftWest, Psi3ShiftFull);
	    mTpcEpManager->fillTpcResolution(Psi1ShiftEast, Psi1ShiftWest, Psi2ShiftEast, Psi2ShiftWest, Psi3ShiftEast, Psi3ShiftWest);
	  }
	}
	if(mMode == 5) // calculate charged hadron v1 from ZDC & EPD and charged hadron v2 and v3 from TPC
	{
	  TVector2 vQ1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 vQ1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 vQ1ZdcFull = mZdcEpManager->getQ1VecFull(vQ1ZdcEast,vQ1ZdcWest,mMode); // TVector2 vQ1ZdcFull = vQ1ZdcWest-vQ1ZdcEast;
	  if( mAnaCut->passQVecZdc(vQ1ZdcEast,vQ1ZdcWest,vQ1ZdcFull) ) // ZDC EP
	  { // charged hadron v1 from ZDC
	    const double Psi1ZdcFull = TMath::ATan2(vQ1ZdcFull.Y(),vQ1ZdcFull.X()); // -pi to pi
	    for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack)	  
	    { // calculate charged hadron v1 from ZDC and charged hadron v2 & v3 from TPC
	      StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	      if(!picoTrack) continue;

	      const TVector3 primMom = picoTrack->pMom(); // primary Momentum
	      const double pt  = primMom.Pt();
	      const double eta = primMom.PseudoRapidity();
	      const double phi = primMom.Phi(); // -pi to pi

	      if(mAnaCut->passTrkTpcFlowFull(picoTrack, primVtx) && mZdcEpManager->getZdcFullEpResVal(cent9) > 0.0)
	      {
		const double v1Zdc = TMath::Cos(1.0*(phi-Psi1ZdcFull))/mZdcEpManager->getZdcFullEpResVal(cent9);
		mZdcEpManager->fillZdcFullEpDFlow(eta, pt, v1Zdc, reweight);
	      }
	    }
	  }

	  TVector2 vQ1EpdSideEast = mEpdEpManager->getQ1VecSideShiftEast(); // get Q1Vector from EPD
	  TVector2 vQ1EpdSideWest = mEpdEpManager->getQ1VecSideShiftWest();
	  TVector2 vQ1EpdSideFull = mEpdEpManager->getQ1VecSideShiftFull();
	  if( mAnaCut->passQVecEpdSide(vQ1EpdSideEast,vQ1EpdSideWest,vQ1EpdSideFull) ) // EPD EP
	  {
	    const double Psi1SideEpdEast = mEpdEpManager->getPsi1SideShiftEast(); // EPD EP
	    const double Psi1SideEpdWest = mEpdEpManager->getPsi1SideShiftWest();
	    for(unsigned int iEpdHit = 0; iEpdHit < nEpdHits; ++iEpdHit)
	    { // charged hadron v1 from EPD
	      StPicoEpdHit *picoEpdHit = (StPicoEpdHit*)mPicoDst->epdHit(iEpdHit); 
	      if(!picoEpdHit) continue;

	      TVector3 EpdVector = mEpdEpManager->getEpdRanVec(picoEpdHit,primVtx);
	      const double eta = EpdVector.PseudoRapidity();
	      const double phi = EpdVector.Phi(); // -pi to pi

	      if( mAnaCut->passHitEpdFlowEast(picoEpdHit) && mEpdEpManager->getEpdSubEp1SideResVal(cent9) > 0.0 ) // negative eta
	      {
		const double v1Epd = TMath::Cos(1.0*(phi-Psi1SideEpdWest))/mEpdEpManager->getEpdSubEp1SideResVal(cent9);
		mEpdEpManager->fillEpdSubEpDFlow(eta, v1Epd, reweight);
	      }
	      if( mAnaCut->passHitEpdFlowWest(picoEpdHit) && mEpdEpManager->getEpdSubEp1SideResVal(cent9) > 0.0  ) // positive eta
	      {
		const double v1Epd = TMath::Cos(1.0*(phi-Psi1SideEpdEast))/mEpdEpManager->getEpdSubEp1SideResVal(cent9);
		mEpdEpManager->fillEpdSubEpDFlow(eta, v1Epd, reweight);
	      }
	    }
	    for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack)	  
	    { // calculate charged hadron v1 from ZDC and charged hadron v2 & v3 from TPC
	      StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	      if(!picoTrack) continue;

	      const TVector3 primMom = picoTrack->pMom(); // primary Momentum
	      const double pt  = primMom.Pt();
	      const double eta = primMom.PseudoRapidity();
	      const double phi = primMom.Phi(); // -pi to pi

	      if(mAnaCut->passTrkTpcFlowEast(picoTrack, primVtx) && mEpdEpManager->getEpdSubEp1SideResVal(cent9) > 0.0) // neg
	      { // correlate track from East to EP from West EPD
		const double v1Epd = TMath::Cos(1.0*(phi-Psi1SideEpdWest))/mEpdEpManager->getEpdSubEp1SideResVal(cent9);
		mEpdEpManager->fillEpdSubEpDFlow(eta, v1Epd, reweight);
	      }
	      if(mAnaCut->passTrkTpcFlowWest(picoTrack, primVtx) && mEpdEpManager->getEpdSubEp1SideResVal(cent9) > 0.0) // neg
	      { // correlate track from West to EP from East EPD
		const double v1Epd = TMath::Cos(1.0*(phi-Psi1SideEpdEast))/mEpdEpManager->getEpdSubEp1SideResVal(cent9);
		mEpdEpManager->fillEpdSubEpDFlow(eta, v1Epd, reweight);
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
	    { // calculate charged hadron v1 from ZDC and charged hadron v2 & v3 from TPC
	      StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	      if(!picoTrack) continue;

	      const TVector3 primMom = picoTrack->pMom(); // primary Momentum
	      const double pt  = primMom.Pt();
	      // const double eta = primMom.PseudoRapidity();
	      const double phi = primMom.Phi(); // -pi to pi

	      if(mAnaCut->passTrkTpcFlowEast(picoTrack, primVtx)) // neg
	      { // correlate track from East to EP from West
		if(mTpcEpManager->getTpcSubEp1ResVal(cent9) > 0.0)
		{
		  const double v1Tpc = TMath::Cos(1.0*(phi-Psi1TpcWest))/mTpcEpManager->getTpcSubEp1ResVal(cent9);
		  mTpcEpManager->fillTpcSubEpV1(pt, v1Tpc, reweight);
		}
		if(mTpcEpManager->getTpcSubEp2ResVal(cent9) > 0.0)
		{
		  const double v2Tpc = TMath::Cos(2.0*(phi-Psi2TpcWest))/mTpcEpManager->getTpcSubEp2ResVal(cent9);
		  mTpcEpManager->fillTpcSubEpV2(pt, v2Tpc, reweight);
		}
		if(mTpcEpManager->getTpcSubEp3ResVal(cent9) > 0.0)
		{
		  const double v3Tpc = TMath::Cos(3.0*(phi-Psi3TpcWest))/mTpcEpManager->getTpcSubEp3ResVal(cent9);
		  mTpcEpManager->fillTpcSubEpV3(pt, v3Tpc, reweight);
		}
	      }
	      if(mAnaCut->passTrkTpcFlowWest(picoTrack, primVtx)) // neg
	      { // correlate track from West to EP from East
		if(mTpcEpManager->getTpcSubEp1ResVal(cent9) > 0.0)
		{
		  const double v1Tpc = TMath::Cos(1.0*(phi-Psi1TpcEast))/mTpcEpManager->getTpcSubEp1ResVal(cent9);
		  mTpcEpManager->fillTpcSubEpV1(pt, v1Tpc, reweight);
		}
		if(mTpcEpManager->getTpcSubEp2ResVal(cent9) > 0.0)
		{
		  const double v2Tpc = TMath::Cos(2.0*(phi-Psi2TpcEast))/mTpcEpManager->getTpcSubEp2ResVal(cent9);
		  mTpcEpManager->fillTpcSubEpV2(pt, v2Tpc, reweight);
		}
		if(mTpcEpManager->getTpcSubEp3ResVal(cent9) > 0.0)
		{
		  const double v3Tpc = TMath::Cos(3.0*(phi-Psi3TpcEast))/mTpcEpManager->getTpcSubEp3ResVal(cent9);
		  mTpcEpManager->fillTpcSubEpV3(pt, v3Tpc, reweight);
		}
	      }
	    }
	  }
	}
      }
      mZdcEpManager->clearZdcEpManager();
      mEpdEpManager->clearEpdEpManager();
      mTpcEpManager->clearTpcEpManager();
    }
  }

  return kStOK;
}
