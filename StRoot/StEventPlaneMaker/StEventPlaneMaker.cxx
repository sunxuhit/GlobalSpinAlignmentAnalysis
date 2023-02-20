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
    str_mOutPutGainCorr = Form("./file_%s_GainCorr_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 1) // fill Re-Center Correction for ZDC & EPD & TPC Sub EP
  {
    str_mOutPutReCenterPar = Form("./file_%s_ReCenterPar_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 2) // fill Shift Correction for ZDC & EPD & TPC Sub EP
  {
    str_mOutPutShiftPar = Form("./file_%s_ShiftPar_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 3) // fill Shift Correction for ZDC & EPD Full EP
  {
    str_mOutPutShiftPar = Form("./file_%s_ShiftParFull_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 4) // fill Event Plane Resolution for ZDC & EPD & TPC Sub EP
  {
    str_mOutPutResolution = Form("./file_%s_EpResolution_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 5) // fill Charged Hadron v1 for ZDC and EPD & v2 and v3 for TPC
  {
    str_mOutPutFlow = Form("./file_%s_ChargedFlow_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
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
  }

  if(mMode == 0)
  { // fill Gain Correction Factors for ZDC and phi Weight for EPD
    file_mOutPutGainCorr = new TFile(str_mOutPutGainCorr.c_str(),"RECREATE");
    mZdcEpManager->initZdcGain(); // ZDC

    mEpdEpManager->initEpdPhiWgt(); // EPD
    mEpdEpManager->initEpdSubEpRaw();
  }
  if(mMode == 1)
  { // fill ReCenter Correction Parameters for ZDC-SMD & EPD & TPC Sub EP
    file_mOutPutReCenterPar = new TFile(str_mOutPutReCenterPar.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->initZdcReCenter();
    mZdcEpManager->initZdcSubEpRaw();

    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->initEpdReCtr();
    mEpdEpManager->initEpdSubEpWgt();

    mTpcEpManager->initTpcReCenter(); // TPC
    mTpcEpManager->initTpcSubEpRaw();
  }
  if(mMode == 2)
  { // fill Shift Correction Parameters for ZDC-SMD & EPD & TPC Sub EP
    file_mOutPutShiftPar = new TFile(str_mOutPutShiftPar.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->readZdcReCenter();
    mZdcEpManager->initZdcShift();
    mZdcEpManager->initZdcSubEpReCenter();

    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->readEpdReCtr();
    mEpdEpManager->initEpdShift();
    mEpdEpManager->initEpdSubEpReCtr();

    mTpcEpManager->readTpcReCenter(); // TPC
    mTpcEpManager->initTpcShift();
    mTpcEpManager->initTpcSubEpReCenter();
  }
  if(mMode == 3)
  { // fill Shift Correction Parameters for ZDC-SMD & EPD Full EP
    file_mOutPutShiftPar = new TFile(str_mOutPutShiftPar.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->readZdcReCenter();
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
  { // fill Event Plane Resolution for ZDC-SMD & EPD & TPC Sub EP
    file_mOutPutResolution = new TFile(str_mOutPutResolution.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->readZdcReCenter();
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

    mTpcEpManager->readTpcReCenter(); // TPC
    mTpcEpManager->readTpcShift();
    mTpcEpManager->initTpcResolution();
    mTpcEpManager->initTpcSubEpShift();
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
      mZdcEpManager->writeZdcReCenter(); // ZDC
      mZdcEpManager->writeZdcSubEpRaw();

      mEpdEpManager->writeEpdReCtr(); // EPD
      mEpdEpManager->writeEpdSubEpWgt();

      mTpcEpManager->writeTpcReCenter(); // TPC
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
      mZdcEpManager->writeZdcSubEpReCenter();

      mEpdEpManager->writeEpdShift(); // EPD
      mEpdEpManager->writeEpdSubEpReCtr();

      mTpcEpManager->writeTpcShift(); // TPC
      mTpcEpManager->writeTpcSubEpReCenter();
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
    const int grefMult     = mPicoEvent->grefMult();
    const TVector3 primVtx = mPicoEvent->primaryVertex();
    const double vx        = mPicoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
    const double vy        = mPicoEvent->primaryVertex().y();
    const double vz        = mPicoEvent->primaryVertex().z();
    const double vzVpd     = mPicoEvent->vzVpd();
    const double zdcX      = mPicoEvent->ZDCx();
    // const unsigned short nBTofHits  = mPicoEvent->btofTrayMultiplicity();
    const unsigned int nBTofHits    = mPicoDst->numberOfBTofHits(); // get number of tof hits
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
    const int triggerBin  = mAnaUtils->getTriggerBin(mPicoEvent);
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

    if(!isPileUpEvent && mAnaCut->isGoodCentrality(cent9) && mAnaCut->passEventCut(mPicoEvent))
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
	      // cout << "iEastWest = " << iEastWest << ", iVertHori = " << iVertHori << ", iSlat = " << iSlat << ", zdc = " << mZdcSmdCorrection->getZdcSmd(iEastWest,iVertHori,iSlat) << endl;
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
	    mEpdEpManager->fillEpdPhiWgtEast(picoEpdHit, primVtx);
	  }
	  if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	  {
	    mEpdEpManager->addHitRawWest(picoEpdHit, primVtx);
	    mEpdEpManager->fillEpdPhiWgtWest(picoEpdHit, primVtx);
	  }
	}

	TVector2 Q1EpdEast = mEpdEpManager->getQ1VecRawEast(); // get Q1Vector from EPD
	TVector2 Q1EpdWest = mEpdEpManager->getQ1VecRawWest();
	TVector2 Q1EpdFull = mEpdEpManager->getQ1VecRawFull();
	if( !(Q1EpdEast.Mod() < 1e-10 || Q1EpdWest.Mod() < 1e-10 || Q1EpdFull.Mod() < 1e-10) ) // EPD EP
	{
	  const double Psi1RawEast = mEpdEpManager->getPsi1RawEast();
	  const double Psi1RawWest = mEpdEpManager->getPsi1RawWest();
	  const double Psi1RawFull = mEpdEpManager->getPsi1RawFull();
	  mEpdEpManager->fillEpdSubEpRaw(Psi1RawEast,Psi1RawWest,Psi1RawFull);
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
	      mEpdEpManager->fillEpdReCtrEast(picoEpdHit, primVtx); // fill EPD ReCenter Parameters East
	    }
	    if( mAnaCut->passHitEpdEpWest(picoEpdHit) ) // positive eta
	    {
	      mEpdEpManager->addHitWgtWest(picoEpdHit, primVtx);
	      mEpdEpManager->fillEpdReCtrWest(picoEpdHit, primVtx); // fill EPD ReCenter Parameters West
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

	  if(mMode == 1) // calculate raw Q2Vector and Q3Vector from TPC
	  {
	    if( mAnaCut->passTrackTpcEpEast(picoTrack, primVtx) ) // negative eta
	    {
	      mTpcEpManager->addTrackRawEast(picoTrack);
	      mTpcEpManager->fillTpcReCenterEast(picoTrack); // fill TPC ReCenter Parameters East
	    }
	    if( mAnaCut->passTrackTpcEpWest(picoTrack, primVtx) ) // positive eta
	    {
	      mTpcEpManager->addTrackRawWest(picoTrack);
	      mTpcEpManager->fillTpcReCenterWest(picoTrack); // fill TPC ReCenter Parameters West
	    }
	  }
	  if(mMode == 2 || mMode == 4) // calculate recentered Q2Vector and Q3Vector from TPC
	  {
	    if( mAnaCut->passTrackTpcEpEast(picoTrack, primVtx) ) // negative eta
	    {
	      mTpcEpManager->addTrackReCenterEast(picoTrack);
	    }
	    if( mAnaCut->passTrackTpcEpWest(picoTrack, primVtx) ) // positive eta
	    {
	      mTpcEpManager->addTrackReCenterWest(picoTrack);
	    }
	  }
	}

	if(mMode == 1) // fill recenter correction parameter for ZDC & EPD & TPC Sub EP
	{
	  TVector2 Q1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 Q1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 Q1ZdcFull = mZdcEpManager->getQ1VecFull(Q1ZdcEast,Q1ZdcWest,mMode); // TVector2 Q1ZdcFull = Q1ZdcWest-Q1ZdcEast;
	  if( !(Q1ZdcEast.Mod() < 1e-10 || Q1ZdcWest.Mod() < 1e-10 || Q1ZdcFull.Mod() < 1e-10) ) // ZDC EP
	  {
	    mZdcEpManager->fillZdcSubEpRaw(Q1ZdcEast,Q1ZdcWest,Q1ZdcFull);
	    mZdcEpManager->fillZdcReCenterEast(Q1ZdcEast);
	    mZdcEpManager->fillZdcReCenterWest(Q1ZdcWest);
	  }

	  TVector2 Q1EpdEast = mEpdEpManager->getQ1VecWgtEast(); // get Q1Vector from EPD
	  TVector2 Q1EpdWest = mEpdEpManager->getQ1VecWgtWest();
	  TVector2 Q1EpdFull = mEpdEpManager->getQ1VecWgtFull();
	  if( !(Q1EpdEast.Mod() < 1e-10 || Q1EpdWest.Mod() < 1e-10 || Q1EpdFull.Mod() < 1e-10) ) // EPD EP
	  {
	    const double Psi1WgtEast = mEpdEpManager->getPsi1WgtEast(); // EPD EP
	    const double Psi1WgtWest = mEpdEpManager->getPsi1WgtWest();
	    const double Psi1WgtFull = mEpdEpManager->getPsi1WgtFull();
	    mEpdEpManager->fillEpdSubEpWgt(Psi1WgtEast, Psi1WgtWest, Psi1WgtFull);
	  }

	  const int numTrackRawEast = mTpcEpManager->getNumTrkRawEast(); // TPC EP
	  const int numTrackRawWest = mTpcEpManager->getNumTrkRawWest();
	  if(mAnaCut->passNumTrackTpcSubEpRaw(numTrackRawEast, numTrackRawWest))
	  {
	    const double Psi2RawEast = mTpcEpManager->getPsi2RawEast();
	    const double Psi2RawWest = mTpcEpManager->getPsi2RawWest();
	    const double Psi3RawEast = mTpcEpManager->getPsi3RawEast();
	    const double Psi3RawWest = mTpcEpManager->getPsi3RawWest();
	    mTpcEpManager->fillTpcSubEpRaw(Psi2RawEast, Psi2RawWest, Psi3RawEast, Psi3RawWest);
	  }
	}
	if(mMode == 2) // fill shift correction parameter for ZDC & EPD & TPC Sub EP
	{
	  TVector2 Q1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 Q1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 Q1ZdcFull = mZdcEpManager->getQ1VecFull(Q1ZdcEast,Q1ZdcWest,mMode); // TVector2 Q1ZdcFull = Q1ZdcWest-Q1ZdcEast;
	  if( !(Q1ZdcEast.Mod() < 1e-10 || Q1ZdcWest.Mod() < 1e-10 || Q1ZdcFull.Mod() < 1e-10) ) // ZDC EP
	  {
	    mZdcEpManager->fillZdcSubEpReCenter(Q1ZdcEast,Q1ZdcWest,Q1ZdcFull);
	    mZdcEpManager->fillZdcShiftEast(Q1ZdcEast);
	    mZdcEpManager->fillZdcShiftWest(Q1ZdcWest);
	  }

	  TVector2 Q1EpdEast = mEpdEpManager->getQ1VecReCtrEast(); // get Q1Vector from EPD
	  TVector2 Q1EpdWest = mEpdEpManager->getQ1VecReCtrWest();
	  TVector2 Q1EpdFull = mEpdEpManager->getQ1VecReCtrFull();
	  if( !(Q1EpdEast.Mod() < 1e-10 || Q1EpdWest.Mod() < 1e-10 || Q1EpdFull.Mod() < 1e-10) ) // EPD EP
	  {
	    const double Psi1ReCtrEast = mEpdEpManager->getPsi1ReCtrEast(); // EPD EP
	    const double Psi1ReCtrWest = mEpdEpManager->getPsi1ReCtrWest();
	    const double Psi1ReCtrFull = mEpdEpManager->getPsi1ReCtrFull();
	    mEpdEpManager->fillEpdSubEpReCtr(Psi1ReCtrEast, Psi1ReCtrWest, Psi1ReCtrFull);
	    mEpdEpManager->fillEpdShiftEast();
	    mEpdEpManager->fillEpdShiftWest();
	  }

	  const int numTrackReCenterEast = mTpcEpManager->getNumTrkReCenterEast(); // TPC EP
	  const int numTrackReCenterWest = mTpcEpManager->getNumTrkReCenterWest();
	  if(mAnaCut->passNumTrackTpcSubEpReCenter(numTrackReCenterEast, numTrackReCenterWest))
	  {
	    const double Psi2ReCenterEast = mTpcEpManager->getPsi2ReCenterEast();
	    const double Psi2ReCenterWest = mTpcEpManager->getPsi2ReCenterWest();
	    const double Psi3ReCenterEast = mTpcEpManager->getPsi3ReCenterEast();
	    const double Psi3ReCenterWest = mTpcEpManager->getPsi3ReCenterWest();
	    mTpcEpManager->fillTpcSubEpReCenter(Psi2ReCenterEast, Psi2ReCenterWest, Psi3ReCenterEast, Psi3ReCenterWest);
	    mTpcEpManager->fillTpcShiftEast();
	    mTpcEpManager->fillTpcShiftWest();
	  }
	}
	if(mMode == 3) // fill shift correction parameter for ZDC & EPD Full EP
	{
	  TVector2 Q1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 Q1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 Q1ZdcFull = mZdcEpManager->getQ1VecFull(Q1ZdcEast,Q1ZdcWest,mMode); // TVector2 Q1ZdcFull = Q1ZdcWest-Q1ZdcEast;
	  if( !(Q1ZdcEast.Mod() < 1e-10 || Q1ZdcWest.Mod() < 1e-10 || Q1ZdcFull.Mod() < 1e-10) ) // ZDC EP
	  {
	    mZdcEpManager->fillZdcSubEpShift(Q1ZdcEast,Q1ZdcWest,Q1ZdcFull);
	    mZdcEpManager->fillZdcShiftFull(Q1ZdcFull);
	  }

	  TVector2 Q1EpdEast = mEpdEpManager->getQ1VecShiftEast(); // get Q1Vector from EPD
	  TVector2 Q1EpdWest = mEpdEpManager->getQ1VecShiftWest();
	  TVector2 Q1EpdFull = mEpdEpManager->getQ1VecShiftFull();
	  if( !(Q1EpdEast.Mod() < 1e-10 || Q1EpdWest.Mod() < 1e-10 || Q1EpdFull.Mod() < 1e-10) ) // EPD EP
	  {
	    const double Psi1ShiftEast = mEpdEpManager->getPsi1ShiftEast(); // EPD EP
	    const double Psi1ShiftWest = mEpdEpManager->getPsi1ShiftWest();
	    const double Psi1ShiftFull = mEpdEpManager->getPsi1ShiftFull();
	    mEpdEpManager->fillEpdSubEpShift(Psi1ShiftEast, Psi1ShiftWest, Psi1ShiftFull);
	    mEpdEpManager->fillEpdShiftFull();
	  }
	}
	if(mMode == 4) // fill event plane resolution for ZDC-SMD & EPD & TPC Sub EP
	{
	  TVector2 Q1ZdcEast = mZdcEpManager->getQ1VecEast(mMode); // get Q1Vector from ZDC
	  TVector2 Q1ZdcWest = mZdcEpManager->getQ1VecWest(mMode);
	  TVector2 Q1ZdcFull = mZdcEpManager->getQ1VecFull(Q1ZdcEast,Q1ZdcWest,mMode); // TVector2 Q1ZdcFull = Q1ZdcWest-Q1ZdcEast;
	  if( !(Q1ZdcEast.Mod() < 1e-10 || Q1ZdcWest.Mod() < 1e-10 || Q1ZdcFull.Mod() < 1e-10) ) // ZDC EP
	  {
	    mZdcEpManager->fillZdcSubEpShift(Q1ZdcEast,Q1ZdcWest,Q1ZdcFull);
	    mZdcEpManager->fillZdcFullEpShift(Q1ZdcFull);
	    mZdcEpManager->fillZdcResolution(Q1ZdcEast,Q1ZdcWest);
	  }

	  TVector2 Q1EpdEast     = mEpdEpManager->getQ1VecShiftEast(); // get Q1Vector from EPD
	  TVector2 Q1EpdWest     = mEpdEpManager->getQ1VecShiftWest();
	  TVector2 Q1EpdFull     = mEpdEpManager->getQ1VecShiftFull();
	  TVector2 Q1EpdFullCorr = mEpdEpManager->getQ1VecShiftFullCorr();
	  if( !(Q1EpdEast.Mod() < 1e-10 || Q1EpdWest.Mod() < 1e-10 || Q1EpdFull.Mod() < 1e-10) ) // EPD EP
	  {
	    const double Psi1ShiftEast     = mEpdEpManager->getPsi1ShiftEast(); // EPD EP
	    const double Psi1ShiftWest     = mEpdEpManager->getPsi1ShiftWest();
	    const double Psi1ShiftFull     = mEpdEpManager->getPsi1ShiftFull();
	    const double Psi1ShiftFullCorr = mEpdEpManager->getPsi1ShiftFullCorr();
	    mEpdEpManager->fillEpdSubEpShift(Psi1ShiftEast, Psi1ShiftWest, Psi1ShiftFull);
	    mEpdEpManager->fillEpdFullEpShift(Psi1ShiftFullCorr);
	    mEpdEpManager->fillEpdResolution(Psi1ShiftEast, Psi1ShiftWest);
	  }

	  const int numTrackReCenterEast = mTpcEpManager->getNumTrkReCenterEast(); // TPC EP
	  const int numTrackReCenterWest = mTpcEpManager->getNumTrkReCenterWest();
	  if(mAnaCut->passNumTrackTpcSubEpReCenter(numTrackReCenterEast, numTrackReCenterWest))
	  {
	    const double Psi2ShiftEast = mTpcEpManager->getPsi2ShiftEast();
	    const double Psi2ShiftWest = mTpcEpManager->getPsi2ShiftWest();
	    const double Psi3ShiftEast = mTpcEpManager->getPsi3ShiftEast();
	    const double Psi3ShiftWest = mTpcEpManager->getPsi3ShiftWest();
	    mTpcEpManager->fillTpcSubEpShift(Psi2ShiftEast, Psi2ShiftWest, Psi3ShiftEast, Psi3ShiftWest);
	    mTpcEpManager->fillTpcResolution(Psi2ShiftEast, Psi2ShiftWest, Psi3ShiftEast, Psi3ShiftWest);
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

