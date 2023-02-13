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
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
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

  if(mMode == 0) // Gain Correction for ZDC-SMD & EPD
  {
    str_mOutPutGainCorr = Form("./file_%s_GainCorr_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 1) // Re-Center Correction for ZDC-SMD & EPD & TPC Sub EP
  {
    str_mOutPutReCenterPar = Form("./file_%s_ReCenterPar_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 2) // Shift Correction for ZDC-SMD & EPD & TPC Sub EP
  {
    str_mOutPutShiftPar = Form("./file_%s_ShiftPar_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 3) // Shift Correction for ZDC-SMD & EPD Full EP
  {
    str_mOutPutShiftPar = Form("./file_%s_ShiftParFull_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 4) // Event Plane Resolution for ZDC-SMD & EPD & TPC Sub EP
  {
    str_mOutPutResolution = Form("./file_%s_EpResolution_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
}

StEventPlaneMaker::~StEventPlaneMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Init() 
{
  mZdcEpManager = new StZdcEpManager(mType); // initialize ZDC EP Manager
  mTpcEpManager = new StTpcEpManager(mType); // initialize TPC EP Manager
  mAnaCut       = new StAnalysisCut(mType);
  mAnaUtils     = new StAnalysisUtils(mType);
  mAnaUtils->initRunIndex(); // initialize std::map for run index

  if(!mRefMultCorr)
  {
    if( mAnaCut->isIsobar() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr_Isobar();
  }

  if(mMode == 0)
  { // fill Gain Correction Factors for ZDC-SMDa & EPD
    file_mOutPutGainCorr = new TFile(str_mOutPutGainCorr.c_str(),"RECREATE");
    mZdcEpManager->initZdcGain(); // ZDC
  }
  if(mMode == 1)
  { // fill ReCenter Correction Parameters for ZDC-SMD & EPD & TPC Sub EP
    file_mOutPutReCenterPar = new TFile(str_mOutPutReCenterPar.c_str(),"RECREATE");
    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->initZdcReCenter();
    mZdcEpManager->initZdcSubEpRaw();

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
    const int runId    = mPicoEvent->runId();
    const int refMult  = mPicoEvent->refMult();
    const int grefMult = mPicoEvent->grefMult();
    const double vx    = mPicoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
    const double vy    = mPicoEvent->primaryVertex().y();
    const double vz    = mPicoEvent->primaryVertex().z();
    const double vzVpd = mPicoEvent->vzVpd();
    const double zdcX  = mPicoEvent->ZDCx();
    // const unsigned short nBTofHits  = mPicoEvent->btofTrayMultiplicity();
    const unsigned int nBTofHits    = mPicoDst->numberOfBTofHits(); // get number of tof hits
    const unsigned short nBTofMatch = mPicoEvent->nBTOFMatch(); // get number of tof match points
    const unsigned int nTracks      = mPicoDst->numberOfTracks(); // get number of tracks

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

    if(!isPileUpEvent && mAnaCut->passEventCut(mPicoDst))
    { // apply Event Cuts for anlaysis 
      mZdcEpManager->initZdcEpManager(cent9,runIndex,vzBin); // initialize ZDC EP Manager
      mTpcEpManager->initTpcEpManager(cent9,runIndex,vzBin); // initialize TPC EP Manager
      if(mMode == 0)
      { // fill Gain Correction Factors for ZDC-SMD & EPD (TODO)
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
      }
      if(mMode > 0)
      {
	for(int iSlat = 0; iSlat < 8; ++iSlat) // read in raw ADC value from ZDC-SMD & EPD and apply gain correction
	{
	  mZdcEpManager->setZdcSmdGainCorr(0,0,iSlat,mPicoEvent->ZdcSmdEastVertical(iSlat));
	  mZdcEpManager->setZdcSmdGainCorr(0,1,iSlat,mPicoEvent->ZdcSmdEastHorizontal(iSlat));
	  mZdcEpManager->setZdcSmdGainCorr(1,0,iSlat,mPicoEvent->ZdcSmdWestVertical(iSlat));
	  mZdcEpManager->setZdcSmdGainCorr(1,1,iSlat,mPicoEvent->ZdcSmdWestHorizontal(iSlat));
	}
	TVector2 QEast = mZdcEpManager->getQEast(mMode); // get QVector from ZDC
	TVector2 QWest = mZdcEpManager->getQWest(mMode);
	TVector2 QFull = mZdcEpManager->getQFull(QEast,QWest,mMode); // TVector2 QFull = QWest-QEast;

	for(unsigned int iTrack = 0; iTrack < nTracks; ++iTrack) // get Q2Vector and Q3Vector from TPC
	{
	  StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	  if(!picoTrack) continue;

	  if(mMode == 1) // get raw Q2Vector and Q3Vector from TPC
	  {
	    if( mAnaCut->passTrackTpcEpEast(picoTrack, mPicoEvent) ) // negative eta
	    {
	      mTpcEpManager->addTrackRawEast(picoTrack);
	      mTpcEpManager->fillTpcReCenterEast(picoTrack); // fill TPC ReCenter Parameters East
	    }
	    if( mAnaCut->passTrackTpcEpWest(picoTrack, mPicoEvent) ) // positive eta
	    {
	      mTpcEpManager->addTrackRawWest(picoTrack);
	      mTpcEpManager->fillTpcReCenterWest(picoTrack); // fill TPC ReCenter Parameters West
	    }
	  }
	  if(mMode == 2 || mMode == 4) // get recentered Q2Vector and Q3Vector from TPC
	  {
	    if( mAnaCut->passTrackTpcEpEast(picoTrack, mPicoEvent) ) // negative eta
	    {
	      mTpcEpManager->addTrackReCenterEast(picoTrack);
	    }
	    if( mAnaCut->passTrackTpcEpWest(picoTrack, mPicoEvent) ) // positive eta
	    {
	      mTpcEpManager->addTrackReCenterWest(picoTrack);
	    }
	  }
	}

	if(mMode == 1) // fill recenter correction parameter for ZDC-SMD & EPD & TPC Sub EP
	{
	  if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
	  {
	    mZdcEpManager->fillZdcReCenterEast(QEast);
	    mZdcEpManager->fillZdcReCenterWest(QWest);
	    mZdcEpManager->fillZdcSubEpRaw(QEast,QWest,QFull);
	  }

	  const int numTrackRawEast = mTpcEpManager->getNumTrkRawEast();
	  const int numTrackRawWest = mTpcEpManager->getNumTrkRawWest();
	  const double Psi2RawEast = mTpcEpManager->getPsi2RawEast();
	  const double Psi2RawWest = mTpcEpManager->getPsi2RawWest();
	  const double Psi3RawEast = mTpcEpManager->getPsi3RawEast();
	  const double Psi3RawWest = mTpcEpManager->getPsi3RawWest();
	  if(mAnaCut->passNumTrackTpcSubEpRaw(numTrackRawEast, numTrackRawWest))
	  {
	    mTpcEpManager->fillTpcSubEpRaw(Psi2RawEast, Psi2RawWest, Psi3RawEast, Psi3RawWest);
	  }
	}
	if(mMode == 2) // fill shift correction parameter for ZDC-SMD & EPD & TPC Sub EP
	{
	  if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
	  {
	    mZdcEpManager->fillZdcShiftEast(QEast);
	    mZdcEpManager->fillZdcShiftWest(QWest);
	    mZdcEpManager->fillZdcSubEpReCenter(QEast,QWest,QFull);
	  }
	}
	if(mMode == 3) // fill shift correction parameter for ZDC-SMD & EPD Full EP
	{
	  if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
	  {
	    mZdcEpManager->fillZdcShiftFull(QFull);
	    mZdcEpManager->fillZdcSubEpShift(QEast,QWest,QFull);
	  }
	}
	if(mMode == 4) // fill event plane resolution for ZDC-SMD & EPD & TPC Sub EP
	{
	  if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
	  {
	    mZdcEpManager->fillZdcResolution(QEast,QWest);
	    mZdcEpManager->fillZdcSubEpShift(QEast,QWest,QFull);
	    mZdcEpManager->fillZdcFullEpShift(QFull);
	  }
	}
      }
      mZdcEpManager->clearZdcEpManager();
      mTpcEpManager->clearTpcEpManager();
    }
  }

  return kStOK;
}

