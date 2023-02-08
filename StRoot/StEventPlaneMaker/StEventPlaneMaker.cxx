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
// #include "StThreeVectorF.hh"
// #include "StMessMgr.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneMaker.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"

ClassImp(StEventPlaneMaker)

//-----------------------------------------------------------------------------
StEventPlaneMaker::StEventPlaneMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int mode, const int beamType) : StMaker(name), mMode(mode), mType(beamType)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;
  // mMode = Mode;
  // mType = beamType;

  if(mMode == 0) // Gain Correction for ZDC-SMD
  {
    str_mOutPutGainCorr = Form("./file_%s_GainCorr_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
  if(mMode == 1) // Re-Center Correction for ZDC-SMD & TPC
  {
    str_mOutPutReCenterPar = Form("./file_%s_ReCenterParameter_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
  }
}

//----------------------------------------------------------------------------- 
StEventPlaneMaker::~StEventPlaneMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Init() 
{
  mZdcEpManager = new StZdcEpManager(mType); // initialize ZDC EP Manager
  mAnaCut       = new StAnalysisCut(mType);
  mAnaUtils     = new StAnalysisUtils(mType);
  mAnaUtils->initRunIndex(); // initialize std::map for run index

  if(!mRefMultCorr)
  {
    if( mAnaCut->isIsobar() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr_Isobar();
  }

  if(mMode == 0)
  { // fill Gain Correction Factors for ZDC-SMD
    file_mOutPutGainCorr = new TFile(str_mOutPutGainCorr.c_str(),"RECREATE");
    mZdcEpManager->initZdcGainCorr();
  }
  if(mMode == 1)
  { // fill ReCenter Correction Parameters for ZDC-SMD & TPC
    file_mOutPutReCenterPar = new TFile(str_mOutPutReCenterPar.c_str(),"RECREATE");
    mZdcEpManager->initZdcReCenter();
    mZdcEpManager->initZdcRawEP();
    mZdcEpManager->readGainCorr();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Finish() 
{
  if(mMode == 0)
  {
    if(str_mOutPutGainCorr != "")
    {
      file_mOutPutGainCorr->cd();
      mZdcEpManager->writeZdcGainCorr();
      file_mOutPutGainCorr->Close();
    }
  }
  if(mMode == 1)
  {
    if(str_mOutPutReCenterPar != "")
    {
      file_mOutPutReCenterPar->cd();
      mZdcEpManager->writeZdcReCenter();
      mZdcEpManager->writeZdcRawEP();
      file_mOutPutGainCorr->Close();
    }
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
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
    // if(mRefMultCorr->isBadRun(runId))
    // {
    //   LOG_ERROR << "Bad Run from StRefMultCorr! Skip!" << endm;
    //   return kStErr;
    // }

    const int cent9       = mRefMultCorr->getCentralityBin9(); // get Centrality9
    const double reweight = mRefMultCorr->getWeight(); // get weight
    const int runIndex    = mAnaUtils->findRunIndex(runId); // find run index for a specific run
    const int triggerBin  = mAnaUtils->getTriggerBin(mPicoEvent);
    const int vzBin       = mAnaUtils->getVzBin(vz); // 0 for -vz || 1 for +vz
    // cout << "runId = " << runId << ", runIndex = " << runIndex << endl;
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StEventPlaneUtility! Skip!" << endm;
      return kStErr;
    }

    bool isPileUpEventStAnalysisCut = mAnaCut->isPileUpEvent(1.0*refMult, 1.0*nBTofMatch,vz); // alway return false for Isobar
    bool isPileUpEventStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut(1.0*refMult, 1.0*nBTofMatch,vz); // valid for Isobar
    bool isPileUpEvent = isPileUpEventStAnalysisCut || isPileUpEventStRefMultCorr;
    // cout << "isPileUpEvent = " << isPileUpEvent << ", isPileUpEventStAnalysisCut = " << isPileUpEventStAnalysisCut << ", isPileUpEventStRefMultCorr = " << isPileUpEventStRefMultCorr << endl;

    if(!isPileUpEvent && mAnaCut->passEventCut(mPicoDst))
    { // apply Event Cuts for anlaysis 
      // ZDC-SMD EP
      mZdcEpManager->initZdcEp(cent9,runIndex,vzBin);
      if(mMode == 0)
      { // fill Gain Correction Factors for BBC & ZDC
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
	      mZdcEpManager->fillZdcGainCorr(iEastWest,iVertHori,iSlat,mZdcEpManager->getZdcSmd(iEastWest,iVertHori,iSlat));
	      // cout << "iEastWest = " << iEastWest << ", iVertHori = " << iVertHori << ", iSlat = " << iSlat << ", zdc = " << mZdcSmdCorrection->getZdcSmd(iEastWest,iVertHori,iSlat) << endl;
	    }
	  }
	}
      }
      if(mMode > 0)
      {
	for(int iSlat = 0; iSlat < 8; ++iSlat) // read in raw ADC value from ZDC-SMD
	{
	  mZdcEpManager->setZdcSmdGainCorr(0,0,iSlat,mPicoEvent->ZdcSmdEastVertical(iSlat));
	  mZdcEpManager->setZdcSmdGainCorr(0,1,iSlat,mPicoEvent->ZdcSmdEastHorizontal(iSlat));
	  mZdcEpManager->setZdcSmdGainCorr(1,0,iSlat,mPicoEvent->ZdcSmdWestVertical(iSlat));
	  mZdcEpManager->setZdcSmdGainCorr(1,1,iSlat,mPicoEvent->ZdcSmdWestHorizontal(iSlat));
	}

	if(mMode == 1) // apply gain correction and fill recenter correction parameter
	{
	  TVector2 QEast = mZdcEpManager->getQEast(mMode);
	  TVector2 QWest = mZdcEpManager->getQWest(mMode);
	  TVector2 QFull = QWest-QEast;
	  if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
	  {
	    mZdcEpManager->fillZdcReCenterEast(QEast);
	    mZdcEpManager->fillZdcReCenterWest(QWest);
	    mZdcEpManager->fillZdcRawEP(QEast,QWest,QFull);
	  }
	}
      }
      mZdcEpManager->clearZdcEp();
    }
  }

  return kStOK;
}

