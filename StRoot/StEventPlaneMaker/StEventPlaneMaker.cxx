#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneMaker.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneUtility.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCut.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneHistoManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneProManager.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"

#include <algorithm>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

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
  mEventPlaneCut = new StEventPlaneCut(mType);
  mEventPlaneHistoManager = new StEventPlaneHistoManager(mType);
  mEventPlaneUtility = new StEventPlaneUtility(mType);
  mEventPlaneUtility->initRunIndex(); // initialize std::map for run index
  mEventPlaneProManager = new StEventPlaneProManager(mType);
  mZdcEpManager = new StZdcEpManager(mType); // initialize ZDC EP Manager

  if(!mRefMultCorr)
  {
    // if(!mEventPlaneCut->isBES()) mRefMultCorr = CentralityMaker::instance()->getgRefMultCorr_Run14_AuAu200_VpdMB5_P16id(); // 200GeV_2014
    // if(mEventPlaneCut->isBES()) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr(); // BESII
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr(); // BESII
  }

  if(mMode == 0)
  { // fill Gain Correction Factors for ZDC-SMD
    file_mOutPutGainCorr = new TFile(str_mOutPutGainCorr.c_str(),"RECREATE");
    mEventPlaneHistoManager->initZdcGainCorr();
  }
  if(mMode == 1)
  { // fill ReCenter Correction Parameters for ZDC-SMD & TPC
    file_mOutPutReCenterPar = new TFile(str_mOutPutReCenterPar.c_str(),"RECREATE");
    mEventPlaneProManager->initZdcReCenter();
    mEventPlaneHistoManager->initZdcRawEP();
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
      mEventPlaneHistoManager->writeZdcGainCorr();
      file_mOutPutGainCorr->Close();
    }
  }
  if(mMode == 1)
  {
    if(str_mOutPutReCenterPar != "")
    {
      file_mOutPutReCenterPar->cd();
      mEventPlaneProManager->writeZdcReCenter();
      mEventPlaneHistoManager->writeZdcRawEP();
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
  if( mEventPlaneCut->isMinBias(mPicoEvent) )
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
    const unsigned int numOfBTofHits = mPicoDst->numberOfBTofHits(); // get number of tof hits
    // const unsigned short numOfBTofHits = mPicoEvent->btofTrayMultiplicity();
    const unsigned short numOfBTofMatch = mPicoEvent->nBTOFMatch(); // get number of tof match points
    const unsigned int nTracks = mPicoDst->numberOfTracks(); // get number of tracks

    // StRefMultCorr Cut & centrality
    if(!mRefMultCorr)
    {
      LOG_WARN << " No mRefMultCorr! Skip! " << endl;
      return kStErr;
    }

    mRefMultCorr->init(runId);
    // if(!mEventPlaneCut->isBES()) mRefMultCorr->initEvent(grefMult,vz,zdcX); // 200GeV_2014
    // if(mEventPlaneCut->isBES()) mRefMultCorr->initEvent(refMult,vz,zdcX); // BES-II might need Luminosity corrections
    mRefMultCorr->initEvent(refMult,vz,zdcX); // BES-II might need Luminosity corrections

    // if(mRefMultCorr->isBadRun(runId))
    // {
    //   LOG_ERROR << "Bad Run from StRefMultCorr! Skip!" << endm;
    //   return kStErr;
    // }

    // vz sign
    int vz_sign = 0; // 0 for -vz || 1 for vz
    vz > 0.0 ? vz_sign = 1 : vz_sign = 0;

    const int cent9 = mRefMultCorr->getCentralityBin9(); // get Centrality9
    const double reweight = mRefMultCorr->getWeight(); // get weight
    const int runIndex = mEventPlaneUtility->findRunIndex(runId); // find run index for a specific run
    const int triggerBin = mEventPlaneCut->getTriggerBin(mPicoEvent);
    // cout << "runId = " << runId << ", runIndex = " << runIndex << endl;
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StEventPlaneUtility! Skip!" << endm;
      return kStErr;
    }

    // bool isPileUpEventStEventPlaneCut = mEventPlaneCut->isPileUpEvent(grefMult,numOfBTofMatch,numOfBTofHits); // 200GeV
    // if(mEventPlaneCut->isBES()) isPileUpEventStEventPlaneCut = mEventPlaneCut->isPileUpEvent(refMult,numOfBTofMatch,numOfBTofHits); // 54 GeV | always return false for 27 GeV
    // bool isPileUpEventStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut(1.0*refMult, 1.0*numOfBTofMatch); // 27 GeV | always return !true for other energies
    // bool isPileUpEvent = isPileUpEventStEventPlaneCut || isPileUpEventStRefMultCorr;
    // cout << "isPileUpEvent = " << isPileUpEvent << ", isPileUpEventStEventPlaneCut = " << isPileUpEventStEventPlaneCut << ", isPileUpEventStRefMultCorr = " << isPileUpEventStRefMultCorr << endl;
    bool isPileUpEvent = false;

    if(mEventPlaneCut->passEventCut(mPicoDst) && !isPileUpEvent)
    { // apply Event Cuts for anlaysis 
      // ZDC-SMD EP
      mZdcEpManager->initZdcEp(cent9,runIndex,vz_sign);
      if(mMode == 0)
      { // fill Gain Correction Factors for BBC & ZDC
	for(int i_slat = 0; i_slat < 8; ++i_slat) // read in raw ADC value from ZDC-SMD
	{
	  mZdcEpManager->setZdcSmd(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	  mZdcEpManager->setZdcSmd(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	  mZdcEpManager->setZdcSmd(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	  mZdcEpManager->setZdcSmd(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
	}
	for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest) // fill ZDC Gain Correction Histograms
	{
	  for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
	  {
	    for(int i_slat = 0; i_slat < 8; ++i_slat)
	    {
	      mEventPlaneHistoManager->fillZdcGainCorr(i_eastwest,i_verthori,i_slat,runIndex,mZdcEpManager->getZdcSmd(i_eastwest,i_verthori,i_slat));
	      // cout << "i_eastwest = " << i_eastwest << ", i_verthori = " << i_verthori << ", i_slat = " << i_slat << ", zdc = " << mZdcSmdCorrection->getZdcSmd(i_eastwest,i_verthori,i_slat) << endl;
	    }
	  }
	}
      }
      if(mMode > 0)
      {
	for(int i_slat = 0; i_slat < 8; ++i_slat) // read in raw ADC value from ZDC-SMD
	{
	  mZdcEpManager->setZdcSmdGainCorr(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	  mZdcEpManager->setZdcSmdGainCorr(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	  mZdcEpManager->setZdcSmdGainCorr(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	  mZdcEpManager->setZdcSmdGainCorr(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
	}

	if(mMode == 1) // apply gain correction and fill recenter correction parameter
	{
	  TVector2 QEast = mZdcEpManager->getQEast(mMode);
	  TVector2 QWest = mZdcEpManager->getQWest(mMode);
	  TVector2 QFull = QWest-QEast;
	  if( !(QEast.Mod() < 1e-10 || QWest.Mod() < 1e-10 || QFull.Mod() < 1e-10) )
	  {
	    mEventPlaneProManager->fillZdcReCenterEast(QEast,cent9,runIndex,vz_sign);
	    mEventPlaneProManager->fillZdcReCenterWest(QWest,cent9,runIndex,vz_sign);
	    mEventPlaneHistoManager->fillZdcRawEP(QEast,QWest,QFull,cent9,runIndex);
	  }
	}
      }
    }
  }

  return kStOK;
}

