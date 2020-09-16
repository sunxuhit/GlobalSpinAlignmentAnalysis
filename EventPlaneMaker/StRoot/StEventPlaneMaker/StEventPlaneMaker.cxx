#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"

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
StEventPlaneMaker::StEventPlaneMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int Mode, const int energy) : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;
  mMode = Mode;
  mEnergy = energy;

  if(mMode == 0)
  {
    mOutPut_GainCorr = Form("./file_%s_GainCorr_%s.root",recoEP::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
}

//----------------------------------------------------------------------------- 
StEventPlaneMaker::~StEventPlaneMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Init() 
{
  mEventPlaneCut = new StEventPlaneCut(mEnergy);
  mEventPlaneHistoManager = new StEventPlaneHistoManager();
  mEventPlaneUtility = new StEventPlaneUtility(mEnergy);
  mEventPlaneUtility->initRunIndex(); // initialize std::map for run index
  mEventPlaneProManager = new StEventPlaneProManager();
  mZdcEpManager = new StZdcEpManager(mEnergy); // initialize ZDC EP Manager

  if(!mRefMultCorr)
  {
    if(!mEventPlaneCut->isBES(mEnergy)) mRefMultCorr = CentralityMaker::instance()->getgRefMultCorr_Run14_AuAu200_VpdMB5_P16id(); // 200GeV_2014
    if(mEventPlaneCut->isBES(mEnergy)) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr(); // BESII
  }

  if(mMode == 0)
  { // fill Gain Correction Factors for BBC & ZDC
    mFile_GainCorr= new TFile(mOutPut_GainCorr.c_str(),"RECREATE");
    mEventPlaneHistoManager->initZdcGainCorr();
    // mEventPlaneProManager->InitReCenter();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StEventPlaneMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_GainCorr != "")
    {
      mFile_GainCorr->cd();
      mEventPlaneHistoManager->writeZdcGainCorr();
      mFile_GainCorr->Close();
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
    const int runId = mPicoEvent->runId();
    const int refMult = mPicoEvent->refMult();
    const int grefMult = mPicoEvent->grefMult();
    const float vx    = mPicoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
    const float vy    = mPicoEvent->primaryVertex().y();
    const float vz    = mPicoEvent->primaryVertex().z();
    const float vzVpd = mPicoEvent->vzVpd();
    const float zdcX  = mPicoEvent->ZDCx();
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
    if(!mEventPlaneCut->isBES(mEnergy)) mRefMultCorr->initEvent(grefMult,vz,zdcX); // 200GeV_2014
    if(mEventPlaneCut->isBES(mEnergy)) mRefMultCorr->initEvent(refMult,vz,zdcX); // BES-II might need Luminosity corrections

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
    // cout << "runId = " << runId << ", runIndex = " << runIndex << endl;
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StEventPlaneUtility! Skip!" << endm;
      return kStErr;
    }

    if(mMode == 0)
    { // fill Gain Correction Factors for BBC & ZDC
      if(mEventPlaneCut->passEventCut(mPicoDst))
      { // apply Event Cuts for anlaysis 
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
    }
  }

  return kStOK;
}

