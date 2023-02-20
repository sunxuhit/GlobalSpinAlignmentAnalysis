#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TVector3.h"
#include "TString.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StRunQAMaker/StRunQAMaker.h"
#include "StRoot/StRunQAMaker/StRunQAHistoManager.h"
#include "StRoot/StRunQAMaker/StRunQAProManager.h"

ClassImp(StRunQAMaker)

//-----------------------------------------------------------------------------
StRunQAMaker::StRunQAMaker(const char* name, StPicoDstMaker *picoMaker, string jobId, int beamType) : StMaker(name), mType(beamType)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;

  // mType = beamType;
  str_mOutPutRunQA = Form("./file_%s_RunQA_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
}

//----------------------------------------------------------------------------- 
StRunQAMaker::~StRunQAMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StRunQAMaker::Init() 
{
  mRunQAHistoManager = new StRunQAHistoManager(mType);
  mRunQAProManager   = new StRunQAProManager(mType);
  mAnaCut            = new StAnalysisCut(mType);
  mAnaUtils          = new StAnalysisUtils(mType);
  mAnaUtils->initRunIndex(); // initialize std::map for run index

  if(!mRefMultCorr)
  {
    if( mAnaCut->isIsobar() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr_Isobar();
  }

  file_mOutPutRunQA= new TFile(str_mOutPutRunQA.c_str(),"RECREATE");
  file_mOutPutRunQA->cd();
  mRunQAHistoManager->initEventQA();
  mRunQAHistoManager->initTrackQA();
  mRunQAProManager->initRunQA();

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StRunQAMaker::Finish() 
{
  if(str_mOutPutRunQA != "")
  {
    file_mOutPutRunQA->cd();
    mRunQAHistoManager->writeEventQA();
    mRunQAHistoManager->writeTrackQA();
    mRunQAProManager->writeRunQA();
    file_mOutPutRunQA->Close();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StRunQAMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
int StRunQAMaker::Make() 
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
    /*
    if(mRefMultCorr->isBadRun(runId))
    {
      LOG_ERROR << "Bad Run from StRefMultCorr! Skip!" << endm;
      return kStErr;
    }
    */

    const int cent9       = mRefMultCorr->getCentralityBin9(); // get Centrality9
    const double reweight = mRefMultCorr->getWeight(); // get Centrality reweight
    const int runIndex    = mAnaUtils->findRunIndex(runId); // find run index for a specific run
    const int triggerBin  = mAnaUtils->getTriggerBin(mPicoEvent);
    // const int vzBin       = mAnaUtils->getVzBin(vz); // 0 for -vz || 1 for +vz
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

    // fill QA before event cuts
    mRunQAProManager->fillRunQA_Event(triggerBin,runIndex,refMult,grefMult,zdcX,vx,vy,vz,0);
    mRunQAHistoManager->fillEventQA_RefMult(triggerBin,refMult,grefMult,cent9,reweight,nBTofHits,nBTofMatch,0); // wo event cut
    mRunQAHistoManager->fillEventQA_Vertex(triggerBin,vx,vy,vz,vzVpd,0);
    mRunQAHistoManager->fillEventQA_Trigger(triggerBin,0);

    if(!isPileUpEvent && mAnaCut->passEventCut(mPicoEvent))
    { // apply Event Cuts for anlaysis 
      mRunQAProManager->fillRunQA_Event(triggerBin,runIndex,refMult,grefMult,zdcX,vx,vy,vz,1);
      mRunQAHistoManager->fillEventQA_RefMult(triggerBin,refMult,grefMult,cent9,reweight,nBTofHits,nBTofMatch,1); // with event cut
      mRunQAHistoManager->fillEventQA_Vertex(triggerBin,vx,vy,vz,vzVpd,1);
      mRunQAHistoManager->fillEventQA_Trigger(triggerBin,1);

      for(unsigned int iTrack = 0; iTrack < nTracks; iTrack++) // track loop
      {
	StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	if(!picoTrack)
	{
	  continue;
	}

	const TVector3 primMom = picoTrack->pMom();
	const TVector3 globMom = picoTrack->gMom();

	const double gDCA    = picoTrack->gDCA(vx,vy,vz);
	const int nHitsFit   = picoTrack->nHitsFit();
	const int nHitsMax   = picoTrack->nHitsMax();
	const int nHitsDEdx  = picoTrack->nHitsDedx();
	const double dEdx    = picoTrack->dEdx();
	const short charge   = picoTrack->charge();
	const double beta    = mAnaUtils->getBeta(mPicoDst,iTrack);
	const double mass2   = mAnaUtils->getPrimaryMass2(mPicoDst,iTrack);

	mRunQAHistoManager->fillTrackQA_Kinematics(triggerBin,primMom,globMom, 0); // wo track cut
	mRunQAHistoManager->fillTrackQA_Quliaty(triggerBin,gDCA,nHitsFit,nHitsMax,nHitsDEdx,0);
	mRunQAHistoManager->fillTrackQA_PID(triggerBin,primMom.Mag(),charge,dEdx,beta,mass2,0);
	mRunQAProManager->fillRunQA_Track(triggerBin,runIndex,gDCA,nHitsFit,primMom,globMom,0);
	if( mAnaCut->passTrackQA(picoTrack,primVtx) ) // apply QA track cut
	{
	  mRunQAHistoManager->fillTrackQA_Kinematics(triggerBin,primMom,globMom, 1); // with track cut
	  mRunQAHistoManager->fillTrackQA_Quliaty(triggerBin,gDCA,nHitsFit,nHitsMax,nHitsDEdx,1);
	  mRunQAHistoManager->fillTrackQA_PID(triggerBin,primMom.Mag(),charge,dEdx,beta,mass2,1);
	  mRunQAProManager->fillRunQA_Track(triggerBin,runIndex,gDCA,nHitsFit,primMom,globMom,1);
	}
      }
    }
  }

  return kStOK;
}

