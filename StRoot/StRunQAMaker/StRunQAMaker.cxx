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
// #include "StThreeVectorF.hh"
// #include "StMessMgr.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StRunQAMaker/StRunQAMaker.h"
#include "StRoot/StRunQAMaker/StRunQAHistoManager.h"
#include "StRoot/StRunQAMaker/StRunQAProManager.h"

ClassImp(StRunQAMaker)

//-----------------------------------------------------------------------------
StRunQAMaker::StRunQAMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int beamType) : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;

  mType = beamType;
  str_mOutPutRunQA = Form("./file_%s_RunQA_%s.root",globCons::mBeamType[mType].c_str(),jobId.c_str());
}

//----------------------------------------------------------------------------- 
StRunQAMaker::~StRunQAMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StRunQAMaker::Init() 
{
  mRunQAHistoManager = new StRunQAHistoManager();
  mRunQAProManager   = new StRunQAProManager();
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
    const int runId    = mPicoEvent->runId();
    const int refMult  = mPicoEvent->refMult();
    const int grefMult = mPicoEvent->grefMult();
    const float vx     = mPicoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
    const float vy     = mPicoEvent->primaryVertex().y();
    const float vz     = mPicoEvent->primaryVertex().z();
    const float vzVpd  = mPicoEvent->vzVpd();
    const float zdcX   = mPicoEvent->ZDCx();
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
    const int cent9 = mRefMultCorr->getCentralityBin9(); // get Centrality9
    const double reweight = mRefMultCorr->getWeight(); // get Centrality reweight

    // vz sign
    int vzSign = 0; // 0 for -vz || 1 for +vz
    vz > 0.0 ? vzSign = 1 : vzSign = 0;

    const int runIndex = mAnaUtils->findRunIndex(runId); // find run index for a specific run
    const int triggerBin = mAnaUtils->getTriggerBin(mPicoEvent);
    // cout << "runId = " << runId << ", runIndex = " << runIndex << endl;
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StAnalysisUtils! Skip!" << endm;
      return kStErr;
    }

    bool isPileUpEventStAnalysisCut = mAnaCut->isPileUpEvent(1.0*refMult, 1.0*nBTofMatch,vz); // alway return false for Isobar
    bool isPileUpEventStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut(1.0*refMult, 1.0*nBTofMatch,vz); // valid for Isobar
    bool isPileUpEvent = isPileUpEventStAnalysisCut || isPileUpEventStRefMultCorr;
    // cout << "isPileUpEvent = " << isPileUpEvent << ", isPileUpEventStAnalysisCut = " << isPileUpEventStAnalysisCut << ", isPileUpEventStRefMultCorr = " << isPileUpEventStRefMultCorr << endl;

    // fill QA before event cuts
    mRunQAProManager->fillRunQA_Event(triggerBin,runIndex,refMult,grefMult,zdcX,vx,vy,vz,0);
    mRunQAHistoManager->fillEventQA_RefMult(triggerBin,refMult,grefMult,cent9,reweight,nBTofHits,nBTofMatch,0); // wo event cut
    mRunQAHistoManager->fillEventQA_Vertex(triggerBin,vx,vy,vz,vzVpd,0);
    mRunQAHistoManager->fillEventQA_Trigger(triggerBin,0);

    if(!isPileUpEvent && mAnaCut->passEventCut(mPicoDst))
    { // apply Event Cuts for anlaysis 
      mRunQAProManager->fillRunQA_Event(triggerBin,runIndex,refMult,grefMult,zdcX,vx,vy,vz,1);
      mRunQAHistoManager->fillEventQA_RefMult(triggerBin,refMult,grefMult,cent9,reweight,nBTofHits,nBTofMatch,1); // with event cut
      mRunQAHistoManager->fillEventQA_Vertex(triggerBin,vx,vy,vz,vzVpd,1);
      mRunQAHistoManager->fillEventQA_Trigger(triggerBin,1);

      for(unsigned int i_track = 0; i_track < nTracks; i_track++) // track loop
      {
	StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(i_track);
	if(!picoTrack)
	{
	  continue;
	}

	// get pico track info
	// TVector3 primMom;
	// float primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	// float primPy    = picoTrack->pMom().y();
	// float primPz    = picoTrack->pMom().z();
	// primMom.SetXYZ(primPx,primPy,primPz);

	// TVector3 globMom;
	// float globPx     = picoTrack->gMom().x(); // x works for both TVector3 and StThreeVectorF
	// float globPy     = picoTrack->gMom().y();
	// float globPz     = picoTrack->gMom().z();
	// globMom.SetXYZ(globPx,globPy,globPz);

	const TVector3 primMom = picoTrack->pMom();
	const TVector3 globMom = picoTrack->gMom();

	const float gDCA    = picoTrack->gDCA(vx,vy,vz);
	const int nHitsFit  = picoTrack->nHitsFit();
	const int nHitsMax  = picoTrack->nHitsMax();
	const int nHitsDEdx = picoTrack->nHitsDedx();
	const float dEdx    = picoTrack->dEdx();
	const short charge  = picoTrack->charge();
	const float beta    = mAnaUtils->getBeta(mPicoDst,i_track);
	const float mass2   = mAnaUtils->getPrimaryMass2(mPicoDst,i_track);

	mRunQAHistoManager->fillTrackQA_Kinematics(triggerBin,primMom,globMom, 0); // wo track cut
	mRunQAHistoManager->fillTrackQA_Quliaty(triggerBin,gDCA,nHitsFit,nHitsMax,nHitsDEdx,0);
	mRunQAHistoManager->fillTrackQA_PID(triggerBin,primMom.Mag(),charge,dEdx,beta,mass2,0);
	mRunQAProManager->fillRunQA_Track(triggerBin,runIndex,gDCA,nHitsFit,primMom,globMom,0);
	if( mAnaCut->passTrackQA(picoTrack,mPicoEvent) ) // apply QA track cut
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

