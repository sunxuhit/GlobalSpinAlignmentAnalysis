#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"

#include "StRoot/StRunQAMaker/StRunQAMaker.h"
#include "StRoot/StRunQAMaker/StRunQACut.h"
#include "StRoot/StRunQAMaker/StRunQAHistoManager.h"
#include "StRoot/StRunQAMaker/StRunQAUtility.h"
#include "StRoot/StRunQAMaker/StRunQAProManager.h"
#include "StRoot/StRunQAMaker/StRunQACons.h"

#include <algorithm>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

ClassImp(StRunQAMaker)

//-----------------------------------------------------------------------------
StRunQAMaker::StRunQAMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int Mode, const int energy) : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;
  mMode = Mode;
  mEnergy = energy;

  if(mMode == 0)
  {
    // mOutPut_RunQA = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/QA/file_%s_QA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),jobId.c_str());
    mOutPut_RunQA = Form("./file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
}

//----------------------------------------------------------------------------- 
StRunQAMaker::~StRunQAMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StRunQAMaker::Init() 
{
  mRunQACut = new StRunQACut(mEnergy);
  mRunQAHistoManager = new StRunQAHistoManager();
  mRunQAUtility = new StRunQAUtility(mEnergy);
  mRunQAUtility->initRunIndex(); // initialize std::map for run index
  mRunQAProManager = new StRunQAProManager();

  if(!mRefMultCorr)
  {
    if(!mRunQACut->isBES()) mRefMultCorr = CentralityMaker::instance()->getgRefMultCorr_Run14_AuAu200_VpdMB5_P16id(); // 200GeV_2014
    if(mRunQACut->isBES()) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr(); // BESII
  }

  if(mMode == 0)
  { // QA
    mFile_QA= new TFile(mOutPut_RunQA.c_str(),"RECREATE");
    mFile_QA->cd();
    mRunQAHistoManager->initEventQA();
    mRunQAHistoManager->initTrackQA();
    mRunQAProManager->initRunQA();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StRunQAMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_RunQA != "")
    {
      mFile_QA->cd();
      mRunQAHistoManager->writeEventQA();
      mRunQAHistoManager->writeTrackQA();
      mRunQAProManager->writeRunQA();
      mFile_QA->Close();
    }
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
  if( mRunQACut->isMinBias(mPicoEvent) )
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
    if(!mRunQACut->isBES()) mRefMultCorr->initEvent(grefMult,vz,zdcX); // 200GeV_2014
    if(mRunQACut->isBES()) mRefMultCorr->initEvent(refMult,vz,zdcX); // BES-II might need Luminosity corrections

    /*
    if(mRefMultCorr->isBadRun(runId))
    {
      LOG_ERROR << "Bad Run from StRefMultCorr! Skip!" << endm;
      return kStErr;
    }
    */

    // vz sign
    int vz_sign = 0; // 0 for -vz || 1 for vz
    vz > 0.0 ? vz_sign = 1 : vz_sign = 0;

    const int cent9 = mRefMultCorr->getCentralityBin9(); // get Centrality9
    const double reweight = mRefMultCorr->getWeight(); // get Centrality9
    const int runIndex = mRunQAUtility->findRunIndex(runId); // find run index for a specific run
    const int triggerBin = mRunQACut->getTriggerBin(mPicoEvent);
    // cout << "runId = " << runId << ", runIndex = " << runIndex << endl;
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StRunQAUtility! Skip!" << endm;
      return kStErr;
    }

    bool isPileUpEventStRunQACut = mRunQACut->isPileUpEvent(grefMult,numOfBTofMatch,numOfBTofHits); // 200GeV
    if(mRunQACut->isBES()) isPileUpEventStRunQACut = mRunQACut->isPileUpEvent(refMult,numOfBTofMatch,numOfBTofHits); // 54 GeV | always return false for 27 GeV
    bool isPileUpEventStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut(1.0*refMult, 1.0*numOfBTofMatch); // 27 GeV | always return !true for other energies
    bool isPileUpEvent = isPileUpEventStRunQACut || isPileUpEventStRefMultCorr;
    // cout << "isPileUpEvent = " << isPileUpEvent << ", isPileUpEventStRunQACut = " << isPileUpEventStRunQACut << ", isPileUpEventStRefMultCorr = " << isPileUpEventStRefMultCorr << endl;

    if(mMode == 0)
    { // fill QA before event cuts
      mRunQAProManager->fillRunQA_Event(triggerBin,runIndex,refMult,grefMult,zdcX,vx,vy,vz,0);
      mRunQAHistoManager->fillEventQA_RefMult(triggerBin,refMult,grefMult,cent9,reweight,numOfBTofHits,numOfBTofMatch,0); // wo event cut
      mRunQAHistoManager->fillEventQA_Vertex(triggerBin,vx,vy,vz,vzVpd,0);
      mRunQAHistoManager->fillEventQA_Trigger(triggerBin,0);

      if(mRunQACut->passEventCut(mPicoDst) && !isPileUpEvent)
      { // apply Event Cuts for anlaysis 
	mRunQAProManager->fillRunQA_Event(triggerBin,runIndex,refMult,grefMult,zdcX,vx,vy,vz,1);
	mRunQAHistoManager->fillEventQA_RefMult(triggerBin,refMult,grefMult,cent9,reweight,numOfBTofHits,numOfBTofMatch,1); // with event cut
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
	  // TVector3 primMom = picoTrack->pMom();
	  // TVector3 globMom = picoTrack->gMom();
	  TVector3 primMom;
	  float primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
	  float primPy    = picoTrack->pMom().y();
	  float primPz    = picoTrack->pMom().z();
	  primMom.SetXYZ(primPx,primPy,primPz);

	  TVector3 globMom;
	  float globPx     = picoTrack->gMom().x(); // x works for both TVector3 and StThreeVectorF
	  float globPy     = picoTrack->gMom().y();
	  float globPz     = picoTrack->gMom().z();
	  globMom.SetXYZ(globPx,globPy,globPz);

	  float gDCA       = picoTrack->gDCA(vx,vy,vz);
	  int nHitsFit     = picoTrack->nHitsFit();
	  int nHitsMax     = picoTrack->nHitsMax();
	  int nHitsDEdx    = picoTrack->nHitsDedx();
	  float dEdx       = picoTrack->dEdx();
	  short charge     = picoTrack->charge();
	  float beta       = mRunQACut->getBeta(mPicoDst,i_track);
	  float mass2      = mRunQACut->getPrimaryMass2(mPicoDst,i_track);

	  mRunQAHistoManager->fillTrackQA_Kinematics(triggerBin,primMom,globMom, 0); // wo track cut
	  mRunQAHistoManager->fillTrackQA_Quliaty(triggerBin,gDCA,nHitsFit,nHitsMax,nHitsDEdx,0);
	  mRunQAHistoManager->fillTrackQA_PID(triggerBin,primMom.Mag(),charge,dEdx,beta,mass2,0);
	  mRunQAProManager->fillRunQA_Track(triggerBin,runIndex,gDCA,nHitsFit,primMom,globMom,0);
	  if( mRunQACut->passTrackQA(picoTrack,mPicoEvent) ) // apply QA track cut
	  {
	    mRunQAHistoManager->fillTrackQA_Kinematics(triggerBin,primMom,globMom, 1); // with track cut
	    mRunQAHistoManager->fillTrackQA_Quliaty(triggerBin,gDCA,nHitsFit,nHitsMax,nHitsDEdx,1);
	    mRunQAHistoManager->fillTrackQA_PID(triggerBin,primMom.Mag(),charge,dEdx,beta,mass2,1);
	    mRunQAProManager->fillRunQA_Track(triggerBin,runIndex,gDCA,nHitsFit,primMom,globMom,1);
	  }
	}
      }
    }
  }

  return kStOK;
}

