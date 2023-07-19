#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TVector3.h"
#include "TString.h"

#include "StMessMgr.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StPileupUtil/StPileupUtil.h"

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
  mPicoDstMaker      = picoMaker;
  mPicoDst           = NULL;
  mRefMultCorr       = NULL;
  mPileupUtilFxt3p85 = NULL;

  str_mOutPutRunQA = Form("./file_RunQA_%s_%s.root",globCons::str_mBeamType[mType].c_str(),jobId.c_str());
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
    if( mAnaCut->isFxt3p85GeV_2018() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  if(mAnaCut->isFxt3p85GeV_2018() && !mPileupUtilFxt3p85)
  {
    mPileupUtilFxt3p85 = new StPileupUtil();
    mPileupUtilFxt3p85->init();
  }

  file_mOutPutRunQA = new TFile(str_mOutPutRunQA.c_str(),"RECREATE");
  file_mOutPutRunQA->cd();
  mRunQAHistoManager->initEvtQA();
  mRunQAHistoManager->initTrkQA();
  mRunQAProManager->initRunQA();

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StRunQAMaker::Finish() 
{
  if(str_mOutPutRunQA != "")
  {
    file_mOutPutRunQA->cd();
    mRunQAHistoManager->writeEvtQA();
    mRunQAHistoManager->writeTrkQA();
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

  // Event Information
  const int runId        = mPicoEvent->runId();
  int refMult            = mPicoEvent->refMult(); 
  const int grefMult     = mPicoEvent->grefMult();
  const TVector3 primVtx = mPicoEvent->primaryVertex();
  const double vx        = mPicoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const double vy        = mPicoEvent->primaryVertex().y();
  const double vz        = mPicoEvent->primaryVertex().z();
  const double vzVpd     = mPicoEvent->vzVpd();
  const double zdcX      = mPicoEvent->ZDCx();
  const double vxReCtr   = mAnaUtils->getVxReCtr(vx);
  const double vyReCtr   = mAnaUtils->getVyReCtr(vy);
  const unsigned short nBTofMatch = mPicoEvent->nBTOFMatch(); // get number of tof match points
  const unsigned int nBTofHits    = mPicoDst->numberOfBTofHits(); // get number of tof hits
  const unsigned int nTracks      = mPicoDst->numberOfTracks(); // get number of tracks

  const int triggerBin  = mAnaUtils->getTriggerBin(mPicoEvent);
  const int vzBin       = mAnaUtils->getVzBin(vz); // 0 for [vzMin,vzCtr) || 1 for [vzCtr,vzMax]
  const int runIndex    = mAnaUtils->findRunIndex(runId); // find run index for a specific run
  // cout << "runId = " << runId << ", runIndex = " << runIndex << endl;
  mRunQAHistoManager->fillEvtQaVertexAllTrigs(vx,vy,vz);

  // MinBias trigger
  if( mAnaCut->isMinBias(mPicoEvent) )
  {
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StAnalysisUtils! Skip!" << endm;
      return kStErr;
    }

    if(!mRefMultCorr) // StRefMultCorr Cut & centrality
    {
      LOG_WARN << " No mRefMultCorr! Skip! " << endl;
      return kStErr;
    }
    mRefMultCorr->init(runId);
    mRefMultCorr->initEvent(refMult,vz,zdcX);
    int cent9     = mRefMultCorr->getCentralityBin9(); // get Centrality9
    double refWgt = mRefMultCorr->getWeight(); // get Centrality reweight

    if(mAnaCut->isFxt3p85GeV_2018()) // only for Fxt3p85GeV_2018
    { 
      mPileupUtilFxt3p85->initEvent(mPicoDst);
      refMult = mPileupUtilFxt3p85->get_refMultPrim();
      cent9   = mPileupUtilFxt3p85->get_centrality9();
      refWgt  = mPileupUtilFxt3p85->get_centralityWeight();
    }

    bool isPileUpStAnalysisCut = mAnaCut->isPileUpEvent((double)refMult, (double)nBTofMatch,vz); // alway return false for Isobar & Fxt3p85GeV_2018
    bool isPileUpStRefMultCorr = !mRefMultCorr->passnTofMatchRefmultCut((double)refMult, (double)nBTofMatch,vz); // valid for Isobar
    if(mAnaCut->isFxt3p85GeV_2018()) isPileUpStRefMultCorr = mPileupUtilFxt3p85->isPileupEPD(); // valid only for Fxt3p85GeV_2018
    bool isPileUpEvent = isPileUpStAnalysisCut || isPileUpStRefMultCorr;
    // cout << "isPileUpEvent = " << isPileUpEvent << ", isPileUpStAnalysisCut = " << isPileUpStAnalysisCut << ", isPileUpStRefMultCorr = " << isPileUpStRefMultCorr << endl;

    // fill QA before event cuts
    mRunQAProManager->fillRunQAEvt(triggerBin,runIndex,refMult,grefMult,zdcX,vxReCtr,vyReCtr,vz,0);
    mRunQAHistoManager->fillEvtQaRefMult(triggerBin,refMult,grefMult,cent9,refWgt,nBTofHits,nBTofMatch,0); // wo event cut
    mRunQAHistoManager->fillEvtQaVertex(triggerBin,vxReCtr,vyReCtr,vz,vzVpd,vzBin,0);
    mRunQAHistoManager->fillEvtQaTrigger(triggerBin,0);

    if(!isPileUpEvent && mAnaCut->isGoodCent9(cent9) && mAnaCut->passEventCut(mPicoEvent))
    { // apply Event Cuts for anlaysis 
      mRunQAProManager->fillRunQAEvt(triggerBin,runIndex,refMult,grefMult,zdcX,vxReCtr,vyReCtr,vz,1);
      mRunQAHistoManager->fillEvtQaRefMult(triggerBin,refMult,grefMult,cent9,refWgt,nBTofHits,nBTofMatch,1); // with event cut
      mRunQAHistoManager->fillEvtQaVertex(triggerBin,vxReCtr,vyReCtr,vz,vzVpd,vzBin,1);
      mRunQAHistoManager->fillEvtQaTrigger(triggerBin,1);

      for(unsigned int iTrack = 0; iTrack < nTracks; iTrack++) // track loop
      {
	StPicoTrack *picoTrack = (StPicoTrack*)mPicoDst->track(iTrack);
	if(!picoTrack) continue;

	const TVector3 primMom = picoTrack->pMom();
	const TVector3 globMom = picoTrack->gMom();

	const double gDCA     = picoTrack->gDCA(vx,vy,vz);
	const int nHitsFit    = picoTrack->nHitsFit();
	const int nHitsMax    = picoTrack->nHitsMax();
	const int nHitsDEdx   = picoTrack->nHitsDedx();
	const double dEdx     = picoTrack->dEdx();
	const short charge    = picoTrack->charge();
	const double nSigKaon = picoTrack->nSigmaKaon();
	const double beta     = mAnaUtils->getBeta(mPicoDst,iTrack);
	const double mass2    = mAnaUtils->getPrimMass2(mPicoDst,iTrack);

	mRunQAHistoManager->fillTrkQaKinematics(triggerBin,primMom,globMom, 0); // wo track cut
	mRunQAHistoManager->fillTrkQaQuliaty(triggerBin,gDCA,nHitsFit,nHitsMax,nHitsDEdx,0);
	mRunQAHistoManager->fillTrkQaPID(triggerBin,primMom.Mag(),charge,dEdx,beta,mass2,0);
	mRunQAProManager->fillRunQATrk(triggerBin,runIndex,gDCA,nHitsFit,primMom,globMom,0);
	if( mAnaCut->passTrkQA(picoTrack,primVtx) ) // apply QA track cut
	{
	  mRunQAHistoManager->fillTrkQaKinematics(triggerBin,primMom,globMom, 1); // with track cut
	  mRunQAHistoManager->fillTrkQaQuliaty(triggerBin,gDCA,nHitsFit,nHitsMax,nHitsDEdx,1);
	  mRunQAHistoManager->fillTrkQaPID(triggerBin,primMom.Mag(),charge,dEdx,beta,mass2,1);
	  mRunQAProManager->fillRunQATrk(triggerBin,runIndex,gDCA,nHitsFit,primMom,globMom,1);
	}

	bool isEpFull = mAnaCut->passTrkTpcEpFull(picoTrack,primVtx);
	bool isEpEast = mAnaCut->passTrkTpcEpEast(picoTrack,primVtx);
	bool isEpWest = mAnaCut->passTrkTpcEpWest(picoTrack,primVtx);
	bool isFlowFull = mAnaCut->passTrkTpcFlowFull(picoTrack,primVtx);
	bool isFlowEast = mAnaCut->passTrkTpcFlowEast(picoTrack,primVtx);
	bool isFlowWest = mAnaCut->passTrkTpcFlowWest(picoTrack,primVtx);
	bool isKaonFull = mAnaCut->passTrkTpcKaonFull(picoTrack,primVtx);
	bool isKaonEast = mAnaCut->passTrkTpcKaonEast(picoTrack,primVtx);
	bool isKaonWest = mAnaCut->passTrkTpcKaonWest(picoTrack,primVtx);

	mRunQAHistoManager->fillTrkQaEpCut(triggerBin, primMom, isEpFull, isEpEast, isEpWest, 0);
	mRunQAHistoManager->fillTrkQaFlowCut(triggerBin, primMom, isFlowFull, isFlowEast, isFlowWest, 0);
	mRunQAHistoManager->fillTrkQaKaonCut(triggerBin, primMom, nSigKaon, isKaonFull, isKaonEast, isKaonWest, 0);

	mRunQAHistoManager->fillTrkQaEpCut(triggerBin, primMom, isEpFull, isEpEast, isEpWest, 1);
	mRunQAHistoManager->fillTrkQaFlowCut(triggerBin, primMom, isFlowFull, isFlowEast, isFlowWest, 1);
	mRunQAHistoManager->fillTrkQaKaonCut(triggerBin, primMom, nSigKaon, isKaonFull, isKaonEast, isKaonWest, 1);

	if( mAnaCut->passTrkTpcKaonFull(picoTrack, primVtx) )
	{
	  if(charge > 0) // K+
	  {
	    const double yLab = mAnaUtils->getRapidityLab(picoTrack,321);
	    const double yCms = mAnaUtils->getRapidityCMS(yLab);
	    if(mAnaCut->passTrkTofKaonMass(primMom,charge,mass2))
	    {
	      mRunQAHistoManager->fillTrkQaKaonAcptTree(triggerBin, cent9, charge, yLab, yCms, primMom.Pt(), refWgt);
	    }
	    if( mass2 > 0.16 && mass2 < 0.36 ) // temp fix
	    {
	      mRunQAHistoManager->fillTrkQaKaonAcptSpin(triggerBin, cent9, charge, yLab, yCms, primMom.Pt(), refWgt);
	    }
	  }
	  if(charge < 0) // K-
	  {
	    const double yLab = mAnaUtils->getRapidityLab(picoTrack,-321);
	    const double yCms = mAnaUtils->getRapidityCMS(yLab);
	    if(mAnaCut->passTrkTofKaonMass(primMom,charge,mass2))
	    {
	      mRunQAHistoManager->fillTrkQaKaonAcptTree(triggerBin, cent9, charge, yLab, yCms, primMom.Pt(), refWgt);
	    }
	    if( mass2 > 0.16 && mass2 < 0.36 ) // temp fix
	    {
	      mRunQAHistoManager->fillTrkQaKaonAcptSpin(triggerBin, cent9, charge, yLab, yCms, primMom.Pt(), refWgt);
	    }
	  }
	}
      }
    }
  }

  return kStOK;
}

