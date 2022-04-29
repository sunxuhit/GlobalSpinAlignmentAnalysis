#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"

#include "StRoot/StVecMesonMaker/StVecMesonMaker.h"
#include "StRoot/StVecMesonMaker/StVecMesonCut.h"
#include "StRoot/StVecMesonMaker/StVecMesonHistoManager.h"
#include "StRoot/StVecMesonMaker/StVecMesonUtility.h"
#include "StRoot/StVecMesonMaker/StVecMesonProManager.h"
#include "StRoot/StVecMesonMaker/StVecMesonZdcEpManager.h"
// #include "StRoot/StVecMesonMaker/StVecMesonTree.h"
#include "StRoot/StVecMesonMaker/StVecMesonCons.h"

#include <algorithm>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

ClassImp(StVecMesonMaker)

//-----------------------------------------------------------------------------
StVecMesonMaker::StVecMesonMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int Mode, const int energy, const int flag_ME) : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;
  mMode = Mode;
  mEnergy = energy;
  mFlag_ME = flag_ME;

  if(mMode == 0)
  {
    // mOutPut_QA = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/QA/file_%s_QA_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobId.c_str());
    mOutPut_QA = Form("./file_%s_QA_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 1)
  {
    // mOutPut_ZdcGainCorr = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/GainParameter/file_%s_ZdcGainCorr_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobId.c_str());
    mOutPut_ZdcGainCorr = Form("./file_%s_ZdcGainCorr_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  /*
  if(mMode == 2)
  {
    // mOutPut_ReCenterPar = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/ReCenterParameter/file_%s_ReCenterPar_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobId.c_str());
    mOutPut_ReCenterPar = Form("./file_%s_ReCenterPar_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
  if(mMode == 3)
  {
    // mOutPut_ShiftPar = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/ShiftParameter/file_%s_ShiftPar_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobId.c_str()); 
    mOutPut_ShiftPar = Form("./file_%s_ShiftPar_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobId.c_str()); 
  }
  if(mMode == 4)
  {
    // mOutPut_Resolution = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Resolution/file_%s_Resolution_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobId.c_str()); 
    mOutPut_Resolution = Form("./file_%s_Resolution_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobId.c_str()); 
  }
  if(mMode == 5)
  {
    // mOutPut_Phi = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/Phi/Forest/file_%s_Phi_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[mFlag_ME].Data(),jobId.c_str()); 
    mOutPut_Phi = Form("./file_%s_Phi_%s_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[mFlag_ME].Data(),jobId.c_str()); 
  }
  */
}

//----------------------------------------------------------------------------- 
StVecMesonMaker::~StVecMesonMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StVecMesonMaker::Init() 
{
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mVecMesonCut = new StVecMesonCut(mEnergy);
  mVecMesonHistoManager = new StVecMesonHistoManager();
  mVecMesonUtility = new StVecMesonUtility(mEnergy);
  mVecMesonUtility->initRunIndex(); // initialize std::map for run index
  mVecMesonProManager = new StVecMesonProManager();
  mVecMesonZdcEpManager = new StVecMesonZdcEpManager(mEnergy); // initialize ZDC EP Manager

  if(mMode == 0)
  { // QA
    mFile_QA= new TFile(mOutPut_QA.c_str(),"RECREATE");
    mFile_QA->cd();
    mVecMesonHistoManager->initEventQA();
    mVecMesonHistoManager->initTrackQA();
    mVecMesonProManager->initRunQA();
  }

  if(mMode == 1)
  { // fill Gain Correction Factors for BBC & ZDC
    mFile_ZdcGainCorr= new TFile(mOutPut_ZdcGainCorr.c_str(),"RECREATE");
    mVecMesonHistoManager->initZdcGainCorr();
    // mVecMesonProManager->InitReCenter();
  }

  /*
  if(mMode == 2)
  {
    mFile_ReCenterPar = new TFile(mOutPut_ReCenterPar.c_str(),"RECREATE");
    mFile_ReCenterPar->cd();
    // mVecMesonProManager->InitReCenter();
  }

  if(mMode == 3)
  {
    mFile_ShiftPar = new TFile(mOutPut_ShiftPar.c_str(),"RECREATE");
    mUsedTrackCounter = 0;
    // mVecMesonCorrection->InitReCenterCorrection();
    // mVecMesonProManager->InitShift();
  }

  if(mMode == 4)
  {
    mFile_Resolution = new TFile(mOutPut_Resolution.c_str(),"RECREATE");
    // mVecMesonCorrection->InitReCenterCorrection();
    // mVecMesonCorrection->InitShiftCorrection();
    // mVecMesonProManager->InitResolution();
    // mVecMesonHistoManager->InitEP();
  }

  if(mMode == 5)
  {
    // mVecMesonTree = new StVecMesonTree(mEnergy);
    mFile_Phi = new TFile(mOutPut_Phi.c_str(),"RECREATE");
    mFile_Phi->cd();
    // mVecMesonTree->InitPhi();
    // mVecMesonCorrection->InitReCenterCorrection();
    // mVecMesonCorrection->InitShiftCorrection();
    // mVecMesonCorrection->InitResolutionCorr();
  }
  */

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StVecMesonMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_QA != "")
    {
      mFile_QA->cd();
      mVecMesonHistoManager->writeEventQA();
      mVecMesonHistoManager->writeTrackQA();
      mVecMesonProManager->writeRunQA();
      mFile_QA->Close();
    }
  }
  if(mMode == 1)
  {
    if(mOutPut_ZdcGainCorr != "")
    {
      mFile_ZdcGainCorr->cd();
      mVecMesonHistoManager->writeZdcGainCorr();
      mFile_ZdcGainCorr->Close();
    }
  }
  /*
  if(mMode == 2)
  {
    if(mOutPut_ReCenterPar != "")
    {
      mFile_ReCenterPar->cd();
      // mVecMesonProManager->WriteReCenter();
      mVecMesonHistoManager->WriteQA();
      mFile_ReCenterPar->Close();
    }
  }
  if(mMode == 3)
  {
    if(mOutPut_ShiftPar != "")
    {
      mFile_ShiftPar->cd();
      // mVecMesonProManager->WriteShift();
      mFile_ShiftPar->Close();
    }
  }
  if(mMode == 4)
  {
    if(mOutPut_Resolution != "")
    {
      mFile_Resolution->cd();
      // mVecMesonHistoManager->WriteEP();
      // mVecMesonProManager->WriteResolution();
      mFile_Resolution->Close();
    }
  }
  if(mMode == 5)
  {
    if(mOutPut_Phi != "")
    {
      mFile_Phi->cd();
      // mVecMesonTree->WritePhiMass2();
      mFile_Phi->Close();
    }
  }
  */

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StVecMesonMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
int StVecMesonMaker::Make() 
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
  if( mVecMesonCut->isMinBias(mPicoEvent) )
  {
    // Event Information
    const int runId = mPicoEvent->runId();
    const int refMult = mPicoEvent->refMult();
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
    mRefMultCorr->initEvent(refMult,vz,zdcX); // BES-II might need Luminosity corrections
    // if( !mVecMesonCut->isBES(mEnergy) ) mRefMultCorr->initEvent(refMult,vz,zdcX); // for 200 GeV
    // if( mVecMesonCut->isBES(mEnergy) ) mRefMultCorr->initEvent(refMult,vz,0.0); // for BES Energy

    if(mRefMultCorr->isBadRun(runId))
    {
      LOG_ERROR << "Bad Run from StRefMultCorr! Skip!" << endm;
      return kStErr;
    }

    // vz sign
    int vz_sign = 0; // 0 for -vz || 1 for vz
    vz > 0.0 ? vz_sign = 1 : vz_sign = 0;

    const int cent9 = mRefMultCorr->getCentralityBin9(); // get Centrality9
    const int runIndex = mVecMesonUtility->findRunIndex(runId); // find run index for a specific run
    // cout << "runId = " << runId << ", runIndex = " << runIndex << endl;
    if(runIndex < 0)
    {
      LOG_ERROR << "Could not find this run Index from StVecMesonUtility! Skip!" << endm;
      return kStErr;
    }

    if(mMode == 0)
    { // fill QA before event cuts
      mVecMesonHistoManager->fillEventQA_RefMult(refMult,cent9,numOfBTofHits,numOfBTofMatch,0); // wo event cut
      mVecMesonHistoManager->fillEventQA_Vertex(vx,vy,vz,vzVpd,0);
      mVecMesonProManager->fillRunQA_Event(runIndex,refMult,zdcX,vx,vy,vz,0);

      if(mVecMesonCut->passEventCut(mPicoDst))
      { // apply Event Cuts for anlaysis 
	mVecMesonHistoManager->fillEventQA_RefMult(refMult,cent9,numOfBTofHits,numOfBTofMatch,1); // with event cut
	mVecMesonHistoManager->fillEventQA_Vertex(vx,vy,vz,vzVpd,1);
	mVecMesonProManager->fillRunQA_Event(runIndex,refMult,zdcX,vx,vy,vz,1);

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
	  float beta       = mVecMesonCut->getBeta(mPicoDst,i_track);
	  float mass2      = mVecMesonCut->getPrimaryMass2(mPicoDst,i_track);

	  mVecMesonHistoManager->fillTrackQA_Kinematics(primMom,globMom, 0); // wo track cut
	  mVecMesonHistoManager->fillTrackQA_Quliaty(gDCA,nHitsFit,nHitsMax,nHitsDEdx,0);
	  mVecMesonHistoManager->fillTrackQA_PID(primMom.Mag(),charge,dEdx,beta,mass2,0);
	  mVecMesonProManager->fillRunQA_Track(runIndex,gDCA,nHitsFit,primMom,globMom,0);
	  if( mVecMesonCut->passTrackQA(picoTrack) ) // apply QA track cut
	  {
	    mVecMesonHistoManager->fillTrackQA_Kinematics(primMom,globMom, 1); // with track cut
	    mVecMesonHistoManager->fillTrackQA_Quliaty(gDCA,nHitsFit,nHitsMax,nHitsDEdx,1);
	    mVecMesonHistoManager->fillTrackQA_PID(primMom.Mag(),charge,dEdx,beta,mass2,1);
	    mVecMesonProManager->fillRunQA_Track(runIndex,gDCA,nHitsFit,primMom,globMom,1);
	  }
	}
      }
    }
    if(mMode == 1)
    { // fill Gain Correction Factors for BBC & ZDC
      if(mVecMesonCut->passEventCut(mPicoDst))
      { // apply Event Cuts for anlaysis 
	for(int i_slat = 0; i_slat < 8; ++i_slat) // read in raw ADC value from ZDC-SMD
	{
	  mVecMesonZdcEpManager->setZdcSmd(0,0,i_slat,mPicoEvent->ZdcSmdEastVertical(i_slat));
	  mVecMesonZdcEpManager->setZdcSmd(0,1,i_slat,mPicoEvent->ZdcSmdEastHorizontal(i_slat));
	  mVecMesonZdcEpManager->setZdcSmd(1,0,i_slat,mPicoEvent->ZdcSmdWestVertical(i_slat));
	  mVecMesonZdcEpManager->setZdcSmd(1,1,i_slat,mPicoEvent->ZdcSmdWestHorizontal(i_slat));
	}
	for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest) // fill ZDC Gain Correction Histograms
	{
	  for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
	  {
	    for(int i_slat = 0; i_slat < 8; ++i_slat)
	    {
	      mVecMesonHistoManager->fillZdcGainCorr(i_eastwest,i_verthori,i_slat,runIndex,mVecMesonZdcEpManager->getZdcSmd(i_eastwest,i_verthori,i_slat));
	      // cout << "i_eastwest = " << i_eastwest << ", i_verthori = " << i_verthori << ", i_slat = " << i_slat << ", zdc = " << mZdcSmdCorrection->getZdcSmd(i_eastwest,i_verthori,i_slat) << endl;
	    }
	  }
	}
      }
    }
  }

#if 0
  if(mVecMesonCut->passEventCut(mPicoDst))
  {
    const int nTracks = mPicoDst->numberOfTracks();
    const int cent9 = mRefMultCorr->getCentralityBin9();
//    if(cent9 < 0) cout << cent9 << endl;
    const Double_t reweight = mRefMultCorr->getWeight();
    const int nToFMatched = mVecMesonCut->getMatchedToF();

    for(int i = 0; i < nTracks; i++) // track loop
    {
      StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

      if(mMode == 0)
      {
	// float eta = track->pMom().pseudoRapidity();
	float eta = track->pMom().Eta();
	if(fabs(eta) < vmsa::mEtaMax && track->dca() < 3.0)
	{
	  float Mass2 = mVecMesonCut->getMass2(track);
	  float dEdx = track->dEdx();
	  // float p = track->pMom().mag();
	  float p = track->pMom().Mag();
	  mVecMesonHistoManager->FillQA_Detector(dEdx,Mass2,p);
	}
      }
      if(mVecMesonCut->passTrackEP(track)) // track cut
      {
	if(mMode == 0) // fill re-center parameter
	{
	  // float pt = track->pMom().perp();
	  float pt = track->pMom().Perp();

	  if(mVecMesonCorrection->passTrackFull(track))
	  {
	    TVector2 q2Vector_Full = mVecMesonCorrection->calq2Vector(track);
	    mVecMesonProManager->FillTrackFull(q2Vector_Full,cent9,runIndex,vz_sign,pt);
	    mVecMesonCorrection->addTrack_FullRaw(track,cent9,runIndex);
	  }

	  if(mVecMesonCorrection->passTrackEtaEast(track)) // neg eta sub
	  {
	    TVector2 q2Vector_East = mVecMesonCorrection->calq2Vector(track);
	    mVecMesonProManager->FillTrackEast(q2Vector_East,cent9,runIndex,vz_sign,pt);
	    mVecMesonCorrection->addTrack_EastRaw(track,cent9,runIndex);
	  }
	  if(mVecMesonCorrection->passTrackEtaWest(track)) // pos eta sub
	  {
	    TVector2 q2Vector_West = mVecMesonCorrection->calq2Vector(track);
	    mVecMesonProManager->FillTrackWest(q2Vector_West,cent9,runIndex,vz_sign,pt);
	    mVecMesonCorrection->addTrack_WestRaw(track,cent9,runIndex);
	  }
	}
	else // calculate Q Vector after recentering for full event and eta sub event
	{
	  if(mVecMesonCorrection->passTrackFull(track))
	  {
	    mVecMesonCorrection->addTrack_Full(track,cent9,runIndex,vz_sign);
	    mUsedTrackCounter++;
	  }
	  if(mVecMesonCorrection->passTrackEtaEast(track)) // neg eta sub
	  {
	    mVecMesonCorrection->addTrack_East(track,cent9,runIndex,vz_sign);
	  }
	  if(mVecMesonCorrection->passTrackEtaWest(track)) // pos eta sub
	  {
	    mVecMesonCorrection->addTrack_West(track,cent9,runIndex,vz_sign);
	  }
	}
      }
    }

    if(mMode == 0) // fill raw EP
    {
      mVecMesonHistoManager->FillQA_Event(vz,refMult);
      if(mVecMesonCorrection->passTrackEtaNumRawCut())
      {
	TVector2 Q2East = mVecMesonCorrection->getQVectorRaw(0); // 0 = eta_gap, 1 = east/west
	float Psi2_East = 0.5*TMath::ATan2(Q2East.Y(),Q2East.X());
	TVector2 Q2West = mVecMesonCorrection->getQVectorRaw(1); // 0 = eta_gap, 1 = east/west
	float Psi2_West = 0.5*TMath::ATan2(Q2West.Y(),Q2West.X());
	mVecMesonHistoManager->FillEP_Eta(Psi2_East,Psi2_West);
      }
      if(mVecMesonCorrection->passTrackFullNumRawCut())
      {
	TVector2 Q2Full = mVecMesonCorrection->getQVectorRaw(2);
	float Psi2_Full = 0.5*TMath::ATan2(Q2Full.Y(),Q2Full.X());
	mVecMesonHistoManager->FillEP_Full(Psi2_Full);
      }
    }

    if(mMode == 1)
    {
      // full event shift parameter
      if(mVecMesonCorrection->passTrackFullNumCut())
      {
	for(int k = 0; k < 5; k++) // ShiftOrder loop
	{
	  TVector2 Psi2Vector_Full_EP = mVecMesonCorrection->calPsi2_Full_EP(k);
	  mVecMesonProManager->FillEventFull_EP(Psi2Vector_Full_EP,cent9,runIndex,vz_sign,k);
	}
      }

      // eta sub event shift parameter
      if(mVecMesonCorrection->passTrackEtaNumCut())
      {
	for(int k = 0; k < 5; k++)
	{
	  TVector2 Psi2Vector_East_EP = mVecMesonCorrection->calPsi2_East_EP(k);
	  mVecMesonProManager->FillEventEast_EP(Psi2Vector_East_EP,cent9,runIndex,vz_sign,k);

	  TVector2 Psi2Vector_West_EP = mVecMesonCorrection->calPsi2_West_EP(k);
	  mVecMesonProManager->FillEventWest_EP(Psi2Vector_West_EP,cent9,runIndex,vz_sign,k);
	}
      }
    }

    if(mMode == 2) // calculate resolution for eta_sub and random sub event plane
    {
      // calculate Q vector after recentering for Random Sub Event
      int iTrack[mUsedTrackCounter];
      float ranCounter = (float)mUsedTrackCounter/2.0 - 1;
      for(int i = 0; i < mUsedTrackCounter; i++)
      {
        iTrack[i] = i;
      }
      std::srand(time(0));
      std::random_shuffle(iTrack,iTrack+mUsedTrackCounter);
      mUsedTrackCounter = 0;
      for(int i = 0; i < nTracks; i++) // track loop
      {
	StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
	if(mVecMesonCut->passTrackEP(track)) // track cut
	{
	  if(mVecMesonCorrection->passTrackFull(track))
	  {
	    if((float)iTrack[mUsedTrackCounter] > ranCounter) // Sub Event A
	    {
	      mVecMesonCorrection->addTrack_A(track,cent9,runIndex,vz_sign);
	    }
	    else // Sub Event B
	    {
	      mVecMesonCorrection->addTrack_B(track,cent9,runIndex,vz_sign);
	    }
	    mUsedTrackCounter++;
	  }
	}
      }
      mUsedTrackCounter = 0;

      // calculate resolution
      TVector2 QVecEast = mVecMesonCorrection->getQVector(0);
      float Psi2East_ReCenter = 0.5*TMath::ATan2(QVecEast.Y(),QVecEast.X());
      float Psi2East_Shift = mVecMesonCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign);

      TVector2 QVecWest = mVecMesonCorrection->getQVector(1);
      float Psi2West_ReCenter = 0.5*TMath::ATan2(QVecWest.Y(),QVecWest.X());
      float Psi2West_Shift = mVecMesonCorrection->calShiftAngle2West_EP(runIndex,cent9,vz_sign);

      TVector2 QVecFull = mVecMesonCorrection->getQVector(2);
      float Psi2Full_ReCenter = 0.5*TMath::ATan2(QVecFull.Y(),QVecFull.X());
      float Psi2Full_Shift = mVecMesonCorrection->calShiftAngle2Full_EP(runIndex,cent9,vz_sign);

      TVector2 QVecRanA = mVecMesonCorrection->getQVector(3);
      float Psi2RanA_ReCenter = 0.5*TMath::ATan2(QVecRanA.Y(),QVecRanA.X());
      float Psi2RanA_Shift = mVecMesonCorrection->calShiftAngle2A_EP(runIndex,cent9,vz_sign);

      TVector2 QVecRanB = mVecMesonCorrection->getQVector(4);
      float Psi2RanB_ReCenter = 0.5*TMath::ATan2(QVecRanB.Y(),QVecRanB.X());
      float Psi2RanB_Shift = mVecMesonCorrection->calShiftAngle2B_EP(runIndex,cent9,vz_sign);

      if(mVecMesonCorrection->passTrackEtaNumCut())
      {
	mVecMesonHistoManager->FillEP_Sub(Psi2East_ReCenter,Psi2East_Shift,Psi2West_ReCenter,Psi2West_Shift);
	mVecMesonProManager->FillRes_Sub(cent9,Psi2East_Shift,Psi2West_Shift);
      }

      if(mVecMesonCorrection->passTrackFullNumCut())
      {
	mVecMesonHistoManager->FillEP_Ran(Psi2RanA_ReCenter,Psi2RanA_Shift,Psi2RanB_ReCenter,Psi2RanB_Shift,Psi2Full_ReCenter,Psi2Full_Shift);
	mVecMesonProManager->FillRes_Ran(cent9,Psi2RanA_Shift,Psi2RanB_Shift);
      }
    }

    if(mMode == 3)
    { // phi meson
      if(mVecMesonCorrection->passTrackFullNumCut())
      {
	// get QVector of sub event
	TVector2 Q2East = mVecMesonCorrection->getQVector(0); // east
	TVector2 Q2West = mVecMesonCorrection->getQVector(1); // west
	TVector2 Q2Full = mVecMesonCorrection->getQVector(2); // full 
	int NumTrackEast = mVecMesonCorrection->getNumTrack(0);
	int NumTrackWest = mVecMesonCorrection->getNumTrack(1);
	int NumTrackFull = mVecMesonCorrection->getNumTrack(2);
	int NumTrackFullEast = mVecMesonCorrection->getNumTrack(3);
	int NumTrackFullWest = mVecMesonCorrection->getNumTrack(4);

	float Psi2 = mVecMesonCorrection->calShiftAngle2Full_EP(runIndex,cent9,vz_sign);

	// get N_prim, N_non_prim, N_Tof_match
	int N_prim = mVecMesonCut->getNpirm();
	int N_non_prim = mVecMesonCut->getNnonprim();
	int N_Tof_match = mVecMesonCut->getMatchedToF();

	// pass the event information to StVecMesonTree
	mVecMesonTree->clearEvent();
	mVecMesonTree->passEvent(N_prim, N_non_prim, N_Tof_match);

	// pass re-centered event plane to StVecMesonTree
	mVecMesonTree->passEventPlane(Q2East,Q2West,Q2Full);

	// pass NumOfTrack to StVecMesonTree
	mVecMesonTree->passNumTrack(NumTrackEast,NumTrackWest,NumTrackFull,NumTrackFullEast,NumTrackFullWest);

	mVecMesonTree->MixEvent_Phi(mFlag_ME,mPicoDst,cent9,vz,Psi2);
      }
    }
    mVecMesonCorrection->clear();
  }
#endif

  return kStOK;
}

