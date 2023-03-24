#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"

#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

// #include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonEvent.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonTree.h"

ClassImp(StPhiMesonTree)

//------------------------------------------------------------------------------------------------------------------
StPhiMesonTree::StPhiMesonTree(int beamType) : mType(beamType)
{
  // mType = beamType;
}

StPhiMesonTree::~StPhiMesonTree()
{
  /* */
}

//------------------------------------------------------------------------------------------------------------------

void StPhiMesonTree::initPhiTree()
{
  mAnaCut   = new StAnalysisCut(mType);
  mAnaUtils = new StAnalysisUtils(mType);

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mInvMassPhiCent%d",iCent);
    h_mInvMassPhi[iCent] = new TH2F(histName.c_str(),histName.c_str(),20,0.2,5.0,200,0.98,1.08);
    for(int iVz = 0; iVz < mNumMixVzBin; ++iVz)
    {
      for(int iPsi = 0; iPsi < mNumMixPsiBin; ++iPsi)
      {
	clearPhiMixBuffer(iCent,iVz,iPsi);
      }
    }
  }

  mPhiMesonEvent = new StPhiMesonEvent();
  t_mPhiMesonTree = new TTree("PhiMesonEvent","PhiMesonEvent");
  t_mPhiMesonTree->Branch("phiSpinAlignmentBranch","StPhiMesonEvent",&mPhiMesonEvent);
  t_mPhiMesonTree->SetAutoSave(5000000);
}

//------------------------------------------------------------------------------------------------------------------

void StPhiMesonTree::writePhiTree()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mInvMassPhi[iCent]->Write();
  }
  t_mPhiMesonTree->Write("",TObject::kOverwrite);
}

//------------------------------------------------------------------------------------------------------------------

void StPhiMesonTree::clearPhiMixBuffer(int cent9, int vzBin, int PsiBin)
{
  mEventCounter[cent9][vzBin][PsiBin] = 0; // General Event Info
  vec_mRunId[cent9][vzBin][PsiBin].clear();
  vec_mEvtId[cent9][vzBin][PsiBin].clear();
  vec_mRefMult[cent9][vzBin][PsiBin].clear();
  vec_mNumTofMatch[cent9][vzBin][PsiBin].clear();
  vec_mCent9[cent9][vzBin][PsiBin].clear();
  vec_mCent16[cent9][vzBin][PsiBin].clear();
  vec_mRefWgt[cent9][vzBin][PsiBin].clear();
  vec_mZDCx[cent9][vzBin][PsiBin].clear();
  vec_mBBCx[cent9][vzBin][PsiBin].clear();
  vec_mVzVpd[cent9][vzBin][PsiBin].clear();
  vec_mPrimVtx[cent9][vzBin][PsiBin].clear();

  vec_mFlagZdcEp[cent9][vzBin][PsiBin].clear(); // ZDC EP Info
  vec_mQ1ZdcShiftEast[cent9][vzBin][PsiBin].clear();
  vec_mQ1ZdcShiftWest[cent9][vzBin][PsiBin].clear();
  vec_mQ1ZdcShiftFull[cent9][vzBin][PsiBin].clear();

  vec_mFlagEpdEp[cent9][vzBin][PsiBin].clear(); // EPD EP Info
  vec_mQ1EpdShiftEast[cent9][vzBin][PsiBin].clear();
  vec_mQ1EpdShiftWest[cent9][vzBin][PsiBin].clear();
  vec_mQ1EpdShiftFull[cent9][vzBin][PsiBin].clear();

  vec_mFlagTpcEp[cent9][vzBin][PsiBin].clear(); // TPC EP Info
  vec_mQ2TpcReCtrEast[cent9][vzBin][PsiBin].clear();
  vec_mQ2TpcReCtrWest[cent9][vzBin][PsiBin].clear();
  vec_mQ3TpcReCtrEast[cent9][vzBin][PsiBin].clear();
  vec_mQ3TpcReCtrWest[cent9][vzBin][PsiBin].clear();
  vec_mNumTrkReCtrEast[cent9][vzBin][PsiBin].clear();
  vec_mNumTrkReCtrWest[cent9][vzBin][PsiBin].clear();

  for(int evtBin = 0; evtBin < mNumMixBuffer; ++evtBin)
  {
    int phiMixKey = getPhiMixKey(cent9,vzBin,PsiBin,evtBin);
    map_mMomVecKp[phiMixKey].clear();
    map_mMomVecKm[phiMixKey].clear();
    map_mMass2Kp[phiMixKey].clear();
    map_mMass2Km[phiMixKey].clear();
    map_mNSigKp[phiMixKey].clear();
    map_mNSigKm[phiMixKey].clear();
    map_mDcaKp[phiMixKey].clear();
    map_mDcaKm[phiMixKey].clear();
    map_mChargeKp[phiMixKey].clear();
    map_mChargeKm[phiMixKey].clear();
    map_mNHitsFitKp[phiMixKey].clear();
    map_mNHitsFitKm[phiMixKey].clear();
  }
}
//------------------------------------------------------------------------------------------------------------------
void StPhiMesonTree::fillPhiTree(StPicoDst *picoDst, int flagME)
{
  StPicoEvent *picoEvent = (StPicoEvent*)picoDst->event();

  int vzBin  = getVzMixBin(mVz);
  int PsiBin = getPsi2MixBin(mPsi2ShiftFull);
  int evtBin = mEventCounter[mCent9][vzBin][PsiBin];

  const unsigned int nTracks = picoDst->numberOfTracks();
  TVector3 primVtx = picoEvent->primaryVertex();

  // store Enent Information
  mEventCounter[mCent9][vzBin][PsiBin]++; // General Event Info
  vec_mRunId[mCent9][vzBin][PsiBin].push_back(static_cast<int>(picoEvent->runId()));
  vec_mEvtId[mCent9][vzBin][PsiBin].push_back(static_cast<int>(picoEvent->eventId()));
  vec_mRefMult[mCent9][vzBin][PsiBin].push_back(static_cast<int>(picoEvent->refMult()));
  vec_mNumTofMatch[mCent9][vzBin][PsiBin].push_back(static_cast<int>(picoEvent->nBTOFMatch()));
  vec_mCent9[mCent9][vzBin][PsiBin].push_back(static_cast<int>(mCent9));
  vec_mCent16[mCent9][vzBin][PsiBin].push_back(static_cast<int>(mCent16));
  vec_mRefWgt[mCent9][vzBin][PsiBin].push_back(static_cast<double>(mRefWgt));
  vec_mZDCx[mCent9][vzBin][PsiBin].push_back(static_cast<double>(picoEvent->ZDCx()));
  vec_mBBCx[mCent9][vzBin][PsiBin].push_back(static_cast<double>(picoEvent->BBCx()));
  vec_mVzVpd[mCent9][vzBin][PsiBin].push_back(static_cast<double>(picoEvent->vzVpd()));
  vec_mPrimVtx[mCent9][vzBin][PsiBin].push_back(static_cast<TVector3>(primVtx));

  vec_mFlagZdcEp[mCent9][vzBin][PsiBin].push_back(static_cast<int>(mFlagZdcEp)); // ZDC EP Info
  vec_mQ1ZdcShiftEast[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ1VecZdcShiftEast));
  vec_mQ1ZdcShiftWest[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ1VecZdcShiftWest));
  vec_mQ1ZdcShiftFull[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ1VecZdcShiftFull));

  vec_mFlagEpdEp[mCent9][vzBin][PsiBin].push_back(static_cast<int>(mFlagEpdEp)); // EPD EP Info
  vec_mQ1EpdShiftEast[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ1VecEpdShiftEast));
  vec_mQ1EpdShiftWest[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ1VecEpdShiftWest));
  vec_mQ1EpdShiftFull[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ1VecEpdShiftFull));

  vec_mFlagTpcEp[mCent9][vzBin][PsiBin].push_back(static_cast<int>(mFlagTpcEp)); // TPC EP Info
  vec_mQ2TpcReCtrEast[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ2VecTpcReCtrEast));
  vec_mQ2TpcReCtrWest[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ2VecTpcReCtrWest));
  vec_mQ3TpcReCtrEast[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ3VecTpcReCtrEast));
  vec_mQ3TpcReCtrWest[mCent9][vzBin][PsiBin].push_back(static_cast<TVector2>(v_mQ3VecTpcReCtrWest));
  vec_mNumTrkReCtrEast[mCent9][vzBin][PsiBin].push_back(static_cast<int>(mNumTrkReCtrEast));
  vec_mNumTrkReCtrWest[mCent9][vzBin][PsiBin].push_back(static_cast<int>(mNumTrkReCtrWest));

  // store Track Information
  for(unsigned int iTrk = 0; iTrk < nTracks; ++iTrk) // loop over all particles in event
  {
    StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(iTrk);

    if(mAnaCut->passTrkKaonFull(picoTrack, primVtx))
    {
      double mass2 = mAnaUtils->getPrimMass2(picoDst, iTrk);
      int charge = static_cast<int>(picoTrack->charge());
      TVector3 primMom = picoTrack->pMom();

      if(
	  (primMom.Mag() < 0.65 && ((mass2 > anaUtils::mMass2KaonMin[mType] && mass2 < anaUtils::mMass2KaonMax[mType]) || mass2 < -10.0)) // dE/dx + ToF
	  || (primMom.Mag() >= 0.65 && (mass2 > anaUtils::mMass2KaonMin[mType] && mass2 < anaUtils::mMass2KaonMax[mType])) // dE/dx + ToF(always)
	)
      {
	int phiMixKey = getPhiMixKey(mCent9,vzBin,PsiBin,evtBin);
	if(charge > 0)
	{ // K+ candidate
	  map_mMomVecKp[phiMixKey].push_back(static_cast<TVector3>(picoTrack->pMom()));// primMom 
	  map_mMass2Kp[phiMixKey].push_back(static_cast<double>(mass2)); // mass2
	  map_mNSigKp[phiMixKey].push_back(static_cast<double>(picoTrack->nSigmaKaon())); // nSigmaKaon
	  map_mDcaKp[phiMixKey].push_back(static_cast<double>(picoTrack->gDCA(primVtx.X(),primVtx.Y(),primVtx.Z()))); // dca
	  map_mChargeKp[phiMixKey].push_back(static_cast<int>(charge)); // charge
	  map_mNHitsFitKp[phiMixKey].push_back(static_cast<double>(picoTrack->nHitsFit())); // nHitsFit
	}
	if(charge < 0)
	{ // K- candidate
	  map_mMomVecKm[phiMixKey].push_back(static_cast<TVector3>(picoTrack->pMom()));// primMom 
	  map_mMass2Km[phiMixKey].push_back(static_cast<double>(mass2)); // mass2
	  map_mNSigKm[phiMixKey].push_back(static_cast<double>(picoTrack->nSigmaKaon())); // nSigmaKaon
	  map_mDcaKm[phiMixKey].push_back(static_cast<double>(picoTrack->gDCA(primVtx.X(),primVtx.Y(),primVtx.Z()))); // dca
	  map_mChargeKm[phiMixKey].push_back(static_cast<int>(charge)); // charge
	  map_mNHitsFitKm[phiMixKey].push_back(static_cast<double>(picoTrack->nHitsFit())); // nHitsFit
	}
      }
    }
  }


  if(flagME == 0) // same event
  {
    recoPhi(mCent9,vzBin,PsiBin);
    clearPhiMixBuffer(mCent9,vzBin,PsiBin);
  }

  if(flagME == 1 && mEventCounter[mCent9][vzBin][PsiBin] == mNumMixBuffer) // mix event
  {
    mixPhi(mCent9,vzBin,PsiBin);
    clearPhiMixBuffer(mCent9,vzBin,PsiBin);
  }
}

void StPhiMesonTree::recoPhi(int cent9, int vzBin, int PsiBin) // reconstruct phi meson in the same event
{
  for(int iEvt = 0; iEvt < mEventCounter[cent9][vzBin][PsiBin]; ++iEvt)
  {
    int evtBin = iEvt;
    mPhiMesonEvent->clearTrackList();
    mPhiMesonEvent->setRunId(vec_mRunId[cent9][vzBin][PsiBin][evtBin]); // event header
    mPhiMesonEvent->setEvtId(vec_mEvtId[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setRefMult(vec_mRefMult[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setNumTofMatch(vec_mNumTofMatch[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setCentrality9(vec_mCent9[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setCentrality16(vec_mCent16[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setRefWgt(vec_mRefWgt[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setZDCx(vec_mZDCx[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setBBCx(vec_mBBCx[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setVzVpd(vec_mVzVpd[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setPrimVtx(vec_mPrimVtx[cent9][vzBin][PsiBin][evtBin]);

    mPhiMesonEvent->setFlagZdcEp(vec_mFlagZdcEp[cent9][vzBin][PsiBin][evtBin]); // ZDC EP Info
    mPhiMesonEvent->setQ1VecZdcEast(vec_mQ1ZdcShiftEast[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setQ1VecZdcWest(vec_mQ1ZdcShiftWest[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setQ1VecZdcFull(vec_mQ1ZdcShiftFull[cent9][vzBin][PsiBin][evtBin]);

    mPhiMesonEvent->setFlagEpdEp(vec_mFlagEpdEp[cent9][vzBin][PsiBin][evtBin]); // EPD EP Info
    mPhiMesonEvent->setQ1VecEpdEast(vec_mQ1EpdShiftEast[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setQ1VecEpdWest(vec_mQ1EpdShiftWest[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setQ1VecEpdFull(vec_mQ1EpdShiftFull[cent9][vzBin][PsiBin][evtBin]);

    mPhiMesonEvent->setFlagTpcEp(vec_mFlagTpcEp[cent9][vzBin][PsiBin][evtBin]); // TPC EP Info
    mPhiMesonEvent->setQ2VecTpcEast(vec_mQ2TpcReCtrEast[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setQ2VecTpcWest(vec_mQ2TpcReCtrWest[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setQ3VecTpcEast(vec_mQ3TpcReCtrEast[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setQ3VecTpcWest(vec_mQ3TpcReCtrWest[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setNumTrkReCtrEast(vec_mNumTrkReCtrEast[cent9][vzBin][PsiBin][evtBin]);
    mPhiMesonEvent->setNumTrkReCtrWest(vec_mNumTrkReCtrWest[cent9][vzBin][PsiBin][evtBin]);

    // start to select phi candidate in the same event
    int phiMixKey = getPhiMixKey(cent9,vzBin,PsiBin,evtBin);

    TLorentzVector lTrkKp, lTrkKm;
    for(unsigned int iTrkKp = 0; iTrkKp < map_mMomVecKp[phiMixKey].size(); ++iTrkKp)
    { // first track loop over K+ candidates
      TVector3 primMomKp = map_mMomVecKp[phiMixKey][iTrkKp];
      lTrkKp.SetXYZM(primMomKp.X(),primMomKp.Y(),primMomKp.Z(), anaUtils::mMassKaon);

      for(unsigned int iTrkKm = 0; iTrkKm < map_mMomVecKm[phiMixKey].size(); ++iTrkKm)
      { // second track loop over K- candidates
	TVector3 primMomKm = map_mMomVecKm[phiMixKey][iTrkKm];
	lTrkKm.SetXYZM(primMomKm.X(),primMomKm.Y(),primMomKm.Z(),anaUtils::mMassKaon);

	TLorentzVector lTrkPhi = lTrkKp+lTrkKm;
	double invMassPhi = lTrkPhi.M();
	double ptPhi = lTrkPhi.Perp();

	// fill phi candidate (mass within [0.98, 1.08]) into t_mPhiMesonTree
	if(invMassPhi >= anaUtils::mMassPhiMin && invMassPhi <= anaUtils::mMassPhiMax) 
	{
	  mPhiMesonTrack = mPhiMesonEvent->createTrack();
	  mPhiMesonTrack->setTrkMomKp(map_mMomVecKp[phiMixKey][iTrkKp]); // K+
	  mPhiMesonTrack->setTrkMomKm(map_mMomVecKm[phiMixKey][iTrkKm]); // K-
	  mPhiMesonTrack->setMass2Kp(map_mMass2Kp[phiMixKey][iTrkKp]); // K+
	  mPhiMesonTrack->setMass2Km(map_mMass2Km[phiMixKey][iTrkKm]); // K-
	  mPhiMesonTrack->setNSigKp(map_mNSigKp[phiMixKey][iTrkKp]); // K+
	  mPhiMesonTrack->setNSigKm(map_mNSigKm[phiMixKey][iTrkKm]); // K-
	  mPhiMesonTrack->setDcaKp(map_mDcaKp[phiMixKey][iTrkKp]); // K+
	  mPhiMesonTrack->setDcaKm(map_mDcaKm[phiMixKey][iTrkKm]); // K-
	  mPhiMesonTrack->setChargeKp(map_mChargeKp[phiMixKey][iTrkKp]); // K+
	  mPhiMesonTrack->setChargeKm(map_mChargeKm[phiMixKey][iTrkKm]); // K-
	  mPhiMesonTrack->setNHitsFitKp(map_mNHitsFitKp[phiMixKey][iTrkKp]); // K+
	  mPhiMesonTrack->setNHitsFitKm(map_mNHitsFitKm[phiMixKey][iTrkKm]); // K-
	  mPhiMesonTrack->setFlagKp(iEvt); // K+
	  mPhiMesonTrack->setFlagKm(iEvt); // K-
	  h_mInvMassPhi[cent9]->Fill(ptPhi,invMassPhi); // Fill histogram with InvMassAB information
	}
      }
    }
  }
  t_mPhiMesonTree->Fill();
}

void StPhiMesonTree::mixPhi(int cent9, int vzBin, int PsiBin) // reconstruct phi meson in the mixed event
{
  for(int iEvtA = 0; iEvtA < mEventCounter[cent9][vzBin][PsiBin]-1; iEvtA++)
  {
    int phiMixKeyA = getPhiMixKey(cent9,vzBin,PsiBin,iEvtA);
    for(int iEvtB = iEvtA+1; iEvtB < mEventCounter[cent9][vzBin][PsiBin]; iEvtB++)
    {
      int phiMixKeyB = getPhiMixKey(cent9,vzBin,PsiBin,iEvtB);
      if(iEvtA == 0 && iEvtB == 1)
      { // set event header with the info of first event in the event buffer
	int evtBin = iEvtA;
	mPhiMesonEvent->clearTrackList();
	mPhiMesonEvent->setRunId(vec_mRunId[cent9][vzBin][PsiBin][evtBin]); // event header
	mPhiMesonEvent->setEvtId(vec_mEvtId[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setRefMult(vec_mRefMult[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setNumTofMatch(vec_mNumTofMatch[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setCentrality9(vec_mCent9[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setCentrality16(vec_mCent16[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setRefWgt(vec_mRefWgt[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setZDCx(vec_mZDCx[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setBBCx(vec_mBBCx[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setVzVpd(vec_mVzVpd[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setPrimVtx(vec_mPrimVtx[cent9][vzBin][PsiBin][evtBin]);

	mPhiMesonEvent->setFlagZdcEp(vec_mFlagZdcEp[cent9][vzBin][PsiBin][evtBin]); // ZDC EP Info
	mPhiMesonEvent->setQ1VecZdcEast(vec_mQ1ZdcShiftEast[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setQ1VecZdcWest(vec_mQ1ZdcShiftWest[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setQ1VecZdcFull(vec_mQ1ZdcShiftFull[cent9][vzBin][PsiBin][evtBin]);

	mPhiMesonEvent->setFlagEpdEp(vec_mFlagEpdEp[cent9][vzBin][PsiBin][evtBin]); // EPD EP Info
	mPhiMesonEvent->setQ1VecEpdEast(vec_mQ1EpdShiftEast[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setQ1VecEpdWest(vec_mQ1EpdShiftWest[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setQ1VecEpdFull(vec_mQ1EpdShiftFull[cent9][vzBin][PsiBin][evtBin]);

	mPhiMesonEvent->setFlagTpcEp(vec_mFlagTpcEp[cent9][vzBin][PsiBin][evtBin]); // TPC EP Info
	mPhiMesonEvent->setQ2VecTpcEast(vec_mQ2TpcReCtrEast[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setQ2VecTpcWest(vec_mQ2TpcReCtrWest[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setQ3VecTpcEast(vec_mQ3TpcReCtrEast[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setQ3VecTpcWest(vec_mQ3TpcReCtrWest[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setNumTrkReCtrEast(vec_mNumTrkReCtrEast[cent9][vzBin][PsiBin][evtBin]);
	mPhiMesonEvent->setNumTrkReCtrWest(vec_mNumTrkReCtrWest[cent9][vzBin][PsiBin][evtBin]);
      }

      // start events mixing
      TLorentzVector lTrkKpEvtA, lTrkKmEvtB;
      // mix K+ candidates from event A with K- candidates from event B
      for(unsigned int iTrkKp = 0; iTrkKp < map_mMomVecKp[phiMixKeyA].size(); ++iTrkKp)
      { // first track loop over K+ candidates from event A
	TVector3 primMomKp = map_mMomVecKp[phiMixKeyA][iTrkKp];
	lTrkKpEvtA.SetXYZM(primMomKp.X(),primMomKp.Y(),primMomKp.Z(), anaUtils::mMassKaon); // K+

	for(unsigned int iTrkKm = 0; iTrkKm < map_mMomVecKm[phiMixKeyB].size(); ++iTrkKm)
	{ // second track loop over K- candidates from event B
	  TVector3 primMomKm = map_mMomVecKm[phiMixKeyB][iTrkKm];
	  lTrkKmEvtB.SetXYZM(primMomKm.X(),primMomKm.Y(),primMomKm.Z(),anaUtils::mMassKaon); // K-

	  TLorentzVector lTrkPhi = lTrkKpEvtA+lTrkKmEvtB; // K+ from EvtA + K- from EvtB
	  double invMassPhi = lTrkPhi.M();
	  double ptPhi      = lTrkPhi.Perp();

	  // fill phi background (mass within [0.98, 1.08]) into t_mPhiMesonTree
	  if(invMassPhi >= anaUtils::mMassPhiMin && invMassPhi <= anaUtils::mMassPhiMax) 
	  {
	    mPhiMesonTrack = mPhiMesonEvent->createTrack();
	    mPhiMesonTrack->setTrkMomKp(map_mMomVecKp[phiMixKeyA][iTrkKp]); // K+ from EvtA
	    mPhiMesonTrack->setTrkMomKm(map_mMomVecKm[phiMixKeyB][iTrkKm]); // K- from EvtB
	    mPhiMesonTrack->setMass2Kp(map_mMass2Kp[phiMixKeyA][iTrkKp]); // K+ from EvtA
	    mPhiMesonTrack->setMass2Km(map_mMass2Km[phiMixKeyB][iTrkKm]); // K- from EvtB
	    mPhiMesonTrack->setNSigKp(map_mNSigKp[phiMixKeyA][iTrkKp]); // K+ from EvtA
	    mPhiMesonTrack->setNSigKm(map_mNSigKm[phiMixKeyB][iTrkKm]); // K- from EvtB
	    mPhiMesonTrack->setDcaKp(map_mDcaKp[phiMixKeyA][iTrkKp]); // K+ from EvtA
	    mPhiMesonTrack->setDcaKm(map_mDcaKm[phiMixKeyB][iTrkKm]); // K- from EvtB
	    mPhiMesonTrack->setChargeKp(map_mChargeKp[phiMixKeyA][iTrkKp]); // K+ from EvtA
	    mPhiMesonTrack->setChargeKm(map_mChargeKm[phiMixKeyB][iTrkKm]); // K- from EvtB
	    mPhiMesonTrack->setNHitsFitKp(map_mNHitsFitKp[phiMixKeyA][iTrkKp]); // K+ from EvtA
	    mPhiMesonTrack->setNHitsFitKm(map_mNHitsFitKm[phiMixKeyB][iTrkKm]); // K- from EvtB
	    mPhiMesonTrack->setFlagKp(iEvtA); // K+ from EvtA
	    mPhiMesonTrack->setFlagKm(iEvtB); // K- from EvtB
	    h_mInvMassPhi[cent9]->Fill(ptPhi,invMassPhi); // Fill histogram with InvMassAB information
	  }
	}
      }

      TLorentzVector lTrkKmEvtA, lTrkKpEvtB;
      // mix K- candidates from event A with K+ candidates from event B
      for(unsigned int iTrkKm = 0; iTrkKm < map_mMomVecKm[phiMixKeyA].size(); ++iTrkKm)
      { // first track loop over K- candidates from event A
	TVector3 primMomKm = map_mMomVecKm[phiMixKeyA][iTrkKm];
	lTrkKmEvtA.SetXYZM(primMomKm.X(),primMomKm.Y(),primMomKm.Z(),anaUtils::mMassKaon); // K-

	for(unsigned int iTrkKp = 0; iTrkKp < map_mMomVecKp[phiMixKeyB].size(); ++iTrkKp)
	{ // second track loop over K+ candidates from event B
	  TVector3 primMomKp = map_mMomVecKp[phiMixKeyB][iTrkKp];
	  lTrkKpEvtB.SetXYZM(primMomKp.X(),primMomKp.Y(),primMomKp.Z(), anaUtils::mMassKaon); // K+

	  TLorentzVector lTrkPhi = lTrkKmEvtA+lTrkKpEvtB; // K- from EvtA + K+ from EvtB
	  double invMassPhi = lTrkPhi.M();
	  double ptPhi      = lTrkPhi.Perp();

	  // fill phi background (mass within [0.98, 1.08]) into t_mPhiMesonTree
	  if(invMassPhi >= anaUtils::mMassPhiMin && invMassPhi <= anaUtils::mMassPhiMax) 
	  {
	    mPhiMesonTrack = mPhiMesonEvent->createTrack();
	    mPhiMesonTrack->setTrkMomKp(map_mMomVecKp[phiMixKeyB][iTrkKp]); // K+ from EvtB
	    mPhiMesonTrack->setTrkMomKm(map_mMomVecKm[phiMixKeyA][iTrkKm]); // K- from EvtA
	    mPhiMesonTrack->setMass2Kp(map_mMass2Kp[phiMixKeyB][iTrkKp]); // K+ from EvtB
	    mPhiMesonTrack->setMass2Km(map_mMass2Km[phiMixKeyA][iTrkKm]); // K- from EvtA
	    mPhiMesonTrack->setNSigKp(map_mNSigKp[phiMixKeyB][iTrkKp]); // K+ from EvtB
	    mPhiMesonTrack->setNSigKm(map_mNSigKm[phiMixKeyA][iTrkKm]); // K- from EvtA
	    mPhiMesonTrack->setDcaKp(map_mDcaKp[phiMixKeyB][iTrkKp]); // K+ from EvtB
	    mPhiMesonTrack->setDcaKm(map_mDcaKm[phiMixKeyA][iTrkKm]); // K- from EvtA
	    mPhiMesonTrack->setChargeKp(map_mChargeKp[phiMixKeyB][iTrkKp]); // K+ from EvtB
	    mPhiMesonTrack->setChargeKm(map_mChargeKm[phiMixKeyA][iTrkKm]); // K- from EvtA
	    mPhiMesonTrack->setNHitsFitKp(map_mNHitsFitKp[phiMixKeyB][iTrkKp]); // K+ from EvtB
	    mPhiMesonTrack->setNHitsFitKm(map_mNHitsFitKm[phiMixKeyA][iTrkKm]); // K- from EvtA
	    mPhiMesonTrack->setFlagKp(iEvtB); // K+ from EvtB
	    mPhiMesonTrack->setFlagKm(iEvtA); // K- from EvtA
	    h_mInvMassPhi[cent9]->Fill(ptPhi,invMassPhi); // Fill histogram with InvMassAB information
	  }
	}
      }
    }
  }
  t_mPhiMesonTree->Fill();
}
//------------------------------------------------------------------------------------------------------------------
int StPhiMesonTree::getPhiMixKey(int cent9, int vzBin, int PsiBin, int evtBin)
{ // return 1000*cent9 + 100*vzBin + 10*PsiBin + evtBin
  int phiMixKey = 1000*cent9 + 100*vzBin + 10*PsiBin + evtBin;

  return phiMixKey;
}

int StPhiMesonTree::getVzMixBin(double vz)
{
  int vzBin = -1;

  double vzBinSize = (anaUtils::mVzMax[mType]-anaUtils::mVzMin[mType])/mNumMixVzBin; // 6cm for IsoBar

  if(std::abs(vz-anaUtils::mVzMin[mType]) < std::numeric_limits<double>::epsilon()) vzBin = 0;
  for(int iVz = 0; iVz < mNumMixVzBin; ++iVz)
  {
    if((vz > anaUtils::mVzMin[mType]+iVz*vzBinSize) && (vz <= anaUtils::mVzMin[mType]+(iVz+1)*vzBinSize))
    {
      vzBin = iVz;
    }
  }

  return vzBin;
}

int StPhiMesonTree::getPsi2MixBin(double Psi2)
{
  int Psi2Bin = -1;
  double Psi2Max = TMath::Pi()/2.0;
  double Psi2Min = -1.0*TMath::Pi()/2.0;
  double Psi2BinSize = (Psi2Max-Psi2Min)/mNumMixPsiBin; // pi/5

  if(std::abs(Psi2-Psi2Min) < std::numeric_limits<double>::epsilon()) Psi2Bin = 0;
  for(int iPsi = 0; iPsi < mNumMixPsiBin; ++iPsi)
  {
    if((Psi2 > Psi2Min+iPsi*Psi2BinSize) && (Psi2 <= Psi2Min+(iPsi+1)*Psi2BinSize))
    {
      Psi2Bin = iPsi;
    }
  }

  return Psi2Bin;
}

void StPhiMesonTree::getPhiEvtSize(int cent9, int vzBin, int PsiBin)
{
  std::cout << "Event Buffer: Centrality = " << cent9 << ", VertexZ = " << vzBin << ", Psi2 = " << PsiBin << std::endl;
  std::cout << "Buffer Depth: " << mEventCounter[cent9][vzBin][PsiBin] << std::endl;

  std::cout << "Size of primaryVertex = " << vec_mPrimVtx[cent9][vzBin][PsiBin].size() << std::endl;;
  std::cout << "Size of refMult       = " << vec_mRefMult[cent9][vzBin][PsiBin].size() << std::endl;;
  std::cout << "---------------------------------------------------------------------------" << std::endl;

  for(int evtBin = 0; evtBin < mEventCounter[cent9][vzBin][PsiBin]; evtBin++)
  {
    int phiMixKey = getPhiMixKey(cent9,vzBin,PsiBin,evtBin);
    std::cout << "Event Number " << evtBin << ":" << std::endl; 
    std::cout << "Positive Particle:" << std::endl;
    std::cout << "  Size of MomVecKp     = " << map_mMomVecKp[phiMixKey].size() << std::endl;;
    std::cout << "  Size of Mass2        = " << map_mMass2Kp[phiMixKey].size() << std::endl;
    std::cout << "  Size of nSigmaKaon   = " << map_mNSigKp[phiMixKey].size() << std::endl;
    std::cout << "  Size of dca          = " << map_mDcaKp[phiMixKey].size() << std::endl;
    std::cout << "  Size of charge       = " << map_mChargeKp[phiMixKey].size() << std::endl;
    std::cout << "  Size of nHitsFit     = " << map_mNHitsFitKp[phiMixKey].size() << std::endl;

    std::cout << "Negative Particle:" << std::endl;
    std::cout << "  Size of MomVecKm     = " << map_mMomVecKm[phiMixKey].size() << std::endl;;
    std::cout << "  Size of Mass2        = " << map_mMass2Km[phiMixKey].size() << std::endl;
    std::cout << "  Size of nSigmaKaon   = " << map_mNSigKm[phiMixKey].size() << std::endl;
    std::cout << "  Size of dca          = " << map_mDcaKm[phiMixKey].size() << std::endl;
    std::cout << "  Size of charge       = " << map_mChargeKm[phiMixKey].size() << std::endl;
    std::cout << "  Size of nHitsFit     = " << map_mNHitsFitKm[phiMixKey].size() << std::endl;
    std::cout << "---------------------------------------------------------------------------" << std::endl;
  }
}
//------------------------------------------------------------------------------------------------------------------
// set QVector from StPhiMesonMaker
void StPhiMesonTree::clearEvtInfo()
{
  mCent9         = -1;
  mCent16        = -1;
  mRefWgt        = -1.0;
  mVz            = -999.9;
  mPsi2ShiftFull = -999.9;

  mFlagZdcEp = 0;
  v_mQ1VecZdcShiftEast.Set(0.0,0.0); 
  v_mQ1VecZdcShiftWest.Set(0.0,0.0); 
  v_mQ1VecZdcShiftFull.Set(0.0,0.0);

  mFlagEpdEp = 0;
  v_mQ1VecEpdShiftEast.Set(0.0,0.0); 
  v_mQ1VecEpdShiftWest.Set(0.0,0.0); 
  v_mQ1VecEpdShiftFull.Set(0.0,0.0);

  mFlagTpcEp = 0;
  v_mQ2VecTpcReCtrEast.Set(0.0,0.0);
  v_mQ2VecTpcReCtrWest.Set(0.0,0.0); 
  v_mQ3VecTpcReCtrEast.Set(0.0,0.0); 
  v_mQ3VecTpcReCtrWest.Set(0.0,0.0);

  mNumTrkReCtrEast = 0; 
  mNumTrkReCtrWest = 0;
}

void StPhiMesonTree::setEvtInfo(int cent9, int cent16, double refwgt, double vz, double Psi2ShiftFull)
{
  mCent9         = cent9;
  mCent16        = cent16;
  mRefWgt        = refwgt;
  mVz            = vz;
  mPsi2ShiftFull = Psi2ShiftFull;
}

void StPhiMesonTree::setZdcQ1Vec(int flagEp, TVector2 Q1East, TVector2 Q1West, TVector2 Q1Full)
{
  mFlagZdcEp = flagEp;
  v_mQ1VecZdcShiftEast = Q1East;
  v_mQ1VecZdcShiftWest = Q1West;
  v_mQ1VecZdcShiftFull = Q1Full;
}

void StPhiMesonTree::setEpdQ1Vec(int flagEp, TVector2 Q1East, TVector2 Q1West, TVector2 Q1Full)
{
  mFlagEpdEp = flagEp;
  v_mQ1VecEpdShiftEast = Q1East;
  v_mQ1VecEpdShiftWest = Q1West;
  v_mQ1VecEpdShiftFull = Q1Full;
}

void StPhiMesonTree::setTpcQVec(int flagEp, TVector2 Q2East, TVector2 Q2West, TVector2 Q3East, TVector2 Q3West)
{
  mFlagTpcEp = flagEp;
  v_mQ2VecTpcReCtrEast = Q2East;
  v_mQ2VecTpcReCtrWest = Q2West;
  v_mQ3VecTpcReCtrEast = Q3East;
  v_mQ3VecTpcReCtrWest = Q3West;
}

void StPhiMesonTree::setNumTrks(int numTrkEast, int numTrkWest)
{
  mNumTrkReCtrEast  = numTrkEast;
  mNumTrkReCtrWest  = numTrkWest;
}
//------------------------------------------------------------------------------------------------------------------
