#include <iostream>

#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"

ClassImp(StAnalysisCut)

//---------------------------------------------------------------------------------

StAnalysisCut::StAnalysisCut(int beamType) : mType(beamType)
{
  // mType = beamType;
}

//---------------------------------------------------------------------------------

StAnalysisCut::~StAnalysisCut()
{
  /* */
}

//---------------------------------------------------------------------------------
// Run Cuts
bool StAnalysisCut::isFxt()
{
  if(mType == 0 || mType == 1) return false; // Isobar

  return true; // Fixed Target
}

bool StAnalysisCut::isIsobar()
{
  if(mType == 0 || mType == 1) return true; // Isobar

  return false; // Fixed Target
}
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// Event Cuts
bool StAnalysisCut::isMinBias(StPicoEvent *picoEvent)
{
  // std::cout << "year: " << picoEvent->year() << std::endl;
  // std::cout << "day: " << picoEvent->day() << std::endl;
  // std::cout << "triggerIds: " << picoEvent->triggerIds()[0] << std::endl;
  if( (mType == 0 || mType == 1) && globCons::mBeamYear[mType] == picoEvent->year() && !(picoEvent->isTrigger(600001) || picoEvent->isTrigger(600011) || picoEvent->isTrigger(600021) || picoEvent->isTrigger(600031)) ) return false; // ZrZr200GeV_2018 || RuRu200GeV_2018
  if( mType == 2 && globCons::mBeamYear[mType] == picoEvent->year() && !(picoEvent->isTrigger(620052)) ) return false; // Fxt3p85GeV_2018

  return true;
}

bool StAnalysisCut::isPileUpEvent(double refMult, double numOfBTofMatch, double vz)
{
  if(this->isIsobar()) return false; // use StRefMultCorr for Isobar runs
  if(this->isFxt()) 
  {
    return false; // under development & always return false for now
  }

  return true;
}

bool StAnalysisCut::isGoodCent9(int cent9)
{
  if(cent9 < 0 || cent9 > 8) return false; // not within 0-80%

  return true;
}

bool StAnalysisCut::passEventCut(StPicoEvent *picoEvent)
{
  if(!picoEvent) return false;

  // const int runId   = picoEvent->runId();
  // const int refMult = picoEvent->refMult();
  // cout << "runId = " << runId << ", refMult = " << refMult << endl;

  // event vertex cut
  const double vx    = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const double vy    = picoEvent->primaryVertex().y();
  const double vz    = picoEvent->primaryVertex().z();
  // const double zdcX  = picoEvent->ZDCx();
  const double vzVpd = picoEvent->vzVpd();
  // vz cut
  if(vz < anaUtils::mVzMin[mType] || vz > anaUtils::mVzMax[mType])
  {
    return false;
  }
  // vr cut
  const double vxReCtr = vx - anaUtils::mVxCtr[mType];
  const double vyReCtr = vy - anaUtils::mVyCtr[mType];
  if(sqrt(vxReCtr*vxReCtr+vyReCtr*vyReCtr) > anaUtils::mVrMax[mType])
  {
    return false;
  }
  // vz-vzVpd cut only for ZrZr200GeV_2018 & RuRu200GeV_2018
  if(isIsobar() && fabs(vz-vzVpd) > anaUtils::mVzVpdDiffMax[mType])
  {
    return false;
  }

  // nTofMatch > 2
  const unsigned short numOfBTofMatch = picoEvent->nBTOFMatch(); // get number of tof match points
  if(numOfBTofMatch <= anaUtils::mMatchedToFMin[mType]) return false;

  return true;
}
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// Track Cuts
bool StAnalysisCut::passTrkBasic(StPicoTrack *picoTrack)
{
  if(!picoTrack) return false;

  // nHitsFit cut
  if(picoTrack->nHitsFit() < anaUtils::mHitsFitTpcMin[mType])
  {
    return false;
  }

  // nHitsRatio cut
  if(picoTrack->nHitsMax() <= anaUtils::mHitsMaxTpcMin[mType])
  {
    return false;
  }
  if((double)picoTrack->nHitsFit()/(double)picoTrack->nHitsMax() < anaUtils::mHitsRatioTpcMin[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkQA(StPicoTrack *picoTrack, TVector3 primVtx)
{
  if(!passTrkBasic(picoTrack)) return false;

  const double vx = primVtx.x();
  const double vy = primVtx.y();
  const double vz = primVtx.z();

  // dca cut
  if(picoTrack->gDCA(vx,vy,vz) > anaUtils::mDcaQaMax[mType])
  {
    return false;
  }

  // eta cut: for ZrZr200GeV_2018 & RuRu200GeV_2018 [-1.0,1.0] | for Fxt3p85GeV_2018 [-2.0,0.0]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaQaMin[mType] || eta > anaUtils::mEtaQaMax[mType])
  {
    return false;
  }

  if(primMom.Pt() < anaUtils::mPrimPtQaMin[mType]) // minimum pT cuts
  {
    return false;
  }

  return true;
}
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// Track Cuts for TPC EP
bool StAnalysisCut::passTrkTpcEpFull(StPicoTrack *picoTrack, TVector3 primVtx)
{
  if(!passTrkBasic(picoTrack)) return false;

  const double vx = primVtx.x();
  const double vy = primVtx.y();
  const double vz = primVtx.z();

  // dca cut
  if(picoTrack->gDCA(vx,vy,vz) > anaUtils::mDcaEpMax[mType])
  {
    return false;
  }

  // eta cut: for ZrZr200GeV_2018 & RuRu200GeV_2018 [-1.0,1.0] | for Fxt3p85GeV_2018 [-2.0,0.0]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaEpMin[mType] || eta > anaUtils::mEtaEpMax[mType])
  {
    return false;
  }

  // momentum cut: 0.2 <= pT <= 2.0 && p <= 10.0 GeV/c
  if(primMom.Pt() < anaUtils::mPrimPtEpMin[mType] || primMom.Pt() > anaUtils::mPrimPtEpMax[mType] || primMom.Mag() > anaUtils::mPrimMomEpMax[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTpcEpEast(StPicoTrack *picoTrack, TVector3 primVtx) // neg
{
  if(!passTrkTpcEpFull(picoTrack, primVtx)) return false;

  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaEpMin[mType] || eta > anaUtils::mEtaEpCtr[mType]-anaUtils::mEtaEpGap[mType])
  { // eta cut: [anaUtils::mEtaEpMin[mType], anaUtils::mEtaEpCtr[mType]-anaUtils::mEtaEpGap[mType]]
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTpcEpWest(StPicoTrack *picoTrack, TVector3 primVtx) // pos
{
  if(!passTrkTpcEpFull(picoTrack, primVtx)) return false;

  // eta cut: [anaUtils::mEtaEpCtr[mType]+anaUtils::mEtaEpGap[mType], anaUtils::mEtaEpMax[mType]]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaEpCtr[mType]+anaUtils::mEtaEpGap[mType] || eta > anaUtils::mEtaEpMax[mType])
  { 
    return false;
  }

  return true;
}

bool StAnalysisCut::passNumTrkTpcSubEpRaw(int numTrackEast, int numTrackWest)
{
  if(numTrackEast < anaUtils::mNumTrackEpMin[mType] || numTrackWest < anaUtils::mNumTrackEpMin[mType])
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StAnalysisCut::passNumTrkTpcSubEpReCtr(int numTrackEast, int numTrackWest)
{
  if(numTrackEast < anaUtils::mNumTrackEpMin[mType] || numTrackWest < anaUtils::mNumTrackEpMin[mType])
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------
// Track Cuts for TPC flow
bool StAnalysisCut::passTrkTpcFlowFull(StPicoTrack *picoTrack, TVector3 primVtx) // neg
{
  if(!passTrkBasic(picoTrack)) return false;

  const double vx = primVtx.x();
  const double vy = primVtx.y();
  const double vz = primVtx.z();

  // dca cut
  if(picoTrack->gDCA(vx,vy,vz) > anaUtils::mDcaKaonMax[mType])
  {
    return false;
  }

  // eta cut: [anaUtils::mEtaKaonMin[mType], anaUtils::mEtaKaonMax[mType]]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaKaonMin[mType] || eta > anaUtils::mEtaKaonMax[mType])
  { 
    return false;
  }

  // momentum cut: pT >= 0.2 GeV/c && p <= 10.0 GeV/c
  if(primMom.Pt() < anaUtils::mPrimPtKaonMin[mType] || primMom.Mag() > anaUtils::mPrimMomKaonMax[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTpcFlowEast(StPicoTrack *picoTrack, TVector3 primVtx) // neg
{
  if(!passTrkTpcFlowFull(picoTrack, primVtx)) return false;

  // eta cut: [anaUtils::mEtaKaonMin[mType], anaUtils::mEtaKaonCtr[mType]]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaKaonMin[mType] || eta > anaUtils::mEtaKaonCtr[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTpcFlowWest(StPicoTrack *picoTrack, TVector3 primVtx) // pos
{
  if(!passTrkTpcFlowFull(picoTrack, primVtx)) return false;

  // eta cut: (mEtaKaonCtr[mType], anaUtils::mEtaKaonMax[mType]]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta <= anaUtils::mEtaKaonCtr[mType] || eta > anaUtils::mEtaKaonMax[mType])
  {
    return false;
  }

  return true;
}
//---------------------------------------------------------------------------------
// Track Cuts for Kaon Candidate
bool StAnalysisCut::passTrkKaonFull(StPicoTrack *picoTrack, TVector3 primVtx)
{
  if(!passTrkBasic(picoTrack)) return false;

  const double vx = primVtx.x();
  const double vy = primVtx.y();
  const double vz = primVtx.z();

  // dca cut
  if(picoTrack->gDCA(vx,vy,vz) > anaUtils::mDcaKaonMax[mType])
  {
    return false;
  }

  // eta cut: for ZrZr200GeV_2018 & RuRu200GeV_2018 [-1.0,1.0] | for Fxt3p85GeV_2018 [-2.0,0.0]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaKaonMin[mType] || eta > anaUtils::mEtaKaonMax[mType])
  {
    return false;
  }

  // momentum cut: pT >= 0.2 && p <= 10.0 GeV/c
  if(primMom.Pt() < anaUtils::mPrimPtKaonMin[mType] || primMom.Mag() > anaUtils::mPrimMomKaonMax[mType])
  {
    return false;
  }

  // nSigmaKaon cut: |nSigmaKaon| <= 3.0
  const double nSigKaon = picoTrack->nSigmaKaon();
  if(nSigKaon < anaUtils::mNSigKaonMin[mType] || nSigKaon > anaUtils::mNSigKaonMax[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkKaonEast(StPicoTrack *picoTrack, TVector3 primVtx) // neg
{
  if(!passTrkKaonFull(picoTrack, primVtx)) return false;

  // eta cut: [-anaUtils::mEtaKaonMax[mType], anaUtils::mEtaKaonCtr[mType]]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaKaonMin[mType] || eta > anaUtils::mEtaKaonCtr[mType])
  { 
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkKaonWest(StPicoTrack *picoTrack, TVector3 primVtx) // pos
{
  if(!passTrkKaonFull(picoTrack, primVtx)) return false;

  // eta cut: (anaUtils::mEtaKaonCtr[mType], anaUtils::mEtaKaonMax[mType]]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta <= anaUtils::mEtaKaonCtr[mType] || eta > anaUtils::mEtaKaonMax[mType])
  {
    return false;
  }

  return true;
}
//---------------------------------------------------------------------------------
// Hit Cuts for EPD EP
bool StAnalysisCut::passHitEpdEpFull(StPicoEpdHit *picoEpdHit)
{
  if(!picoEpdHit) return false;

  // EPD threshhold cut: nMip >= 0.3
  if(picoEpdHit->nMIP() < anaUtils::mMipEpdEpMin[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passHitEpdEpEast(StPicoEpdHit *picoEpdHit) // neg
{
  if(!passHitEpdEpFull(picoEpdHit)) return false;

  const int tileId = picoEpdHit->id(); // tileId < 0 for East EPD
  if(tileId > 0) return false;

  return true;
}

bool StAnalysisCut::passHitEpdEpWest(StPicoEpdHit *picoEpdHit) // pos
{
  if(!passHitEpdEpFull(picoEpdHit)) return false;

  const int tileId = picoEpdHit->id(); // tileId > 0 for West EPD
  if(tileId < 0) return false;

  return true;
}

bool StAnalysisCut::passHitEpdFlowEast(StPicoEpdHit *picoEpdHit) // neg
{
  if(!passHitEpdEpFull(picoEpdHit)) return false;

  const int tileId = picoEpdHit->id(); // tileId < 0 for East EPD
  if(tileId > 0) return false;

  return true;
}

bool StAnalysisCut::passHitEpdFlowWest(StPicoEpdHit *picoEpdHit) // pos
{
  if(!passHitEpdEpFull(picoEpdHit)) return false;

  const int tileId = picoEpdHit->id(); // tileId > 0 for West EPD
  if(tileId < 0) return false;

  return true;
}
