#include <iostream>
#include "TMath.h"

#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonEvent.h"
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
bool StAnalysisCut::isIsobar()
{
  if(mType == 0 || mType == 1) return true; // Isobar

  return false; // Fixed Target
}

bool StAnalysisCut::isFxt3p85GeV_2018()
{
  if(mType == 2) return true; // Fxt3p85GeV_2018

  return false; // Isobar
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
  if(this->isFxt3p85GeV_2018()) return false; // use StPileupUtil for Fxt3p85GeV_2018

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
  const double vxReCtr = vx - anaUtils::mVxCtr[mType]; // equal to StAnalysisUtils::getVxReCtr
  const double vyReCtr = vy - anaUtils::mVyCtr[mType]; // equal to StAnalysisUtils::getVyReCtr
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
  if(isIsobar() && numOfBTofMatch <= anaUtils::mMatchedToFMin[mType]) 
  {
    return false;
  }

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
  if(!picoTrack->isPrimary()) return false; // require primary tracks only

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
  if(!picoTrack->isPrimary()) return false; // require primary tracks only

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
  { // eta cut: [anaUtils::mEtaEpCtr[mType]+anaUtils::mEtaEpGap[mType], anaUtils::mEtaEpMax[mType]]
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
bool StAnalysisCut::passTrkTpcFlowFull(StPicoTrack *picoTrack, TVector3 primVtx)
{
  if(!passTrkBasic(picoTrack)) return false;
  if(!picoTrack->isPrimary()) return false; // require primary tracks only

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
// Track Cuts for Kaon Candidate: used in StPhiMesonMaker
bool StAnalysisCut::passTrkTpcKaonFull(StPicoTrack *picoTrack, TVector3 primVtx)
{
  if(!passTrkBasic(picoTrack)) return false;
  if(!picoTrack->isPrimary()) return false; // require primary tracks only

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

bool StAnalysisCut::passTrkTpcKaonEast(StPicoTrack *picoTrack, TVector3 primVtx) // neg
{
  if(!passTrkTpcKaonFull(picoTrack, primVtx)) return false;

  // eta cut: [-anaUtils::mEtaKaonMax[mType], anaUtils::mEtaKaonCtr[mType]]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaKaonMin[mType] || eta > anaUtils::mEtaKaonCtr[mType])
  { 
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTpcKaonWest(StPicoTrack *picoTrack, TVector3 primVtx) // pos
{
  if(!passTrkTpcKaonFull(picoTrack, primVtx)) return false;

  // eta cut: (anaUtils::mEtaKaonCtr[mType], anaUtils::mEtaKaonMax[mType]]
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta <= anaUtils::mEtaKaonCtr[mType] || eta > anaUtils::mEtaKaonMax[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTofKaonMass(TVector3 primMom, int charge, double mass2)
{
  // if(primMom.Mag() < 0.65 && ((mass2 > anaUtils::mMass2KaonTreeMin[mType] && mass2 < anaUtils::mMass2KaonTreeMax[mType]) || mass2 < -10.0))
  if(primMom.Mag() <= 0.5 && ((mass2 > anaUtils::mMass2KaonTreeMin[mType] && mass2 < anaUtils::mMass2KaonTreeMax[mType]) || mass2 < -10.0))
  { // ToF when valid
    return true;
  }
  // if(primMom.Mag() >= 0.65 && (mass2 > anaUtils::mMass2KaonTreeMin[mType] && mass2 < anaUtils::mMass2KaonTreeMax[mType]))
  if(primMom.Mag() > 0.5 && (mass2 > anaUtils::mMass2KaonTreeMin[mType] && mass2 < anaUtils::mMass2KaonTreeMax[mType]))
  { // ToF always
    return true;
  }

  return false;
}

bool StAnalysisCut::passTrkTofKaonBeta(TVector3 primMom, int charge, double beta)
{
  // Guannan's cuts
  const double eExp = TMath::Sqrt(primMom.Mag2()+anaUtils::mMassKaon*anaUtils::mMassKaon);
  const double betaExp = primMom.Mag()/eExp;
  const double delBeta = TMath::Abs(1.0/beta - 1.0/betaExp);
  if(charge > 0 && primMom.Mag() <= 0.5 && (delBeta < 0.05 || beta < -10.0))
  { // K+: require ToF when valid when mom <= 0.5 GeV/c
    return true;
  }
  if(charge > 0 && primMom.Mag() > 0.5 && delBeta < 0.05)
  { // K+: always return ToF when mom > 0.5 GeV/c
    return true;
  }
  if(charge < 0 && delBeta < 0.05)
  { // K-: always require ToF
    return true;
  }

  return false;
}
//---------------------------------------------------------------------------------
// Track Cuts for Kaon Candidate: used in StPhiMesonAnalyzer
bool StAnalysisCut::passTrkTpcKaonFull(StPhiMesonTrack *phiTrk)
{ // apply to K+&K- track pairs

  // nHitsFit cut
  if(phiTrk->getNHitsFitKp() < 20 ||
     phiTrk->getNHitsFitKm() < 20 ) // temp cuts for cross check with Guannan Xie
  {
    return false;
  }

  // dca cut
  if( phiTrk->getDcaKp() > anaUtils::mDcaKaonMax[mType] ||
      phiTrk->getDcaKm() > anaUtils::mDcaKaonMax[mType] )
  {
    return false;
  }

  // eta cut: for ZrZr200GeV_2018 & RuRu200GeV_2018 [-1.0,1.0] | for Fxt3p85GeV_2018 [-2.0,0.0]
  const TVector3 primMomKp = phiTrk->getTrkMomKp(); // K+ primary Momentum
  const TVector3 primMomKm = phiTrk->getTrkMomKm(); // K- primary Momentum
  const double etaKp = primMomKp.PseudoRapidity();
  const double etaKm = primMomKm.PseudoRapidity();
  if( etaKp < anaUtils::mEtaKaonMin[mType] || etaKp > anaUtils::mEtaKaonMax[mType] ||
      etaKm < anaUtils::mEtaKaonMin[mType] || etaKm > anaUtils::mEtaKaonMax[mType] )
  {
    return false;
  }

  // momentum cut: pT >= 0.2 && p <= 10.0 GeV/c
  if( primMomKp.Pt() < anaUtils::mPrimPtKaonMin[mType] || primMomKp.Mag() > anaUtils::mPrimMomKaonMax[mType] ||
      primMomKm.Pt() < anaUtils::mPrimPtKaonMin[mType] || primMomKm.Mag() > anaUtils::mPrimMomKaonMax[mType] )
  {
    return false;
  }

  // nSigmaKaon cut: |nSigmaKaon| <= 3.0
  const double nSigKp = phiTrk->getNSigKp();
  const double nSigKm = phiTrk->getNSigKm();
  // if( nSigKp < anaUtils::mNSigKaonMin[mType] || nSigKp > anaUtils::mNSigKaonMax[mType] ||
  //     nSigKm < anaUtils::mNSigKaonMin[mType] || nSigKm > anaUtils::mNSigKaonMax[mType] )
  if( nSigKp < -2.5 || nSigKp > 2.5 ||
      nSigKm < -2.5 || nSigKm > 2.5 ) // temp cuts for cross check with Guannan Xie
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTpcKaonEast(StPhiMesonTrack *phiTrk, int charge) // neg
{
  if(!passTrkTpcKaonFull(phiTrk)) return false;

  // eta cut: [-anaUtils::mEtaKaonMax[mType], anaUtils::mEtaKaonCtr[mType]]
  TVector3 primMom;
  if(charge > 0) primMom = phiTrk->getTrkMomKp(); // K+ primary Momentum
  if(charge < 0) primMom = phiTrk->getTrkMomKm(); // K- primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaKaonMin[mType] || eta > anaUtils::mEtaKaonCtr[mType])
  { 
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTpcKaonWest(StPhiMesonTrack *phiTrk, int charge) // pos
{
  if(!passTrkTpcKaonFull(phiTrk)) return false;

  // eta cut: (anaUtils::mEtaKaonCtr[mType], anaUtils::mEtaKaonMax[mType]]
  TVector3 primMom;
  if(charge > 0) primMom = phiTrk->getTrkMomKp(); // K+ primary Momentum
  if(charge < 0) primMom = phiTrk->getTrkMomKm(); // K- primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta <= anaUtils::mEtaKaonCtr[mType] || eta > anaUtils::mEtaKaonMax[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTofKaonMass(StPhiMesonTrack *phiTrk)
{
  const double mass2Kp = phiTrk->getMass2Kp();
  const double mass2Km = phiTrk->getMass2Km();

  if( (mass2Kp > anaUtils::mMass2KaonSpinMin[mType] && mass2Kp < anaUtils::mMass2KaonSpinMax[mType]) &&
      (mass2Km > anaUtils::mMass2KaonSpinMin[mType] && mass2Km < anaUtils::mMass2KaonSpinMax[mType]) )
  { //  always return ToF for K+K- pairs
    return true;
  }

  return false;
}

bool StAnalysisCut::passTrkTofKaonBeta(StPhiMesonTrack *phiTrk)
{
  // Guannan's cuts
  const double betaKp      = phiTrk->getBetaKp(); // K+
  const TVector3 primMomKp = phiTrk->getTrkMomKp();
  const double eKpExp      = TMath::Sqrt(primMomKp.Mag2()+anaUtils::mMassKaon*anaUtils::mMassKaon);
  const double betaKpExp   = primMomKp.Mag()/eKpExp;
  const double delBetaKp   = TMath::Abs(1.0/betaKp - 1.0/betaKpExp);

  const double betaKm      = phiTrk->getBetaKm(); // K-
  const TVector3 primMomKm = phiTrk->getTrkMomKm();
  const double eKmExp      = TMath::Sqrt(primMomKm.Mag2()+anaUtils::mMassKaon*anaUtils::mMassKaon);
  const double betaKmExp   = primMomKm.Mag()/eKmExp;
  const double delBetaKm   = TMath::Abs(1.0/betaKm - 1.0/betaKmExp);

  if(primMomKp.Mag() <= 0.5 && (delBetaKp < 0.03 || betaKp < -10.0) && delBetaKm < 0.03)
  { // K+: require ToF when valid & K-: always require ToF
    return true;
  }
  if(primMomKp.Mag() > 0.5 && delBetaKp < 0.03 && delBetaKm < 0.03 )
  { // always require ToF for K+ & K-
    return true;
  }

  return false;
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

bool StAnalysisCut::passQVecEpdSide(TVector2 Q1VecEast, TVector2 Q1VecWest, TVector2 Q1VecFull)
{
  if(isIsobar() && Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() > 0.0 && Q1VecFull.Mod() > 0.0)
  {
    return true;
  }
  if(isFxt3p85GeV_2018() && Q1VecEast.Mod() > 0.0)
  { // only require East EPD for FXT
    return true;
  }

  return false;
}

bool StAnalysisCut::passQVecEpdGrp(TVector2 Q1VecEast, TVector2 Q1VecWest, TVector2 Q1VecFull, int grpId)
{
  if(isIsobar() && Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() > 0.0 && Q1VecFull.Mod() > 0.0)
  {
    return true;
  }
  if(isFxt3p85GeV_2018() && grpId == 0 && Q1VecEast.Mod() > 0.0)
  { // only require East EPD for FXT
    return true;
  }
  if(isFxt3p85GeV_2018() && grpId == 1 && Q1VecEast.Mod() > 0.0)
  { // only require East EPD for FXT
    return true;
  }

  return false;
}
//---------------------------------------------------------------------------------
// Hit Cuts for ZDC EP
bool StAnalysisCut::passQVecZdc(TVector2 Q1VecEast, TVector2 Q1VecWest, TVector2 Q1VecFull)
{
  if(isIsobar() && Q1VecEast.Mod() > 0.0 && Q1VecWest.Mod() > 0.0 && Q1VecFull.Mod() > 0.0)
  {
    return true;
  }
  if(isFxt3p85GeV_2018() && Q1VecEast.Mod() > 0.0)
  { // only require East EPD for FXT && NOT really used in FXT
    return true;
  }

  return false;
}

// only used for deuteron flow comparison in Fxt3p85GeV_2018
bool StAnalysisCut::passTrkTpcFlow(StPicoTrack *picoTrack, TVector3 primVtx)
{
  if(!picoTrack) return false;
  if(!picoTrack->isPrimary()) return false; // require primary tracks only

  // nHitsFit cut
  if(picoTrack->nHitsFit() <= 15) return false;

  // nHitsRatio cut
  if(picoTrack->nHitsMax() <= anaUtils::mHitsMaxTpcMin[mType]) return false;
  if((double)picoTrack->nHitsFit()/(double)picoTrack->nHitsMax() <= 0.52) return false;

  // dca cut
  const double vx = primVtx.x();
  const double vy = primVtx.y();
  const double vz = primVtx.z();
  if(picoTrack->gDCA(vx,vy,vz) >= 3) return false;

  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(eta < -2.0 || eta > 0.0) return false;

  return true;
}

bool StAnalysisCut::passTrkProFlow(StPicoTrack *picoTrack)
{
  const double nSigProton = picoTrack->nSigmaProton();
  if(nSigProton < -2.0 || nSigProton > 2.0) return false;

  return true;
}

bool StAnalysisCut::passTrkDeuFlow(double pMag, double deuteronZ, double mass2)
{
  double pBins[51] = {
    0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4,
    1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4,
    2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4,
    3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4,
    4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 200};

  int theBin = -1;

  for(int no=0; no<50; no++) {
    if(pMag>=pBins[no] && pMag<pBins[no+1]) {
      theBin = no;}
  }

  if(theBin == -1) return false;

  double zMeans[50] = {
    -0.0832867, -0.0557943, -0.0302878, -0.0154656, 0.00222155, 0.0147798, 0.0243567, 0.0327773, 0.0387324, 0.0372758  ,
    0.0461688 , 0.061474  , 0.057022  , 0.0480333 , 0.0454707 , 0.0466773, 0.0473922, 0.0470951, 0.0481996, 0.0468892  ,
    0.0464752 , 0.0473581 , 0.0467947 , 0.0465219 , 0.0474517 , 0.0474816, 0.0454535, 0.0349253, 0.0282318, 0.0248662  ,
    0.0246939 , 0.0265303 , 0.0284592 , 0.0294229 , 0.0296879 , 0.0322723, 0.0364428, 0.0406073, 0.0451247, 0.0497799  ,
    0.0550221 , 0.061172  , 0.0669082 , 0.0711066 , 0.0890019 , 0.0902058, 0.0907085, 0.0915537, 0.0914749, 0.0874104};

  double dsigma = deuteronZ - zMeans[theBin];

  /*
  // default z cuts
  double lowZ[30] = {
    -0.27 , -0.27 , -0.27 , -0.27 , -0.27 , -0.27 , -0.27 , -0.27 , -0.27 , -0.24,
    -0.22, -0.20, -0.18, -0.16, -0.13, -0.11, -0.09, -0.07, -0.05, -0.04  ,
    -0.04, -0.03, -0.02, -0.01, -0.02, -0.04, -0.06, -0.08, -0.1 , -0.12};

  double highZ[30] = {
    0.27 , 0.27 , 0.27 , 0.27 , 0.27 , 0.27 , 0.27 , 0.27 , 0.27 , 0.27,
    0.27, 0.27, 0.27, 0.26, 0.25, 0.24, 0.22, 0.21, 0.19, 0.17,
    0.14, 0.13, 0.12, 0.11, 0.13, 0.13, 0.12, 0.11, 0.09, 0.08};

  if(pMag<0.8 && dsigma>=-0.3 && dsigma<0.3)
    return true;
  else if(pMag>=0.8 && pMag<3.2 && dsigma>lowZ[theBin-3] && dsigma<highZ[theBin-3])
    return true;
  else if(pMag>=3.2 && dsigma>=-0.4 && dsigma<0.4 && mass2<=4.8 && mass2>=2.8)
    return true;
    */

  // z cuts
  Float_t lowZ[30] = {
    -0.3 , -0.3 , -0.3 , -0.3 , -0.3 , -0.3 , -0.3 , -0.3 , -0.3 , -0.27  ,
    -0.24, -0.21, -0.18, -0.16, -0.13, -0.11, -0.09, -0.07, -0.05, -0.04  ,
    -0.04, -0.03, -0.02, -0.01, -0.02, -0.04, -0.06, -0.08, -0.1 , -0.12};
  Float_t highZ[30] = {
    0.3 , 0.3 , 0.3 , 0.3 , 0.3 , 0.3 , 0.3 , 0.3 , 0.3 , 0.29,
    0.29, 0.27, 0.27, 0.26, 0.25, 0.24, 0.22, 0.21, 0.19, 0.17,
    0.13, 0.14, 0.12, 0.11, 0.13, 0.13, 0.12, 0.11, 0.09, 0.08};

  if(pMag<0.8 && dsigma>=-0.3 && dsigma<0.3)
    return true;
  else if(pMag>=0.8 && pMag<3.2 && dsigma>lowZ[theBin-3] && dsigma<highZ[theBin-3])
    return true;
  else if(pMag>=3.2 && dsigma>=-0.4 && dsigma<0.4 && mass2<=4.8 && mass2>=2.8)
    return true;

  return false;
}
//---------------------------------------------------------------------------------
// Track Cuts for TPC EP: used in StPhiMesonAnalyzer
bool StAnalysisCut::passTrkTpcEpFull(TVector3 primMom, double gDca)
{
  // dca cut
  if(gDca > anaUtils::mDcaEpMax[mType])
  {
    return false;
  }

  // eta cut: for ZrZr200GeV_2018 & RuRu200GeV_2018 [-1.0,1.0] | for Fxt3p85GeV_2018 [-2.0,0.0]
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

bool StAnalysisCut::passTrkTpcEpEast(TVector3 primMom, double gDca) // neg
{
  if(!passTrkTpcEpFull(primMom,gDca)) return false;

  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaEpMin[mType] || eta > anaUtils::mEtaEpCtr[mType]-anaUtils::mEtaEpGap[mType])
  { // eta cut: [anaUtils::mEtaEpMin[mType], anaUtils::mEtaEpCtr[mType]-anaUtils::mEtaEpGap[mType]]
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkTpcEpWest(TVector3 primMom, double gDca) // pos
{
  if(!passTrkTpcEpFull(primMom,gDca)) return false;

  // eta cut: [anaUtils::mEtaEpCtr[mType]+anaUtils::mEtaEpGap[mType], anaUtils::mEtaEpMax[mType]]
  const double eta = primMom.PseudoRapidity();
  if(eta < anaUtils::mEtaEpCtr[mType]+anaUtils::mEtaEpGap[mType] || eta > anaUtils::mEtaEpMax[mType])
  { // eta cut: [anaUtils::mEtaEpCtr[mType]+anaUtils::mEtaEpGap[mType], anaUtils::mEtaEpMax[mType]]
    return false;
  }

  return true;
}
//---------------------------------------------------------------------------------
// Track Cuts for phi flow: used in StPhiMesonAnalyzer
bool StAnalysisCut::passTrkPhiFlowEast(double yPhiCms) // neg
{
  // rapidity cut: [anaUtils::mRapPhiMin[mType], anaUtils::mRapPhiCtr[mType]]
  if(yPhiCms < anaUtils::mRapPhiMin[mType] || yPhiCms > anaUtils::mRapPhiCtr[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrkPhiFlowWest(double yPhiCms) // pos
{
  // rapidity cut: (anaUtils::mRapPhiCtr[mType], anaUtils::mRapPhiMax[mType]]
  if(yPhiCms <= anaUtils::mRapPhiCtr[mType] || yPhiCms > anaUtils::mRapPhiMax[mType])
  {
    return false;
  }

  return true;
}
