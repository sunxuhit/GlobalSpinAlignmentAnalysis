#include <iostream>

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
// #include "StMessMgr.h"

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
bool StAnalysisCut::isFixedTarget()
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
  // if( mType == 2 && globCons::mBeamYear[mType] == picoEvent->year() && !(picoEvent->isTrigger(600001) || picoEvent->isTrigger(600011) || picoEvent->isTrigger(600021) || picoEvent->isTrigger(600031)) ) return false; // Fixed Target

  return true;
}

bool StAnalysisCut::isPileUpEvent(double refMult, double numOfBTofMatch, double vz)
{
  if(this->isIsobar()) return false; // use StRefMultCorr for Isobar runs

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
  if(sqrt(vx*vx+vy*vy) > anaUtils::mVrMax[mType])
  {
    return false;
  }
  // vz-vzVpd cut only for 200 GeV
  if(!isFixedTarget() && fabs(vz-vzVpd) > anaUtils::mVzVpdDiffMax[mType])
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
bool StAnalysisCut::passTrackBasic(StPicoTrack *picoTrack)
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

bool StAnalysisCut::passTrackQA(StPicoTrack *picoTrack, TVector3 primVtx)
{
  if(!passTrackBasic(picoTrack)) return false;

  const double vx = primVtx.x();
  const double vy = primVtx.y();
  const double vz = primVtx.z();

  // dca cut
  if(picoTrack->gDCA(vx,vy,vz) > anaUtils::mDcaQaMax[mType])
  {
    return false;
  }

  // eta cut
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(fabs(eta) > anaUtils::mEtaQaMax[mType])
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
bool StAnalysisCut::passTrackTpcEpFull(StPicoTrack *picoTrack, TVector3 primVtx)
{
  if(!passTrackBasic(picoTrack)) return false;

  const double vx = primVtx.x();
  const double vy = primVtx.y();
  const double vz = primVtx.z();

  // dca cut
  if(picoTrack->gDCA(vx,vy,vz) > anaUtils::mDcaEpMax[mType])
  {
    return false;
  }

  // eta cut
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();
  if(fabs(eta) > anaUtils::mEtaEpMax[mType])
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

bool StAnalysisCut::passTrackTpcEpEast(StPicoTrack *picoTrack, TVector3 primVtx) // neg
{
  if(!passTrackTpcEpFull(picoTrack, primVtx)) return false;

  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();

  // eta cut: [-anaUtils::mEtaEpMax[mType], -anaUtils::mEtaEpGap[mType]]
  if(eta < -1.0*anaUtils::mEtaEpMax[mType] || eta > -1.0*anaUtils::mEtaEpGap[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrackTpcEpWest(StPicoTrack *picoTrack, TVector3 primVtx) // pos
{
  if(!passTrackTpcEpFull(picoTrack, primVtx)) return false;

  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const double eta = primMom.PseudoRapidity();

  // eta cut: [anaUtils::mEtaEpGap[mType], anaUtils::mEtaEpMax[mType]]
  if(eta < anaUtils::mEtaEpGap[mType] || eta > anaUtils::mEtaEpMax[mType])
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passNumTrackTpcSubEpRaw(int numTrackEast, int numTrackWest)
{
  if(numTrackEast < anaUtils::mNumTrackEpMin[mType] || numTrackWest < anaUtils::mNumTrackEpMin[mType])
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StAnalysisCut::passNumTrackTpcSubEpReCenter(int numTrackEast, int numTrackWest)
{
  if(numTrackEast < anaUtils::mNumTrackEpMin[mType] || numTrackWest < anaUtils::mNumTrackEpMin[mType])
  {
    return kFALSE;
  }

  return kTRUE;
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
