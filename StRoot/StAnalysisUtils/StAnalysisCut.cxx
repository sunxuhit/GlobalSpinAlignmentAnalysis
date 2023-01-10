#include <iostream>

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
// #include "StMessMgr.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"

#include "TVector3.h"

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

bool StAnalysisCut::passEventCut(StPicoDst *picoDst)
{
  StPicoEvent *picoEvent = picoDst->event();
  if(!picoEvent)
  {
    return false;
  }

  // const int runId   = picoEvent->runId();
  // const int refMult = picoEvent->refMult();
  // cout << "runId = " << runId << ", refMult = " << refMult << endl;

  // event vertex cut
  const float vx    = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const float vy    = picoEvent->primaryVertex().y();
  const float vz    = picoEvent->primaryVertex().z();
  // const float zdcX  = picoEvent->ZDCx();
  const float vzVpd = picoEvent->vzVpd();
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
  if(picoTrack->nHitsFit() < anaUtils::mHitsFitTpcMin)
  {
    return false;
  }

  // nHitsRatio cut
  if(picoTrack->nHitsMax() <= anaUtils::mHitsMaxTpcMin)
  {
    return false;
  }
  if((float)picoTrack->nHitsFit()/(float)picoTrack->nHitsMax() < anaUtils::mHitsRatioTpcMin)
  {
    return false;
  }

  return true;
}

bool StAnalysisCut::passTrackQA(StPicoTrack *picoTrack, StPicoEvent *picoEvent)
{
  if(!picoEvent) return false;
  if(!passTrackBasic(picoTrack)) return false;

  const float vx    = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const float vy    = picoEvent->primaryVertex().y();
  const float vz    = picoEvent->primaryVertex().z();

  // dca cut
  if(picoTrack->gDCA(vx,vy,vz) > anaUtils::mDcaQaMax)
  {
    return false;
  }

  // eta cut
  // TVector3 primMom; // temp fix for StThreeVectorF & TVector3
  // float primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
  // float primPy    = picoTrack->pMom().y();
  // float primPz    = picoTrack->pMom().z();
  // primMom.SetXYZ(primPx,primPy,primPz);
  const TVector3 primMom = picoTrack->pMom(); // primary Momentum
  const float eta = primMom.PseudoRapidity();
  if(fabs(eta) > anaUtils::mEtaQaMax)
  {
    return false;
  }

  if(primMom.Pt() < anaUtils::mPrimPtQaMin) // minimum pT cuts
  {
    return false;
  }

  return true;
}
//---------------------------------------------------------------------------------
