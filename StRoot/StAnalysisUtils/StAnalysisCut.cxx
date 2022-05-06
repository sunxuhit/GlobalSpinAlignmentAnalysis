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

StAnalysisCut::StAnalysisCut(int beamType)
{
  mType = beamType;
}

//---------------------------------------------------------------------------------

StAnalysisCut::~StAnalysisCut()
{
  /* */
}

//---------------------------------------------------------------------------------

// Event Cuts
bool StAnalysisCut::isBES()
{
  if(mType == 0 || mType == 1) return false; // Isobar

  return true; // BES
}

bool StAnalysisCut::isIsobar()
{
  if(mType == 0 || mType == 1) return true; // Isobar

  return false; // BES
}

bool StAnalysisCut::isMinBias(StPicoEvent *picoEvent)
{
  // std::cout << "year: " << picoEvent->year() << std::endl;
  // std::cout << "day: " << picoEvent->day() << std::endl;
  // std::cout << "triggerIds: " << picoEvent->triggerIds()[0] << std::endl;
  // if(mType == 0 && globCons::mBeamYear[mType] == picoEvent->year() && !( picoEvent->isTrigger(450005) || picoEvent->isTrigger(450015) || picoEvent->isTrigger(450025) || picoEvent->isTrigger(450050) || picoEvent->isTrigger(450060) )) return false; // 200GeV_2014

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
  if(!isBES() && fabs(vz-vzVpd) > anaUtils::mVzVpdDiffMax[mType])
  {
    return false;
  }

  // nTofMatch > 2
  const unsigned short numOfBTofMatch = picoEvent->nBTOFMatch(); // get number of tof match points
  if(numOfBTofMatch <= anaUtils::mMatchedToFMin) return false;

  return true;
}

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
bool StAnalysisCut::passTrackBasic(StPicoTrack *picoTrack)
{
  if(!picoTrack) return false;

  // nHitsFit cut
  if(picoTrack->nHitsFit() < anaUtils::mHitsFitTPCMin)
  {
    return false;
  }

  // nHitsRatio cut
  if(picoTrack->nHitsMax() <= anaUtils::mHitsMaxTPCMin)
  {
    return false;
  }
  if((float)picoTrack->nHitsFit()/(float)picoTrack->nHitsMax() < anaUtils::mHitsRatioTPCMin)
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
  if(picoTrack->gDCA(vx,vy,vz) > anaUtils::mDcaTrQAMax)
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
  if(fabs(eta) > anaUtils::mEtaMax)
  {
    return false;
  }

  if(primMom.Pt() < anaUtils::mGlobPtMin) // minimum pT cuts
  {
    return false;
  }

  return true;
}

//---------------------------------------------------------------------------------
