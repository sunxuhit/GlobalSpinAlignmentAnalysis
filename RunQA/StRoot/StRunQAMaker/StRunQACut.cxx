#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
// #include "StRefMultCorr/StRefMultCorr.h"
// #include "StRefMultCorr/CentralityMaker.h"
#include "StMessMgr.h"

#include "StRoot/StRunQAMaker/StRunQACut.h"
#include "StRoot/StRunQAMaker/StRunQACons.h"

#include "TVector3.h"

ClassImp(StRunQACut)

//---------------------------------------------------------------------------------

StRunQACut::StRunQACut(int energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StRunQACut::~StRunQACut()
{
  /* */
}

//---------------------------------------------------------------------------------

// Event Cuts
bool StRunQACut::isMinBias(StPicoEvent *picoEvent)
{
  // std::cout << "year: " << picoEvent->year() << std::endl;
  // std::cout << "day: " << picoEvent->day() << std::endl;
  // std::cout << "triggerIds: " << picoEvent->triggerIds()[0] << std::endl;
  if(mEnergy == 0 && runQA::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(450005) || picoEvent->isTrigger(450015) || picoEvent->isTrigger(450025) || picoEvent->isTrigger(450050) || picoEvent->isTrigger(450060) )) return false; // 200GeV_2014
  if(mEnergy == 1 && runQA::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(580001) || picoEvent->isTrigger(580021) )) return false; // 54GeV_2017 | 580011 ?
  if(mEnergy == 2 && runQA::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(610001) || picoEvent->isTrigger(610011) || picoEvent->isTrigger(610021) || picoEvent->isTrigger(610031) || picoEvent->isTrigger(610041) || picoEvent->isTrigger(610051) )) return false; // 27GeV_2018

  return true;
}

bool StRunQACut::isBES(int energy)
{
  if(energy == 0) return false; // 200 GeV

  return true; // BES
}

bool StRunQACut::passEventCut(StPicoDst *picoDst)
{
  StPicoEvent *picoEvent = picoDst->event();
  if(!picoEvent)
  {
    return kFALSE;
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
  if(fabs(vz) > runQA::mVzMaxMap[mEnergy])
  {
    return kFALSE;
  }
  // vr cut
  if(sqrt(vx*vx+vy*vy) > runQA::mVrMax[mEnergy])
  {
    return kFALSE;
  }
  // vz-vzVpd cut only for 200 GeV
  if(!isBES(mEnergy) && fabs(vz-vzVpd) > runQA::mVzVpdDiffMax[mEnergy])
  {
    return kFALSE;
  }

  /*
  // ToF Match & ToF Hits cut
  const int refMult = picoEvent->refMult();
  const unsigned short numOfBTofMatch = picoEvent->nBTOFMatch();
  const unsigned int numOfBTofHits = picoDst->numberOfBTofHits(); // get number of tof hits
  if(numOfBTofMatch <= runQA::mMatchedToFMin[mEnergy])
  {
    return kFALSE;
  }

  if(mEnergy == 1) // ToF Hits vs RefMult cut for 54 GeV
  { // from Shaowei Lan
    float tofHits_low = (float)refMult*2.88 - 155.0;
    if( numOfBTofHits < tofHits_low )
    {
      return kFALSE;
    }
  }

  if(mEnergy == 2) // ToF Match vs. RefMult cut for 27 GeV
  { // from Zaochen Ye
    float tofMatch_up = (float)refMult*1.8 + 15.0;
    float tofMatch_low = (float)refMult*0.75 - 20.0;
    if( numOfBTofMatch > tofMatch_up || numOfBTofMatch < tofMatch_low )
    {
      return kFALSE;
    }
  }
  */

  return kTRUE;
}

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
bool StRunQACut::passTrackBasic(StPicoTrack *picoTrack)
{
  if(!picoTrack) return kFALSE;

  // nHitsFit cut
  if(picoTrack->nHitsFit() < runQA::mHitsFitTPCMin)
  {
    return kFALSE;
  }

  // nHitsRatio cut
  if(picoTrack->nHitsMax() <= runQA::mHitsMaxTPCMin)
  {
    return kFALSE;
  }
  if((float)picoTrack->nHitsFit()/(float)picoTrack->nHitsMax() < runQA::mHitsRatioTPCMin)
  {
    return kFALSE;
  }

  // eta cut
  // float eta = picoTrack->pMom().pseudoRapidity();
  // float eta = picoTrack->pMom().PseudoRapidity();
  TVector3 primMom; // temp fix for StThreeVectorF & TVector3
  float primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
  float primPy    = picoTrack->pMom().y();
  float primPz    = picoTrack->pMom().z();
  primMom.SetXYZ(primPx,primPy,primPz);
  float eta = primMom.PseudoRapidity();
  if(fabs(eta) > runQA::mEtaMax)
  {
    return kFALSE;
  }

  if(primMom.Pt() < runQA::mGlobPtMin) // minimum pT cuts
  {
    return kFALSE;
  }

  return kTRUE;
}

bool StRunQACut::passTrackQA(StPicoTrack *picoTrack, StPicoEvent *picoEvent)
{
  if(!passTrackBasic(picoTrack)) return kFALSE;
  if(!picoEvent) return kFALSE;

  const float vx    = picoEvent->primaryVertex().x(); // x works for both TVector3 and StThreeVectorF
  const float vy    = picoEvent->primaryVertex().y();
  const float vz    = picoEvent->primaryVertex().z();

  if(picoTrack->gDCA(vx,vy,vz) > runQA::mDcaTrQAMax)
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------
float StRunQACut::getBeta(StPicoDst *picoDst, int i_track)
{
  float Beta = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    Beta = tofTrack->btofBeta();
  }

  return Beta;
}

float StRunQACut::getPrimaryMass2(StPicoDst *picoDst, int i_track)
{
  float Mass2 = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    float Beta = tofTrack->btofBeta();
    // float Momentum = picoTrack->pMom().mag(); // primary momentum for 54GeV_2017
    // float Momentum = picoTrack->pMom().Mag(); // primary momentum for 27GeV_2018
    TVector3 primMom; // temp fix for StThreeVectorF & TVector3
    float primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
    float primPy    = picoTrack->pMom().y();
    float primPz    = picoTrack->pMom().z();
    primMom.SetXYZ(primPx,primPy,primPz);
    float Momentum = primMom.Mag(); // primary momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && Beta != 0)
    {
      Mass2 = Momentum*Momentum*(1.0/(Beta*Beta) - 1.0);
    }
  }

  return Mass2;
}

float StRunQACut::getGlobalMass2(StPicoDst *picoDst, int i_track)
{
  float Mass2 = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    float Beta = tofTrack->btofBeta();
    // float Momentum = picoTrack->gMom().mag(); // global momentum for 54GeV_2017
    // float Momentum = picoTrack->gMom().Mag(); // global momentum for 27GeV_2018
    TVector3 globMom; // temp fix for StThreeVectorF & TVector3
    float globPx     = picoTrack->gMom().x(); // x works for both TVector3 and StThreeVectorF
    float globPy     = picoTrack->gMom().y();
    float globPz     = picoTrack->gMom().z();
    globMom.SetXYZ(globPx,globPy,globPz);
    float Momentum = globMom.Mag(); // global momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && Beta != 0)
    {
      Mass2 = Momentum*Momentum*(1.0/(Beta*Beta) - 1.0);
    }
  }

  return Mass2;
}

int StRunQACut::getTriggerBin(StPicoEvent *picoEvent)
{
  // std::cout << "year: " << picoEvent->year() << std::endl;
  // std::cout << "day: " << picoEvent->day() << std::endl;
  // std::cout << "triggerIds: " << picoEvent->triggerIds()[0] << std::endl;
  if( mEnergy == 0 && runQA::mBeamYear[mEnergy] == picoEvent->year() )
  { // 200GeV_2014
    if( picoEvent->isTrigger(450005) ) return 0; // VPDMB-5-p-nobsmd
    if( picoEvent->isTrigger(450015) ) return 1; // VPDMB-5-p-nobsmd
    if( picoEvent->isTrigger(450025) ) return 2; // VPDMB-5-p-nobsmd
    if( picoEvent->isTrigger(450050) ) return 3; // VPDMB-5-p-nobsmd-hlt
    if( picoEvent->isTrigger(450060) ) return 4; // VPDMB-5-p-nobsmd-hlt
  }
  if( mEnergy == 1 && runQA::mBeamYear[mEnergy] == picoEvent->year() )
  { // 54GeV_2017
    if( picoEvent->isTrigger(580001) ) return 0; // minBias
    if( picoEvent->isTrigger(580021) ) return 1; // minBias
  }
  if( mEnergy == 2 && runQA::mBeamYear[mEnergy] == picoEvent->year() )
  { // 27GeV_2018
    if( picoEvent->isTrigger(610001) ) return 0; // mb
    if( picoEvent->isTrigger(610011) ) return 1; // mb
    if( picoEvent->isTrigger(610021) ) return 2; // mb
    if( picoEvent->isTrigger(610031) ) return 3; // mb
    if( picoEvent->isTrigger(610041) ) return 4; // mb
    if( picoEvent->isTrigger(610051) ) return 5; // mb
  }

  return -1;
}
//---------------------------------------------------------------------------------
