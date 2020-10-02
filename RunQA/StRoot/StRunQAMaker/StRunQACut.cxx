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

bool StRunQACut::isBES()
{
  if(mEnergy == 0) return false; // 200 GeV

  return true; // BES
}

bool StRunQACut::isPileUpEvent(int refMult, int numOfBTofMatch, int numOfBTofHits)
{
  if(mEnergy == 0)
  {
    // numOfBTofHits Cuts
    const double h0Lower = -10.6691; // Lower Band
    const double h1Lower = 0.131997;
    const double h2Lower = 0.000108037;
    const double h3Lower = -8.93329e-08;
    const double h4Lower = 3.23579e-11;
    const double h5Lower = -4.66919e-15;
    const double h0LowerExt = 210.652; // pol1 when numOfBTofHits > 2650.0
    const double h1LowerExt = 0.0789303;
    const double h0Upper = 9.03081; // Upper Band
    const double h1Upper = 0.387506;
    const double h2Upper = -0.000141825;
    const double h3Upper = 9.97792e-08;
    const double h4Upper = -3.97333e-11;
    const double h5Upper = 5.42501e-15;
    const double h0UpperExt = 361.012; // pol1 when numOfBTofHits > 2650.0
    const double h1UpperExt = 0.107772;

    double refmultTofHitsLower = h0Lower+h1Lower*(numOfBTofHits)+h2Lower*pow(numOfBTofHits,2)+h3Lower*pow(numOfBTofHits,3)+h4Lower*pow(numOfBTofHits,4)+h5Lower*pow(numOfBTofHits,5);
    double refmultTofHitsUpper = h0Upper+h1Upper*(numOfBTofHits)+h2Upper*pow(numOfBTofHits,2)+h3Upper*pow(numOfBTofHits,3)+h4Upper*pow(numOfBTofHits,4)+h5Upper*pow(numOfBTofHits,5);
    if(numOfBTofHits > 2650)
    {
      refmultTofHitsLower = h0LowerExt+h1LowerExt*(numOfBTofHits);
      refmultTofHitsUpper = h0UpperExt+h1UpperExt*(numOfBTofHits);
    }

    // numOfBTofMatch Cuts
    const double m0Lower = -6.35447;
    const double m1Lower = 0.471171;
    const double m2Lower = 0.00242705;
    const double m3Lower = -7.70005e-06;
    const double m4Lower = 8.47891e-09;
    const double m5Lower = 1.18207e-12;
    const double m0LowerExt = -175.356;
    const double m1LowerExt = 1.20386;
    const double m0Upper = 10.7104;
    const double m1Upper = 1.2943;
    const double m2Upper = 0.00132705;
    const double m3Upper = -1.87584e-05;
    const double m4Upper = 7.54871e-08;
    const double m5Upper = -8.80916e-11;
    const double m0UpperExt = 378.139;
    const double m1UpperExt = 0.518849;

    double refmultTofMatchLower = m0Lower+m1Lower*(numOfBTofMatch)+m2Lower*pow(numOfBTofMatch,2)+m3Lower*pow(numOfBTofMatch,3)+m4Lower*pow(numOfBTofMatch,4)+m5Lower*pow(numOfBTofMatch,5);
    double refmultTofMatchUpper = m0Upper+m1Upper*(numOfBTofMatch)+m2Upper*pow(numOfBTofMatch,2)+m3Upper*pow(numOfBTofMatch,3)+m4Upper*pow(numOfBTofMatch,4)+m5Upper*pow(numOfBTofMatch,5);
    if(numOfBTofMatch > 420)
    {
      refmultTofMatchLower = m0LowerExt+m1LowerExt*(numOfBTofMatch);
      refmultTofMatchUpper = m0UpperExt+m1UpperExt*(numOfBTofMatch);
    }
    if(numOfBTofMatch > 800)
    {
      refmultTofMatchLower = m0LowerExt+m1LowerExt*(numOfBTofMatch);
      refmultTofMatchUpper = refmultTofMatchLower + 50; // only consider the lower limit when numOfBTofMatch > 800 || just for completeness
    }

    // good events: numOfBTofMatch > 2 && refMult within numOfBTofHits Cuts && refMult within numOfBTofMatch Cuts
    if( (numOfBTofMatch > runQA::mMatchedToFMin[mEnergy]) && (refMult >= refmultTofHitsLower && refMult <= refmultTofHitsUpper) && (refMult >= refmultTofMatchLower && refMult <= refmultTofMatchUpper) ) return kFALSE;
  }

  if(mEnergy == 1) // ToF Hits vs RefMult cut for 54 GeV
  { // from Shaowei Lan
    double tofHits_low = (double)refMult*2.88 - 155.0;

    // good events: numOfBTofMatch > 2 && refMult within numOfBTofHits Cuts
    if( (numOfBTofMatch > runQA::mMatchedToFMin[mEnergy]) && (numOfBTofHits >= tofHits_low) ) return kFALSE;
  }

  if(mEnergy == 2)
  { // will use StRefMultCorr
    return kFALSE;
  }

  return kTRUE;
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
  if(!isBES() && fabs(vz-vzVpd) > runQA::mVzVpdDiffMax[mEnergy])
  {
    return kFALSE;
  }

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
