#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
// #include "StRefMultCorr/StRefMultCorr.h"
// #include "StRefMultCorr/CentralityMaker.h"
#include "StMessMgr.h"

#include "StRoot/StVecMesonMaker/StVecMesonCut.h"
#include "StRoot/StVecMesonMaker/StVecMesonCons.h"

ClassImp(StVecMesonCut)

//---------------------------------------------------------------------------------

StVecMesonCut::StVecMesonCut(int energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StVecMesonCut::~StVecMesonCut()
{
  /* */
}

//---------------------------------------------------------------------------------

// Event Cuts
bool StVecMesonCut::isMinBias(StPicoEvent *picoEvent)
{
  // std::cout << "year: " << picoEvent->year() << std::endl;
  // std::cout << "day: " << picoEvent->day() << std::endl;
  // std::cout << "triggerIds: " << picoEvent->triggerIds()[0] << std::endl;
  if(mEnergy == 1 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(580001) || picoEvent->isTrigger(580011) || picoEvent->isTrigger(580021) )) return false; // 54GeV_2017
  if(mEnergy == 2 && vmsa::mBeamYear[mEnergy] == picoEvent->year() && !( picoEvent->isTrigger(610001) || picoEvent->isTrigger(610011) || picoEvent->isTrigger(610021) || picoEvent->isTrigger(610031) || picoEvent->isTrigger(610041) || picoEvent->isTrigger(610051) )) return false; // 27GeV_2018

  return true;
}

bool StVecMesonCut::isBES(int energy)
{
  if(energy == 0) return false; // 200 GeV

  return true; // BES
}

bool StVecMesonCut::passEventCut(StPicoDst *picoDst)
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
  if(fabs(vz) > vmsa::mVzMaxMap[mEnergy])
  {
    return kFALSE;
  }
  // vr cut
  if(sqrt(vx*vx+vy*vy) > vmsa::mVrMax)
  {
    return kFALSE;
  }
  // vz-vzVpd cut
  if(fabs(vz-vzVpd) > vmsa::mVzVpdDiffMax)
  {
    return kFALSE;
  }

  // initialize mMatchedToF
  mMatchedToF = 0;
  mN_prim = 0;
  mN_non_prim = 0;

  // ToF matched points cut
  int nMatchedToF = 0;
  int nN_prim = 0;
  int nN_non_prim = 0;
  const int nTracks = picoDst->numberOfTracks();
  for(int i_track = 0; i_track < nTracks; ++i_track)
  {
    StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track);
    if(!picoTrack)
    {
      continue;
    }
    if(picoTrack->gDCA(vx,vy,vz) > 3) // global track
    {
      nN_non_prim++;
    }
    else
    {
      nN_prim++;
      int tofIndex = picoTrack->bTofPidTraitsIndex();
      if(tofIndex >= 0)
      {
	StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
	if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && tofTrack->btofBeta() != 0)
	{
	  nMatchedToF++;
	}
      }
    }
  }

  mMatchedToF = nMatchedToF;
  mN_prim = nN_prim;
  mN_non_prim = nN_non_prim;

  const unsigned short numOfBTofMatch = picoEvent->nBTOFMatch();
  if(numOfBTofMatch < vmsa::mMatchedToFMin)
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

int StVecMesonCut::getMatchedToF()
{
  return mMatchedToF;
}

int StVecMesonCut::getNpirm()
{
  return mN_prim;
}

int StVecMesonCut::getNnonprim()
{
  return mN_non_prim;
}
//---------------------------------------------------------------------------------
float StVecMesonCut::getBeta(StPicoDst *picoDst, int i_track)
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

float StVecMesonCut::getPrimaryMass2(StPicoDst *picoDst, int i_track)
{
  float Mass2 = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    float Beta = tofTrack->btofBeta();
    float Momentum = picoTrack->pMom().Mag(); // primary momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && Beta != 0)
    {
      Mass2 = Momentum*Momentum*(1.0/(Beta*Beta) - 1.0);
    }
  }

  return Mass2;
}

float StVecMesonCut::getGlobalMass2(StPicoDst *picoDst, int i_track)
{
  float Mass2 = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    float Beta = tofTrack->btofBeta();
    float Momentum = picoTrack->gMom().Mag(); // global momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && Beta != 0)
    {
      Mass2 = Momentum*Momentum*(1.0/(Beta*Beta) - 1.0);
    }
  }

  return Mass2;
}

bool StVecMesonCut::passTrackBasic(StPicoTrack *picoTrack)
{
  // nHitsFit cut
  if(picoTrack->nHitsFit() < vmsa::mHitsFitTPCMin)
  {
    return kFALSE;
  }

  // nHitsRatio cut
  if(picoTrack->nHitsMax() <= vmsa::mHitsMaxTPCMin)
  {
    return kFALSE;
  }
  if((float)picoTrack->nHitsFit()/(float)picoTrack->nHitsMax() < vmsa::mHitsRatioTPCMin)
  {
    return kFALSE;
  }

  // eta cut
  // float eta = picoTrack->pMom().pseudoRapidity();
  float eta = picoTrack->pMom().Eta();
  if(fabs(eta) > vmsa::mEtaMax)
  {
    return kFALSE;
  }

  return kTRUE;
}


#if 0
bool StVecMesonCut::passSigPionCut(StPicoTrack* track, float scale_nSigma_factor)
{
  float nSigmaPion = track->nSigmaPion();
  if(fabs(nSigmaPion*scale_nSigma_factor) > vmsa::mNSigmaPionMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StVecMesonCut::passSigKaonCut(StPicoTrack* track, float scale_nSigma_factor)
{
  float nSigmaKaon = track->nSigmaKaon();
  if(fabs(nSigmaKaon*scale_nSigma_factor) > vmsa::mNSigmaKaonMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

bool StVecMesonCut::passSigProntonCut(StPicoTrack* track, float scale_nSigma_factor)
{
  float nSigmaProton = track->nSigmaProton();
  if(fabs(nSigmaProton*scale_nSigma_factor) > vmsa::mNSigmaProtonMax)
  {
    return kFALSE;
  }
  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StVecMesonCut::passTrackEP(StPicoTrack *track)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  // dca cut for event plane reconstruction: 200GeV = 3.0, BES = 1.0
  if(track->dca() > vmsa::mDcaEPMax[mEnergy])
  {
    return kFALSE;
  }

  // pt cut 0.2 - 2.0 GeV/c
  // float pt = track->pMom().perp();
  // float p  = track->pMom().mag();
  float pt = track->pMom().Perp();
  float p  = track->pMom().Mag();
  if(!(pt > vmsa::mPrimPtMin[mEnergy] && pt < vmsa::mPrimPtMax && p < vmsa::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------
bool StVecMesonCut::passTrackPhi(StPicoTrack *track)
{
  if(!track) return kFALSE;

  if(!passTrackBasic(track)) return kFALSE;

  // dca cut for phi-meson TTree production: 3.0
  if(track->dca() > vmsa::mDcaTrMax_phi)
  {
    return kFALSE;
  }

  // primary pt and momentum cut: PtMin = 0.1
  // if(!(track->pMom().perp() > vmsa::mGlobPtMin && track->pMom().mag() < vmsa::mPrimMomMax))
  if(!(track->pMom().Perp() > vmsa::mGlobPtMin && track->pMom().Mag() < vmsa::mPrimMomMax))
  {
    return kFALSE;
  }

  return kTRUE;
}
#endif
