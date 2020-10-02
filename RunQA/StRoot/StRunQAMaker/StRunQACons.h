#ifndef StRunQACons_h
#define StRunQACons_h

#include <string>
#include "TString.h"
// #include "StarClassLibrary/SystemOfUnits.h"

namespace runQA
{
  //--------------------------------------------------
  // used in RunQA
  int const NumBeamEnergy = 3;
  std::string const mBeamEnergy[NumBeamEnergy] = {"200GeV_2014","54GeV_2017","27GeV_2018"};
  float const mEnergyValue[NumBeamEnergy] = {200.0,54.0,27.0};
  int const mBeamYear[NumBeamEnergy] = {2014,2017,2018};

  // event cut
  float const mVzMaxMap[NumBeamEnergy] = {6.0,40.0,70.0}; // 0: 200GeV_2014 | 1: 54GeV_2017 | 2: 27GeV_2018 
  float const mVrMax[NumBeamEnergy] = {2.0,2.0,2.0};
  float const mVzVpdDiffMax[NumBeamEnergy] = {3.0,3.0,3.0}; // 3.0
  int const mMatchedToFMin[NumBeamEnergy] = {2,2,2}; // 2

  // track cut
  float const mSigScaleMap[NumBeamEnergy] = {1.0,1.0,1.0};
  float const mDcaTrQAMax = 3.0; // use primary tracks run-by-run QA 
  int const mHitsDedxMin = 5;
  int const mHitsFitTPCMin = 15;
  int const mHitsMaxTPCMin = 0;
  float const mHitsRatioTPCMin = 0.51;
  float const mEtaMax = 1.0;
  float const mPrimPtMin[NumBeamEnergy] = {0.2,0.2,0.2}; // for event plane reconstruction and for pion, kaon, proton: 0.2 for BES
  float const mPrimPtMax = 2.0;
  float const mPrimPtWeight = 2.0;
  float const mPrimMomMax = 10.0; // also use for gMom
  float const mGlobPtMin = 0.1; // for phi, Lambda, K0s

  int const mNumOfRunIndex = 4000;
}

#endif
