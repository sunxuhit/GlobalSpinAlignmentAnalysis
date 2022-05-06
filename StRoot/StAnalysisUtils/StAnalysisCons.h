#ifndef StRunQACons_h
#define StRunQACons_h

#include <string>
#include "TString.h"
// #include "StarClassLibrary/SystemOfUnits.h"

namespace anaUtils
{
  //--------------------------------------------------
  // used in StAnalysisUtils
  const int mNumBeamType = 2; // 0 for ZrZr200GeV_2018, 1 for RuRu200GeV_2018

  // event cut
  const float mVzMax[mNumBeamType]        = {70.0,70.0};
  const float mVzMin[mNumBeamType]        = {-70.0,-70.0};
  const float mVrMax[mNumBeamType]        = {2.0,2.0};
  const float mVzVpdDiffMax[mNumBeamType] = {3.0,3.0};
  const unsigned short mMatchedToFMin     = 2;

  // track cut
  // float const mSigScaleMap[mNumBeamType] = {1.0,1.0};
  const float mDcaTrQAMax      = 3.0; // use primary tracks run-by-run QA
  const int   mHitsDedxMin     = 5;
  const int   mHitsFitTPCMin   = 15;
  const int   mHitsMaxTPCMin   = 0;
  const float mHitsRatioTPCMin = 0.51;
  const float mEtaMax          = 1.0;

  const float mPrimPtMin[mNumBeamType] = {0.2,0.2}; // for event plane reconstruction and for pion, kaon, proton: 0.2 for BES
  const float mPrimPtMax = 2.0;
  const float mPrimPtWeight = 2.0;
  const float mPrimMomMax = 10.0; // also use for gMom
  const float mGlobPtMin = 0.1; // for phi, Lambda, K0s
}

#endif
