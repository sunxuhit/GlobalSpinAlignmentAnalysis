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

  // event cut | copied from Isobar Blind Analysis
  const float mVzMax[mNumBeamType]        = {25.0,25.0};
  const float mVzMin[mNumBeamType]        = {-35.0,-35.0};
  const float mVrMax[mNumBeamType]        = {2.0,2.0};
  const float mVzVpdDiffMax[mNumBeamType] = {5.0,5.0};
  const unsigned short mMatchedToFMin     = 2;

  // track cut
  const int   mHitsFitTPCMin   = 15;
  const int   mHitsMaxTPCMin   = 0;
  const float mHitsRatioTPCMin = 0.51;
  // const int   mHitsDedxMin     = 5; // basic track cuts

  const float mPrimPtQaMin = 0.1; // cuts for track QA
  const float mDcaQaMax    = 3.0;
  const float mEtaQaMax    = 1.5;

  // const float mPrimPtEpMin    = 0.2; // for event plane reconstruction
  // const float mPrimPtEpMax    = 2.0;
  // const float mPrimPtEpWeight = 2.0;
  // const float mPrimMomEpMax   = 10.0;
  // const float mDcaEpMax       = 3.0;
  // const float mEtaEpMax       = 1.5;
  // const float mEtaEpGap       = 0.05;

  // const float mPrimPtDauMin  = 0.1; // for phi decay daughter K+ and K-
  // const float mPrimMomDauMax = 10.0;
  // const float mDcaDauMax     = 3.0;
  // const float mEtaDauMax     = 1.0;
}

#endif
