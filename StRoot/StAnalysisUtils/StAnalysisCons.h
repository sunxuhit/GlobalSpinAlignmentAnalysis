#ifndef StRunQACons_h
#define StRunQACons_h

#include <string>
#include "TString.h"
// #include "StarClassLibrary/SystemOfUnits.h"

namespace anaUtils
{
  //--------------------------------------------------
  // used in StAnalysisUtils
  const int mNumBeamUtils = 2; // 0 for ZrZr200GeV_2018, 1 for RuRu200GeV_2018, 2 for Fixed Target?

  // event cuts
  const float mVzMax[mNumBeamUtils]                  = {25.0, 25.0};
  const float mVzMin[mNumBeamUtils]                  = {-35.0, -35.0};
  const float mVrMax[mNumBeamUtils]                  = {2.0, 2.0};
  const float mVzVpdDiffMax[mNumBeamUtils]           = {5.0, 5.0};
  const unsigned short mMatchedToFMin[mNumBeamUtils] = {2, 2};

  // track cuts: Basic 
  const int   mHitsDedxMin     = 5;
  const int   mHitsFitTpcMin   = 15;
  const int   mHitsMaxTpcMin   = 0;
  const float mHitsRatioTpcMin = 0.51;

  // track cuts: RunQA
  const float mDcaQaMax    = 3.0; // use primary tracks run-by-run QA
  const float mEtaQaMax    = 1.0;
  const float mPrimPtQaMin = 0.1; 

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
