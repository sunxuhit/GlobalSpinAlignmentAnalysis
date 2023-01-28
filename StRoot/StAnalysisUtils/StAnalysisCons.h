#ifndef StAnalysisCons_h
#define StAnalysisCons_h

#include <string>
// #include "TString.h"

namespace anaUtils
{
  //--------------------------------------------------
  // shared in all analysis module
  const int mNumBeamUtils = 2; // 0 for ZrZr200GeV_2018, 1 for RuRu200GeV_2018, 2 for Fixed Target?

  // event cuts | copied from Isobar CME analysis
  const float mVzMin[mNumBeamUtils]                  = {-35.0, -35.0};
  const float mVzMax[mNumBeamUtils]                  = {25.0, 25.0};
  const float mVrMax[mNumBeamUtils]                  = {2.0, 2.0};
  const float mVzVpdDiffMax[mNumBeamUtils]           = {5.0, 5.0};
  const unsigned short mMatchedToFMin[mNumBeamUtils] = {2, 2};

  // track cuts: Basic 
  const int   mHitsFitTpcMin[mNumBeamUtils]   = {15, 15};
  const int   mHitsMaxTpcMin[mNumBeamUtils]   = {0, 0};
  const float mHitsRatioTpcMin[mNumBeamUtils] = {0.52, 0.52};
  // const int   mHitsDedxTpcMin[mNumBeamUtils]  = {5, 5}; // not used

  // track cuts: RunQA
  const float mDcaQaMax[mNumBeamUtils]    = {3.0,3.0}; // use primary tracks run-by-run QA
  const float mEtaQaMax[mNumBeamUtils]    = {1.0,1.0};
  const float mPrimPtQaMin[mNumBeamUtils] = {0.1,0.1}; 

  // track cuts: TPC Event Plane Maker
  const float mDcaEpMax[mNumBeamUtils]       = {3.0, 3.0};
  const float mEtaEpMax[mNumBeamUtils]       = {1.0, 1.0};
  const float mPrimPtEpMin[mNumBeamUtils]    = {0.2, 0.2};
  const float mPrimPtEpMax[mNumBeamUtils]    = {2.0, 2.0};
  const float mPrimPtEpWeight[mNumBeamUtils] = {2.0, 2.0};
  const float mPrimMomEpMax[mNumBeamUtils]   = {10.0, 10.0};
  const float mEtaEpGap[mNumBeamUtils]       = {0.05, 0.05}; // eta gap between Tracks and EP is 0.05 && eta gap between East and West TPC EP is 0.1

  // track cuts: phi Meson Maker
  // const float mDcaDauMax     = 3.0;
  // const float mEtaDauMax     = 1.0;
  // const float mPrimPtDauMin  = 0.2; // for phi decay daughter K+ and K-
  // const float mPrimMomDauMax = 10.0;
}

#endif
