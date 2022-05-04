#ifndef StSpinAlignmentCons_h
#define StSpinAlignmentCons_h

#include <string>

namespace globCons
{
  const int mNumBeamType = 2; // 0 for ZrZr200GeV_2018, 1 for RuRu200GeV_2018
  const std::string mBeamType[mNumBeamType] = {"ZrZr200GeV_2018", "RuRu200GeV_2018"};
  const float mBeamEnergy[mNumBeamType] = {200.0, 200.0};
  const int mBeamYear[mNumBeamType] = {2018, 2018};
  const int mMaxRunIndex = 4000; // maximum number of runIndex

  // event cut
  float const mVzMaxMap[mNumBeamType] = {40.0,70.0}; // 0: ZrZr200GeV_2018, 1: RuRu200GeV_2018
  float const mVrMax[mNumBeamType] = {2.0,2.0};
  float const mVzVpdDiffMax[mNumBeamType] = {3.0,3.0}; // 3.0
  int const mMatchedToFMin[mNumBeamType] = {2,2}; // 2

  // track cut
  float const mSigScaleMap[mNumBeamType] = {1.0,1.0};
  float const mDcaTrQAMax = 3.0; // use primary tracks run-by-run QA 
  int const mHitsDedxMin = 5;
  int const mHitsFitTPCMin = 15;
  int const mHitsMaxTPCMin = 0;
  float const mHitsRatioTPCMin = 0.51;
  float const mEtaMax = 1.0;
  float const mPrimPtMin[mNumBeamType] = {0.2,0.2}; // for event plane reconstruction and for pion, kaon, proton: 0.2 for BES
  float const mPrimPtMax = 2.0;
  float const mPrimPtWeight = 2.0;
  float const mPrimMomMax = 10.0; // also use for gMom
  float const mGlobPtMin = 0.1; // for phi, Lambda, K0s

  // Event Plane Maker
  int const mTrackMin = 2;
  int const mTrackMin_Full = 4;
  float const mEta_Gap = 0.05;
  float const mShiftOrder[5] = {2.0, 4.0, 6.0, 8.0, 10.0};
}

#endif
