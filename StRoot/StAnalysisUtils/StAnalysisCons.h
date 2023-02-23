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
  const double mVzMin[mNumBeamUtils]                 = {-35.0, -35.0};
  const double mVzMax[mNumBeamUtils]                 = {25.0, 25.0};
  const double mVrMax[mNumBeamUtils]                 = {2.0, 2.0};
  const double mVzVpdDiffMax[mNumBeamUtils]          = {5.0, 5.0};
  const unsigned short mMatchedToFMin[mNumBeamUtils] = {2, 2};
  // const int mNumVzBin[mNumBeamUtils]                 = {2, 2}; // 0 for vz < 0 & 1 for vz >= 0

  // track cuts: Basic 
  const int mHitsFitTpcMin[mNumBeamUtils]      = {15, 15};
  const int mHitsMaxTpcMin[mNumBeamUtils]      = {0, 0};
  const double mHitsRatioTpcMin[mNumBeamUtils] = {0.52, 0.52};
  // const int   mHitsDedxTpcMin[mNumBeamUtils]  = {5, 5}; // not used

  // track cuts: RunQA
  const double mDcaQaMax[mNumBeamUtils]    = {3.0,3.0}; // use primary tracks run-by-run QA
  const double mEtaQaMax[mNumBeamUtils]    = {1.0,1.0};
  const double mPrimPtQaMin[mNumBeamUtils] = {0.1,0.1}; 

  // track cuts: TPC Event Plane Maker
  const double mDcaEpMax[mNumBeamUtils]       = {3.0, 3.0};
  const double mEtaEpMax[mNumBeamUtils]       = {1.0, 1.0};
  const double mPrimPtEpMin[mNumBeamUtils]    = {0.2, 0.2};
  const double mPrimPtEpMax[mNumBeamUtils]    = {2.0, 2.0};
  const double mPrimPtEpWeight[mNumBeamUtils] = {2.0, 2.0};
  const double mPrimMomEpMax[mNumBeamUtils]   = {10.0, 10.0};
  const double mEtaEpGap[mNumBeamUtils]       = {0.05, 0.05}; // eta gap between Tracks and EP is 0.05 && eta gap between East and West TPC EP is 0.1
  const int mNumTrackEpMin[mNumBeamUtils]     = {2, 2};

  // track cuts: Kaon
  const double mDcaKaonMax[mNumBeamUtils]     = {3.0, 3.0};
  const double mEtaKaonMax[mNumBeamUtils]     = {1.0, 1.0};
  const double mPrimPtKaonMin[mNumBeamUtils]  = {0.2, 0.2};
  const double mPrimMomKaonMax[mNumBeamUtils] = {10.0, 10.0};
  const double mNSigKaonMin[mNumBeamUtils]   = {-2.5, -2.5};
  const double mNSigKaonMax[mNumBeamUtils]   = {2.5, 2.5};

  // hit cuts: EPD Event Plane Maker
  const double mMipEpdEpMax[mNumBeamUtils] = {2.0, 2.0};
  const double mMipEpdEpMin[mNumBeamUtils] = {0.3, 0.3};
}

#endif
