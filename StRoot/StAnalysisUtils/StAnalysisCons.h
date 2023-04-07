#ifndef StAnalysisCons_h
#define StAnalysisCons_h

#include <string>
// #include "TString.h"

namespace anaUtils
{
  //--------------------------------------------------
  // shared in all analysis module
  const int mNumBeamUtils = 3; // 0 for ZrZr200GeV_2018, 1 for RuRu200GeV_2018, 2 for Fxt3p85GeV_2018

  // event cuts | copied from Isobar CME analysis & 3 GeV flow analysis
  const double mVzMin[mNumBeamUtils]                 = {-35.0, -35.0, 198.0};
  const double mVzCtr[mNumBeamUtils]                 = {0.0, 0.0, 200.0};
  const double mVzMax[mNumBeamUtils]                 = {25.0, 25.0, 202.0};
  const double mVzVpdDiffMax[mNumBeamUtils]          = {5.0, 5.0, 5.0};
  const double mVxCtr[mNumBeamUtils]                 = {0.0, 0.0, 0.0};
  const double mVyCtr[mNumBeamUtils]                 = {0.0, 0.0, -2.0};
  const double mVrMax[mNumBeamUtils]                 = {2.0, 2.0, 2.0};
  const unsigned short mMatchedToFMin[mNumBeamUtils] = {2, 2, 2};
  // const int mNumVzBin[mNumBeamUtils]                 = {2, 2, 2}; // 0 for vz < 0 & 1 for vz >= 0

  // track cuts: Basic 
  const int mHitsFitTpcMin[mNumBeamUtils]      = {15, 15, 15}; // Default >= 15 | SysError: >= 20
  const int mHitsMaxTpcMin[mNumBeamUtils]      = {0, 0, 0};
  const double mHitsRatioTpcMin[mNumBeamUtils] = {0.52, 0.52, 0.52};
  // const int   mHitsDedxTpcMin[mNumBeamUtils]  = {5, 5}; // not used

  // track cuts: RunQA
  const double mDcaQaMax[mNumBeamUtils]    = {3.0, 3.0, 3.0}; // use primary tracks run-by-run QA
  const double mEtaQaMin[mNumBeamUtils]    = {-1.0, -1.0, -2.0};
  // const double mEtaQaCtr[mNumBeamUtils]    = {0.0, 0.0, -1.05}; // Temp for Fxt3p85GeV_2018
  const double mEtaQaMax[mNumBeamUtils]    = {1.0, 1.0, 0.0};
  const double mPrimPtQaMin[mNumBeamUtils] = {0.2, 0.2, 0.2}; 

  // track cuts: TPC Event Plane Maker
  const double mDcaEpMax[mNumBeamUtils]       = {3.0, 3.0, 3.0};
  const double mEtaEpMin[mNumBeamUtils]       = {-1.0, -1.0, -2.0};
  const double mEtaEpCtr[mNumBeamUtils]       = {0.0, 0.0, -1.05};
  const double mEtaEpMax[mNumBeamUtils]       = {1.0, 1.0, 0.0};
  const double mEtaEpGap[mNumBeamUtils]       = {0.05, 0.05, 0.05}; // eta gap between Tracks and EP is 0.05 && eta gap between East and West TPC EP is 0.1
  const double mPrimPtEpMin[mNumBeamUtils]    = {0.2, 0.2, 0.2};
  const double mPrimPtEpMax[mNumBeamUtils]    = {2.0, 2.0, 2.0};
  const double mPrimPtEpWeight[mNumBeamUtils] = {2.0, 2.0, 2.0};
  const double mPrimMomEpMax[mNumBeamUtils]   = {10.0, 10.0, 10.0};
  const int mNumTrackEpMin[mNumBeamUtils]     = {2, 2, 2};

  // track cuts: Kaon
  const double mDcaKaonMax[mNumBeamUtils]     = {3.0, 3.0, 3.0}; // Default: <= 3.0 | SysError: <= 2.0 & <= 2.5
  const double mEtaKaonMin[mNumBeamUtils]     = {-1.0, -1.0, -2.0};
  const double mEtaKaonCtr[mNumBeamUtils]     = {0.0, 0.0, -1.05};
  const double mEtaKaonMax[mNumBeamUtils]     = {1.0, 1.0, 0.0};
  const double mPrimPtKaonMin[mNumBeamUtils]  = {0.2, 0.2, 0.2};
  const double mPrimMomKaonMax[mNumBeamUtils] = {10.0, 10.0, 10.0};
  const double mNSigKaonMin[mNumBeamUtils]    = {-3.0, -3.0, -3.0}; // Default: <= 2.5 | SysError: <= 3.0 & <= 2.0
  const double mNSigKaonMax[mNumBeamUtils]    = {3.0, 3.0, 3.0};
  const double mMass2KaonMin[mNumBeamUtils]   = {0.1, 0.1, 0.1};
  const double mMass2KaonMax[mNumBeamUtils]   = {0.4, 0.4, 0.4};

  // hit cuts: EPD Event Plane Maker
  const double mMipEpdEpMax[mNumBeamUtils] = {2.0, 2.0, 2.0};
  const double mMipEpdEpMin[mNumBeamUtils] = {0.3, 0.3, 0.3};

  // constants used in StPhiMesonMaker
  const double mMassPion     = 0.139570;
  const double mMassKaon     = 0.493677;
  const double mMassProton   = 0.938272;
  const double mMassDeuteron = 1.875613;
  const double mMassPhi      = 1.019455;
  const double mMassPhiMin   = 0.98;
  const double mMassPhiMax   = 1.08;

  // rapidity of center of mass system
  const double mRapCtrM[mNumBeamUtils] = {0.0, 0.0, -1.045};
}

#endif
