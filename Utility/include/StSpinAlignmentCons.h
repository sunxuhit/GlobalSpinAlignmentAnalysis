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
}

#endif
