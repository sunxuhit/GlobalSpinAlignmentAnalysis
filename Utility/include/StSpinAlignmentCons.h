#ifndef StSpinAlignmentCons_h
#define StSpinAlignmentCons_h

#include <string>

namespace globCons
{
  const int mNumBeamType                        = 3; // 0 for ZrZr200GeV_2018, 1 for RuRu200GeV_2018, 2 for Fxt3p85GeV_2018
  const std::string str_mBeamType[mNumBeamType] = {"ZrZr200GeV_2018", "RuRu200GeV_2018", "Fxt3p85GeV_2018"};
  const double mBeamEnergy[mNumBeamType]        = {200.0, 200.0, 3.85};
  const int mBeamYear[mNumBeamType]             = {2018, 2018, 2018};
  const int mMaxRunIndex[mNumBeamType]          = {1500, 1500, 300}; // maximum number of runIndex
  const int mRunIndexLo[mNumBeamType]           = {0,    2500, 0}; 
  const int mRunIndexHi[mNumBeamType]           = {1500, 4000, 300};
}

#endif
