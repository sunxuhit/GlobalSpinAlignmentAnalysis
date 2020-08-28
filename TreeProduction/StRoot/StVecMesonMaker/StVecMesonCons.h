#ifndef StVecMesonCons_h
#define StVecMesonCons_h

#include <string>
#include "TString.h"
// #include "StarClassLibrary/SystemOfUnits.h"

namespace vmsa
{
  //--------------------------------------------------
  // used in TreeProduction
  int const NumBeamEnergy = 3;
  std::string const mBeamEnergy[NumBeamEnergy] = {"200GeV_2014","54GeV_2017","27GeV_2018"};
  float const mEnergyValue[NumBeamEnergy] = {200.0,54.0,27.0};
  int const mBeamYear[NumBeamEnergy] = {2014,2017,2018};

  // event cut
  float const mVzMaxMap[NumBeamEnergy] = {6.0,40.0,70.0}; // 0: 200GeV_2014 | 1: 54GeV_2017 | 2: 27GeV_2018 
  float const mVrMax[NumBeamEnergy] = {2.0,2.0,2.0};
  float const mVzVpdDiffMax[NumBeamEnergy] = {3.0,3.0,3.0}; // 3.0
  unsigned int const mMatchedToFMin[NumBeamEnergy] = {2,2,2}; // 2

  // track cut
  float const mSigScaleMap[NumBeamEnergy] = {1.0,1.0,1.0};
  float const mDcaTrQAMax = 3.0; // use primary tracks run-by-run QA 
  float const mDcaEPMax[NumBeamEnergy] = {3.0,1.0,1.0}; // for event plane reconstruction: 1.0 for BES
  float const mDcaTrMax = 1.0; // for pion, kaon, proton mDcaTrMax = 1.0 for flow
  float const mDcaTrMax_phi = 3.0; // for phi meson mDcaTrMax = 2.0 to fill a tree and apply an additional cut
  int const mHitsDedxMin = 5;
  int const mHitsFitTPCMin = 15;
  int const mHitsMaxTPCMin = 0;
  float const mHitsRatioTPCMin = 0.51;
  float const mEtaMax = 1.0;
  float const mPrimPtMin[NumBeamEnergy] = {0.2,0.2,0.2}; // for event plane reconstruction and for pion, kaon, proton: 0.2 for BES
  float const mGlobPtMin = 0.1; // for phi, Lambda, K0s
  float const mPrimPtMax = 2.0;
  float const mPrimPtWeight = 2.0;
  float const mPrimMomMax = 10.0; // also use for gMom
  float const mMass2Min = -10.0;
  // double const MAGFIELDFACTOR = kilogauss;
  int const mTrackMin = 2;
  int const mTrackMin_Full = 4;
  float const mToFYLocalMax = 1.8;
  float const mToFZLocalMax = 1.8;
  float const mNSigmaElectronMax = 2.5;
  float const mNSigmaPionMax = 2.5;
  float const mNSigmaKaonMax = 3.0;
  float const mNSigmaProtonMax = 2.5;
  float const mMassPion = 0.13957;
  float const mMassKaon = 0.49368;
  float const mMassProton = 0.93827;
  float const mSigKaon = 2.5;
  float const mNSigmaToF = 0.6;

  // used constant
  float const mEta_Gap = 0.05;
  float const mShiftOrder[5] = {2.0, 4.0, 6.0, 8.0, 10.0};

  int const pt_total = 25; // pt bin
  int const pt_start = 0;
  int const pt_stop  = 25;
  float const ptRawStart[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2};
  float const ptRawStop[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};
  double const pt_bin[pt_total+1] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};

  // mix event
  int const Bin_Centrality = 9;
  int const Bin_VertexZ = 10;
  int const Bin_Phi_Psi = 5;
  int const Buffer_depth = 5;
  TString const MixEvent[2] = {"SE","ME"};

  TString const vm_tree[2]  = {"PhiMesonEvent","KStarEvent"};
  TString const vm_branch[2] = {"phi_SpinAlignment_branch","KStar_SpinAlignment_branch"};
  int const mList_Delta = 20;
  //--------------------------------------------------

  std::string const mPID[3]   = {"Phi","KStar","K0S"};

  // ZDC-SMD constant
  std::string const mEastWest[2] = {"East","West"};
  std::string const mVertHori[2] = {"Vertical","Horizontal"};
}

#endif
