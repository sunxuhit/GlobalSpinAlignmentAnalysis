#ifndef StRunQAProManager_h
#define StRunQAProManager_h

#include <TVector3.h>
// #include "StMessMgr.h"

class TProfile;

class StRunQAProManager
{
  public:
    StRunQAProManager(int beamType);
    virtual ~StRunQAProManager();

    // Run-by-Run QA
    void initRunQA();
    void fillRunQA_Event(int triggerBin, int runIndex, double refMult, double grefMult, double zdcX, double vx, double vy, double vz, int cutSelection);
    void fillRunQA_Track(int triggerBin, int runIndex, double gDca, int nHitsFit, TVector3 pMom, TVector3 gMom, int cutSelection);
    void writeRunQA();

  private:
    static const int mNumCuts = 2; // 0: before cuts | 1: after cuts
    static const int mNumTriggerBins = 10; // 0-8 for different triggerID | 9 for all triggers

    // Run-by-Run QA | x axis is RunIndex
    TProfile *p_mRefMult[mNumCuts][mNumTriggerBins];
    TProfile *p_mGRefMult[mNumCuts][mNumTriggerBins];
    TProfile *p_mZdcX[mNumCuts][mNumTriggerBins];
    TProfile *p_mVz[mNumCuts][mNumTriggerBins];
    TProfile *p_mVr[mNumCuts][mNumTriggerBins];
    TProfile *p_mGDca[mNumCuts][mNumTriggerBins];
    TProfile *p_mNHitsFit[mNumCuts][mNumTriggerBins];
    TProfile *p_mPrimPt[mNumCuts][mNumTriggerBins];
    TProfile *p_mPrimEta[mNumCuts][mNumTriggerBins];
    TProfile *p_mPrimPhi[mNumCuts][mNumTriggerBins];
    TProfile *p_mGlobPt[mNumCuts][mNumTriggerBins];
    TProfile *p_mGlobEta[mNumCuts][mNumTriggerBins];
    TProfile *p_mGlobPhi[mNumCuts][mNumTriggerBins];

    std::string str_mCutStatus[mNumCuts] = {"Bf","Af"};

    const int mType;

    ClassDef(StRunQAProManager,1)
};

#endif
