#ifndef StRunQAProManager_h
#define StRunQAProManager_h

#include <TVector3.h>
// #include "StMessMgr.h"

class TProfile;

class StRunQAProManager
{
  public:
    StRunQAProManager();
    virtual ~StRunQAProManager();

    // Run-by-Run QA
    void initRunQA();
    void fillRunQA_Event(int triggerBin, int runIdenx, float refMult, float grefMult, float zdcX, float vx, float vy, float vz, int cutSelection);
    void fillRunQA_Track(int triggerBin, int runIdenx, float gDca, int nHitsFit, TVector3 pMom, TVector3 gMom, int cutSelection);
    void writeRunQA();

  private:
    const int mNumCuts = 2;
    const int mNumTriggerBins = 10;
    const std::string mCutStatus[mNumCuts] = {"Bf","Af"};

    // Run-by-Run QA | x axis is RunIndex
    TProfile *p_mRefMult[mNumCuts][mNumTriggerBins]; // 0: before cuts | 1: after cuts
    TProfile *p_mGRefMult[mNumCuts][mNumTriggerBins]; // 0-8 for different triggerID | 9 for all triggers
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

    ClassDef(StRunQAProManager,1)
};

#endif
