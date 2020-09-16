#ifndef StRunQAProManager_h
#define StRunQAProManager_h

#include <TVector3.h>
#include <TString.h>
#include "StMessMgr.h"

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
    // Run-by-Run QA | x axis is RunIndex
    TProfile *p_mRefMult[2][10]; // 0: before cuts | 1: after cuts
    TProfile *p_mGRefMult[2][10]; // 0-8 for different triggerID | 9 for all triggers
    TProfile *p_mZdcX[2][10];
    TProfile *p_mVz[2][10];
    TProfile *p_mVr[2][10];
    TProfile *p_mGDca[2][10];
    TProfile *p_mNHitsFit[2][10];
    TProfile *p_mPrimPt[2][10];
    TProfile *p_mPrimEta[2][10];
    TProfile *p_mPrimPhi[2][10];
    TProfile *p_mGlobPt[2][10];
    TProfile *p_mGlobEta[2][10];
    TProfile *p_mGlobPhi[2][10];

    std::string mCutsQA[2] = {"Before","After"};

    ClassDef(StRunQAProManager,1)
};

#endif
