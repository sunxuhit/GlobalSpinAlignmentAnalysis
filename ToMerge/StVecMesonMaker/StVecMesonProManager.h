#ifndef StVecMesonProManager_h
#define StVecMesonProManager_h

#include <TVector3.h>
#include <TString.h>
#include "StMessMgr.h"

class TProfile;

class StVecMesonProManager
{
  public:
    StVecMesonProManager();
    virtual ~StVecMesonProManager();

    // Run-by-Run QA
    void initRunQA();
    void fillRunQA_Event(int runIdenx, float refMult, float zdcX, float vx, float vy, float vz, int cutSelection);
    void fillRunQA_Track(int runIdenx, float gDca, int nHitsFit, TVector3 pMom, TVector3 gMom, int cutSelection);
    void writeRunQA();

  private:
    // Run-by-Run QA | x axis is RunIndex
    TProfile *p_mQA_RefMult[2]; // 0: before cuts | 1: after cuts
    TProfile *p_mQA_ZdcX[2];
    TProfile *p_mQA_Vz[2];
    TProfile *p_mQA_Vr[2];
    TProfile *p_mQA_gDca[2];
    TProfile *p_mQA_nHitsFit[2];
    TProfile *p_mQA_PrimPt[2];
    TProfile *p_mQA_PrimEta[2];
    TProfile *p_mQA_PrimPhi[2];
    TProfile *p_mQA_GlobPt[2];
    TProfile *p_mQA_GlobEta[2];
    TProfile *p_mQA_GlobPhi[2];

    std::string mCutsQA[2] = {"Before","After"};

    ClassDef(StVecMesonProManager,1)
};

#endif
