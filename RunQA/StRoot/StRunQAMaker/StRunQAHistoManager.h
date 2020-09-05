#ifndef StRunQAHistoManager_h
#define StRunQAHistoManager_h

#include "StMessMgr.h"
#include "TVector3.h"

class TH1F;
class TH2F;

class StRunQAHistoManager
{
  public:
    StRunQAHistoManager();
    virtual ~StRunQAHistoManager();

    //--------------QA---------------
    void initEventQA();
    void fillEventQA_RefMult(int refMult, int grefMult, int cent9, double reweight, int tofHits, int tofMatch, int cutSelection);
    void fillEventQA_Vertex(float vx, float vy, float vz, float vzVpd, int cutSelection);
    void fillEventQA_Trigger(int triggerBin, int cutSelection);
    void writeEventQA();

    void initTrackQA();
    void fillTrackQA_Kinematics(TVector3 pMom, TVector3 gMom, int cutSelection);
    void fillTrackQA_Quliaty(float gDca, int nHitsFit, int nHitsMax, int nHitsDEdx, int cutSelection);
    void fillTrackQA_PID(float mom, short charge, float dEdx, float beta, float mass2, int cutSelection);
    void writeTrackQA();
    //--------------QA---------------

  private:
    // QA Histograms
    // Event Level:
    TH1F *h_mRefMult[2]; // 0: before cuts | 1: after cuts
    TH1F *h_mGRefMult[2];
    TH2F *h_mRefMultGRefMult[2];
    TH1F *h_mCentrality9[2];
    TH2F *h_mRefMultTofMatch[2];
    TH2F *h_mRefMultTofHits[2];
    TH2F *h_mGRefMultTofMatch[2];
    TH2F *h_mGRefMultTofHits[2];
    TH2F *h_mVzVzVpd[2];
    TH1F *h_mDiffVzVzVpd[2];
    TH1F *h_mVertexZ[2];
    TH2F *h_mVertexXY[2];
    TH1F *h_mTriggerID[2];
    // Track Level:
    TH1F *h_mPrimPt[2];
    TH1F *h_mPrimEta[2];
    TH1F *h_mPrimPhi[2];
    TH1F *h_mGlobPt[2];
    TH1F *h_mGlobEta[2];
    TH1F *h_mGlobPhi[2];
    TH1F *h_mDca[2];
    TH1F *h_mNHitsFit[2];
    TH1F *h_mNHitsRatio[2];
    TH1F *h_mNHitsDEdx[2];
    TH2F *h_mDEdxMom[2];
    TH2F *h_mBetaMom[2];
    TH2F *h_mMass2Mom[2];

    std::string mCutsQA[2] = {"Before","After"};

  ClassDef(StRunQAHistoManager,1)
};
#endif
