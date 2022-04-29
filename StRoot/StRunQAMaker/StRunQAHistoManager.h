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
    void fillEventQA_RefMult(int triggerBin, int refMult, int grefMult, int cent9, double reweight, int tofHits, int tofMatch, int cutSelection);
    void fillEventQA_Vertex(int triggerBin, float vx, float vy, float vz, float vzVpd, int cutSelection);
    void fillEventQA_Trigger(int triggerBin, int cutSelection);
    void writeEventQA();

    void initTrackQA();
    void fillTrackQA_Kinematics(int triggerBin, TVector3 pMom, TVector3 gMom, int cutSelection);
    void fillTrackQA_Quliaty(int triggerBin, float gDca, int nHitsFit, int nHitsMax, int nHitsDEdx, int cutSelection);
    void fillTrackQA_PID(int triggerBin, float mom, short charge, float dEdx, float beta, float mass2, int cutSelection);
    void writeTrackQA();
    //--------------QA---------------

  private:
    // QA Histograms
    // Event Level:
    TH1F *h_mTriggerID[2];
    TH1F *h_mRefMult[2][10]; // 0: before cuts | 1: after cuts
    TH1F *h_mGRefMult[2][10]; // 0-8 for different triggerID | 9 for all triggers
    TH2F *h_mRefMultGRefMult[2][10];
    TH1F *h_mCentrality9[2][10];
    TH2F *h_mTofMatchRefMult[2][10];
    TH2F *h_mTofHitsRefMult[2][10];
    TH2F *h_mTofMatchGRefMult[2][10];
    TH2F *h_mTofHitsGRefMult[2][10];
    TH2F *h_mVzVzVpd[2][10];
    TH1F *h_mDiffVzVzVpd[2][10];
    TH1F *h_mVertexZ[2][10];
    TH2F *h_mVertexXY[2][10];
    // Track Level:
    TH1F *h_mPrimPt[2][10];
    TH1F *h_mPrimEta[2][10];
    TH1F *h_mPrimPhi[2][10];
    TH1F *h_mGlobPt[2][10];
    TH1F *h_mGlobEta[2][10];
    TH1F *h_mGlobPhi[2][10];
    TH1F *h_mDca[2][10];
    TH1F *h_mNHitsFit[2][10];
    TH1F *h_mNHitsRatio[2][10];
    TH1F *h_mNHitsDEdx[2][10];
    TH2F *h_mDEdxMom[2][10];
    TH2F *h_mBetaMom[2][10];
    TH2F *h_mMass2Mom[2][10];

    std::string mCutsQA[2] = {"Before","After"};

  ClassDef(StRunQAHistoManager,1)
};
#endif
