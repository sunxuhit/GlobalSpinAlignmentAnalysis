#ifndef StRunQAHistoManager_h
#define StRunQAHistoManager_h

#include "TVector3.h"
// #include "StMessMgr.h"

class TH1F;
class TH2F;

class StRunQAHistoManager
{
  public:
    StRunQAHistoManager(int beamType);
    virtual ~StRunQAHistoManager();

    //--------------QA---------------
    void initEventQA();
    void fillEventQA_RefMult(int triggerBin, int refMult, int grefMult, int cent9, double reweight, int tofHits, int tofMatch, int cutSelection);
    void fillEventQA_Vertex(int triggerBin, double vx, double vy, double vz, double vzVpd, int cutSelection);
    void fillEventQA_Trigger(int triggerBin, int cutSelection);
    void writeEventQA();

    void initTrackQA();
    void fillTrackQA_Kinematics(int triggerBin, TVector3 pMom, TVector3 gMom, int cutSelection);
    void fillTrackQA_Quliaty(int triggerBin, double gDca, int nHitsFit, int nHitsMax, int nHitsDEdx, int cutSelection);
    void fillTrackQA_PID(int triggerBin, double mom, short charge, double dEdx, double beta, double mass2, int cutSelection);
    void fillTrackQA_EpCut(int triggerBin, TVector3 pMom, bool isFull, bool isEast, bool isWest, int cutSelection);
    void fillTrackQA_FlowCut(int triggerBin, TVector3 pMom, bool isFull, bool isEast, bool isWest, int cutSelection);
    void fillTrackQA_KaonCut(int triggerBin, TVector3 pMom, double nSigKaon, bool isFull, bool isEast, bool isWest, int cutSelection);
    void writeTrackQA();
    //--------------QA---------------

  private:
    static const int mNumCuts = 2; // 0: before cuts | 1: after cuts
    static const int mNumTriggerBins = 10; // 0-8 for different triggerID | 9 for all triggers
    static const double mVzQaMin[3] = {-100.0, -100.0, 150.0};
    static const double mVzQaMax[3] = {100.0, 100.0, 250.0}

    // QA Histograms
    // Event Level:
    TH1F *h_mTriggerId[mNumCuts];
    TH1F *h_mRefMult[mNumCuts][mNumTriggerBins];
    TH1F *h_mGRefMult[mNumCuts][mNumTriggerBins];
    TH2F *h_mRefMultGRefMult[mNumCuts][mNumTriggerBins];
    TH1F *h_mCentrality9[mNumCuts][mNumTriggerBins];
    TH2F *h_mTofMatchRefMult[mNumCuts][mNumTriggerBins];
    TH2F *h_mTofHitsRefMult[mNumCuts][mNumTriggerBins];
    TH2F *h_mTofMatchGRefMult[mNumCuts][mNumTriggerBins];
    TH2F *h_mTofHitsGRefMult[mNumCuts][mNumTriggerBins];
    TH2F *h_mVzVzVpd[mNumCuts][mNumTriggerBins];
    TH1F *h_mDiffVzVzVpd[mNumCuts][mNumTriggerBins];
    TH1F *h_mVertexZ[mNumCuts][mNumTriggerBins];
    TH2F *h_mVertexXY[mNumCuts][mNumTriggerBins];
    // Track Level:
    TH1F *h_mPrimPt[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimPhi[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEta[mNumCuts][mNumTriggerBins];
    TH1F *h_mGlobPt[mNumCuts][mNumTriggerBins];
    TH1F *h_mGlobEta[mNumCuts][mNumTriggerBins];
    TH1F *h_mGlobPhi[mNumCuts][mNumTriggerBins];
    TH1F *h_mDca[mNumCuts][mNumTriggerBins];
    TH1F *h_mNHitsFit[mNumCuts][mNumTriggerBins];
    TH1F *h_mNHitsRatio[mNumCuts][mNumTriggerBins];
    TH1F *h_mNHitsDEdx[mNumCuts][mNumTriggerBins];
    TH2F *h_mMomDEdx[mNumCuts][mNumTriggerBins];
    TH2F *h_mMomMass2[mNumCuts][mNumTriggerBins];
    TH2F *h_mMomBeta[mNumCuts][mNumTriggerBins];

    TH1F *h_mPrimEtaEpFull[mNumCuts][mNumTriggerBins]; // test StAnalysisCuts
    TH1F *h_mPrimEtaEpEast[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEtaEpWest[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEtaFlowFull[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEtaFlowEast[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEtaFlowWest[mNumCuts][mNumTriggerBins];
    TH2F *h_mPrimEtaNSigKaonFull[mNumCuts][mNumTriggerBins];
    TH2F *h_mPrimEtaNSigKaonEast[mNumCuts][mNumTriggerBins];
    TH2F *h_mPrimEtaNSigKaonWest[mNumCuts][mNumTriggerBins];

    std::string str_mCutStatus[mNumCuts] = {"Bf","Af"};

    const int mType;


  ClassDef(StRunQAHistoManager,1)
};
#endif
