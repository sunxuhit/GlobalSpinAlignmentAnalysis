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
    void initEvtQA();
    void fillEvtQaRefMult(int triggerBin, int refMult, int grefMult, int cent9, double refWgt, int tofHits, int tofMatch, int cutSelection);
    void fillEvtQaVertex(int triggerBin, double vx, double vy, double vz, double vzVpd, int vzBin, int cutSelection);
    void fillEvtQaTrigger(int triggerBin, int cutSelection);
    void writeEvtQA();

    void initTrkQA();
    void fillTrkQaKinematics(int triggerBin, TVector3 pMom, TVector3 gMom, int cutSelection);
    void fillTrkQaQuliaty(int triggerBin, double gDca, int nHitsFit, int nHitsMax, int nHitsDEdx, int cutSelection);
    void fillTrkQaPID(int triggerBin, double mom, short charge, double dEdx, double beta, double mass2, int cutSelection);
    void fillTrkQaEpCut(int triggerBin, TVector3 pMom, bool isFull, bool isEast, bool isWest, int cutSelection);
    void fillTrkQaFlowCut(int triggerBin, TVector3 pMom, bool isFull, bool isEast, bool isWest, int cutSelection);
    void fillTrkQaKaonCut(int triggerBin, TVector3 pMom, double nSigKaon, bool isFull, bool isEast, bool isWest, int cutSelection);
    void fillTrkQaKaonAcptTree(int triggerBin, int cent9, int charge, double yLab, double yCms, double pt, double refWgt);
    void fillTrkQaKaonAcptSpin(int triggerBin, int cent9, int charge, double yLab, double yCms, double pt, double refWgt);
    void writeTrkQA();
    //--------------QA---------------

  private:
    static const int mNumCuts = 2; // 0: before cuts | 1: after cuts
    static const int mNumTriggerBins = 10; // 0-8 for different triggerID | 9 for all triggers
    static const int mNumCentrality = 9; // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    static const int mNumPtBinQA    = 500; // 250 bins from -0.05 to 9.95 GeV/c
    static const int mNumRapBinQA   = 250; // 250 bins from -2.5 to 2.5

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
    TH2F *h_mVzVzBin[mNumCuts][mNumTriggerBins];
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
    // test StAnalysisCuts
    TH1F *h_mPrimEtaEpFull[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEtaEpEast[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEtaEpWest[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEtaFlowFull[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEtaFlowEast[mNumCuts][mNumTriggerBins];
    TH1F *h_mPrimEtaFlowWest[mNumCuts][mNumTriggerBins];
    TH2F *h_mPrimEtaNSigKaonFull[mNumCuts][mNumTriggerBins];
    TH2F *h_mPrimEtaNSigKaonEast[mNumCuts][mNumTriggerBins];
    TH2F *h_mPrimEtaNSigKaonWest[mNumCuts][mNumTriggerBins];

    // Kaon Acceptance
    TH2F *h_mAcptTreeLabKp[mNumCentrality][mNumTriggerBins]; // y vs. pT
    TH2F *h_mAcptTreeCmsKp[mNumCentrality][mNumTriggerBins];
    TH2F *h_mAcptTreeLabKm[mNumCentrality][mNumTriggerBins];
    TH2F *h_mAcptTreeCmsKm[mNumCentrality][mNumTriggerBins];
    TH2F *h_mAcptSpinLabKp[mNumCentrality][mNumTriggerBins]; // y vs. pT
    TH2F *h_mAcptSpinCmsKp[mNumCentrality][mNumTriggerBins];
    TH2F *h_mAcptSpinLabKm[mNumCentrality][mNumTriggerBins];
    TH2F *h_mAcptSpinCmsKm[mNumCentrality][mNumTriggerBins];

    std::string str_mCutStatus[mNumCuts] = {"Bf","Af"};

    const int mType;


  ClassDef(StRunQAHistoManager,1)
};
#endif
