#ifndef StVecMesonHistoManager_h
#define StVecMesonHistoManager_h

#include "StMessMgr.h"
#include "TVector3.h"

class TH1F;
class TH2F;

class StVecMesonHistoManager
{
  public:
    StVecMesonHistoManager();
    virtual ~StVecMesonHistoManager();

    //--------------QA---------------
    void initEventQA();
    void fillEventQA_RefMult(float refMult, float cent9, float tofHits, float tofMatch, int cutSelection);
    void fillEventQA_Vertex(float vx, float vy, float vz, float vzVpd, int cutSelection);
    void writeEventQA();

    void initTrackQA();
    void fillTrackQA_Kinematics(TVector3 pMom, TVector3 gMom, int cutSelection);
    void fillTrackQA_Quliaty(float gDca, int nHitsFit, int nHitsMax, int nHitsDEdx, int cutSelection);
    void fillTrackQA_PID(float mom, short charge, float dEdx, float beta, float mass2, int cutSelection);
    void writeTrackQA();
    //--------------QA---------------

    //--------------ZDC EP---------------
    void initZdcGainCorr();
    void fillZdcGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, float zdcsmd);
    void writeZdcGainCorr();
    //--------------ZDC EP---------------

    /*
    void InitEP();
    void FillEP_Sub(Float_t Psi2East_ReCenter, Float_t Psi2East_Shift, Float_t Psi2West_ReCenter, Float_t Psi2West_Shift);
    void FillEP_Ran(Float_t Psi2RanA_ReCenter, Float_t Psi2RanA_Shift, Float_t Psi2RanB_ReCenter, Float_t Psi2RanB_Shift, Float_t Psi2Full_ReCenter, Float_t Psi2Full_Shift);
    void WriteEP();
    */
    
  private:
    // QA Histograms
    // Event Level:
    TH1F *h_mRefMult[2]; // 0: before cuts | 1: after cuts
    TH1F *h_mCentrality9[2];
    TH2F *h_mRefMultTofMatch[2];
    TH2F *h_mRefMultTofHits[2];
    TH2F *h_mVzVzVpd[2];
    TH1F *h_mDiffVzVzVpd[2];
    TH1F *h_mVertexZ[2];
    TH2F *h_mVertexXY[2];
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

    //--------------ZDC EP---------------
    TH2F *h_mZdcGainCorr[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y); | x-axis: runIndex | y-axis: ADC
    //--------------ZDC EP---------------

    std::string mCutsQA[2] = {"Before","After"};

  ClassDef(StVecMesonHistoManager,1)
};
#endif
