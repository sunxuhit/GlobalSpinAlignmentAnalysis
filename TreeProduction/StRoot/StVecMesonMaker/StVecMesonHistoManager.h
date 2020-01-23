#ifndef StVecMesonHistoManager_h
#define StVecMesonHistoManager_h

#include "StMessMgr.h"

class TH1F;
class TH2F;

class StVecMesonHistoManager
{
  public:
    StVecMesonHistoManager();
    virtual ~StVecMesonHistoManager();

    void Init_EventQA();
    void Fill_EventQA_RefMult(float refMult, float cent9, float tofHits, float tofMatch, int cutSelection);
    void Fill_EventQA_Vertex(float vx, float vy, float vz, float vzVpd, int cutSelection);
    void Write_EventQA();

    void Init_TrackQA();
    void Write_TrackQA();

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
    TH1F *h_mPt[2];
    TH1F *h_mEta[2];
    TH1F *h_mPhi[2];
    TH1F *h_mDca[2];
    TH1F *h_mNHitsFit[2];
    TH1F *h_mNHitsRatio[2];
    TH1F *h_mNHitsDEdx[2];
    TH2F *h_mDEdxMom[2];
    TH2F *h_mNSigmaPionMom[2];
    TH2F *h_mBetaMom[2];
    TH2F *h_mMass2Mom[2];

    /*
    TH1F *h_mEastRaw;
    TH1F *h_mWestRaw;
    TH1F *h_mFullRaw;

    // event plane distribution
    TH1F *h_mEastReCenter;
    TH1F *h_mWestReCenter;
    TH1F *h_mRanAReCenter;
    TH1F *h_mRanBReCenter;
    TH1F *h_mFullReCenter;

    TH1F *h_mEastShift;
    TH1F *h_mWestShift;
    TH1F *h_mRanAShift;
    TH1F *h_mRanBShift;
    TH1F *h_mFullShift;
    */

    std::string mCutsQA[2] = {"Before","After"};

  ClassDef(StVecMesonHistoManager,1)
};
#endif
