#ifndef StPhiMesonTree_h
#define StPhiMesonTree_h

#include <map>
#include <vector>
#include <string>
#include "TObject.h"

class TTree;
class TH2F;
class TVector2;
class TVector3;

class StPicoDst;

class StAnalysisUtils;
class StAnalysisCut;
class StPhiMesonEvent;
class StPhiMesonTrack;

class StPhiMesonTree : public TObject
{
  public:
    StPhiMesonTree(int beamType);
    virtual ~StPhiMesonTree();

    void initPhiTree();
    void clearPhiMixBuffer(int cent9, int vzBin, int PsiBin);

    void fillPhiTree(StPicoDst* picoDst, int flagME);
    void recoPhi(int cent9, int vzBin, int PsiBin); // reconstruct phi meson in the same event
    void mixPhi(int cent9, int vzBin, int PsiBin); // reconstruct phi meson in the mixed event

    void writePhiTree();

    int getPhiMixKey(int cent9, int vzBin, int PsiBin, int evtBin); // return 1000*cent9 + 100*vzBin + 10*PsiBin + evtBin
    int getVzMixBin(double vz);
    int getPsiMixBin(double Psi, int epOrder);
    void getPhiEvtSize(int cent9, int vzBin, int PsiBin);

    // set event info
    void clearEvtInfo();
    void setEvtInfo(int runIdx, int cent9, int cent16, double refwgt, double vz, double PsiShiftFull);
    void setZdcQ1Flag(int flagEp); // ZDC Flag
    void setZdcQ1Vec(TVector2 Q1VecZdcShiftEast,TVector2 Q1VecZdcShiftWest, TVector2 Q1VecZdcShiftFull); // Shift Corrected ZDC Q1Vector: East & West & Full
    void setEpdQ1SideFlag(int flagEp); // EPD Side FLag
    void setEpdQ1SideVec(TVector2 Q1VecEpdSideShiftEast,TVector2 Q1VecEpdSideShiftWest, TVector2 Q1VecEpdSideShiftFull); // Shift Corrected EPD Side Q1Vector: East & West & Full
    void setEpdQ1Grp0Flag(int flagEp); // EPD Grp0 FLag
    void setEpdQ1Grp0Vec(TVector2 Q1VecEpdGrp0ShiftEast,TVector2 Q1VecEpdGrp0ShiftWest, TVector2 Q1VecEpdGrp0ShiftFull); // Shift Corrected EPD Grp0 Q1Vector: East & West & Full
    void setEpdQ1Grp1Flag(int flagEp); // EPD Grp1 FLag
    void setEpdQ1Grp1Vec(TVector2 Q1VecEpdGrp1ShiftEast,TVector2 Q1VecEpdGrp1ShiftWest, TVector2 Q1VecEpdGrp1ShiftFull); // Shift Corrected EPD Grp1 Q1Vector: East & West & Full
    void setTpcQFlag(int flagEp); // TPC Flag
    void setTpcQ1Vec(TVector2 Q1VecTpcReCtrEast,TVector2 Q1VecTpcReCtrWest); // ReCenter Corrected TPC Q1Vector: East & West
    void setTpcQ2Vec(TVector2 Q2VecTpcReCtrEast,TVector2 Q2VecTpcReCtrWest); // ReCenter Corrected TPC Q2Vector: East & West
    void setTpcQ3Vec(TVector2 Q3VecTpcReCtrEast,TVector2 Q3VecTpcReCtrWest); // ReCenter Corrected TPC Q3Vector: East & West
    void setNumTrks(int numTrkEast, int numTrkWest); // Number of Tracks used in ReCenter Correction: East & West

  private:
    // Isobar
    static const int mNumCentrality = 9;  // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    static const int mNumMixVzBin   = 10; // 10 vz bins for Isobar event mixing
    static const int mNumMixPsiBin  = 5;  // 5 TPC Psi2 bins for Isobar event mixing
    static const int mNumMixBuffer  = 5;  // 5 events for event mixing
    // FXT
    // static const int mNumCentrality = 9;  // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    // static const int mNumMixVzBin   = 2;  // 2 vz bins for FXT event mixing
    // static const int mNumMixPsiBin  = 10; // 10 EPD Psi1 bins for FXT event mixing
    // static const int mNumMixBuffer  = 5;  // 5 events for event mixing

    StAnalysisUtils *mAnaUtils; // Analysis Utilities
    StAnalysisCut   *mAnaCut;   // Analysis Cuts

    TTree *t_mPhiMesonTree;
    StPhiMesonEvent *mPhiMesonEvent;
    StPhiMesonTrack *mPhiMesonTrack;

    // event information container 
    // 0 = centrality bin, 1 = vertexZ bin, 2 = Psi bin || push_back->event
    int mEventCounter[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // General Event Info
    std::vector<int> vec_mRunId[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<int> vec_mRunIdx[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<int> vec_mEvtId[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<int> vec_mRefMult[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<int> vec_mNumTofMatch[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<int> vec_mCent9[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<int> vec_mCent16[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<int> vec_mRefWgt[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<double> vec_mZDCx[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<double> vec_mBBCx[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<double> vec_mVzVpd[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector3> vec_mPrimVtx[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];

    std::vector<int> vec_mFlagZdcEp[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // ZDC EP Info
    std::vector<TVector2> vec_mQ1ZdcShiftEast[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ1ZdcShiftWest[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ1ZdcShiftFull[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // shift corrected Full EP

    std::vector<int> vec_mFlagEpdSideEp[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // EPD EP Info
    std::vector<TVector2> vec_mQ1EpdSideShiftEast[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ1EpdSideShiftWest[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ1EpdSideShiftFull[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // shift corrected Full EP

    std::vector<int> vec_mFlagEpdGrp0Ep[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // EPD EP Info
    std::vector<TVector2> vec_mQ1EpdGrp0ShiftEast[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ1EpdGrp0ShiftWest[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // not relavant for FXT
    std::vector<TVector2> vec_mQ1EpdGrp0ShiftFull[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // not relavant for FXT
    std::vector<int> vec_mFlagEpdGrp1Ep[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // EPD EP Info
    std::vector<TVector2> vec_mQ1EpdGrp1ShiftEast[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ1EpdGrp1ShiftWest[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // not relavant for FXT
    std::vector<TVector2> vec_mQ1EpdGrp1ShiftFull[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // not relavant for FXT

    std::vector<int> vec_mFlagTpcEp[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // TPC EP Info
    std::vector<TVector2> vec_mQ1TpcReCtrEast[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ1TpcReCtrWest[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ2TpcReCtrEast[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ2TpcReCtrWest[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ3TpcReCtrEast[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ3TpcReCtrWest[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<int> vec_mNumTrkReCtrEast[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<int> vec_mNumTrkReCtrWest[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];

    // track info container
    // 0 = centrality bin, 1 = vertexZ bin, 2 = Psi bin, 3 = mixed event bin || push_back->track
    std::map<int, std::vector<TVector3> > map_mMomVecKp; // K+
    std::map<int, std::vector<TVector3> > map_mMomVecKm; // K-
    std::map<int, std::vector<double> > map_mMass2Kp;
    std::map<int, std::vector<double> > map_mMass2Km;
    std::map<int, std::vector<double> > map_mBetaKp;
    std::map<int, std::vector<double> > map_mBetaKm;
    std::map<int, std::vector<double> > map_mNSigKp;
    std::map<int, std::vector<double> > map_mNSigKm;
    std::map<int, std::vector<double> > map_mDcaKp;
    std::map<int, std::vector<double> > map_mDcaKm;
    std::map<int, std::vector<double> > map_mChargeKp;
    std::map<int, std::vector<double> > map_mChargeKm;
    std::map<int, std::vector<double> > map_mNHitsFitKp;
    std::map<int, std::vector<double> > map_mNHitsFitKm;

    TH2F *h_mInvMassPhi[mNumCentrality]; // pt vs. invMassPhi
    TH2F *h_mBetaKaon[mNumCentrality]; // p/q vs. 1/beta - 1/betaKaon
    TH2F *h_mMassKaon[mNumCentrality]; // m^2/q^2 vs. 1/beta - 1/betaKaon

    // set QVector
    int mFlagZdcEp; 
    TVector2 v_mQ1VecZdcShiftEast, v_mQ1VecZdcShiftWest, v_mQ1VecZdcShiftFull;

    int mFlagEpdSideEp; 
    TVector2 v_mQ1VecEpdSideShiftEast, v_mQ1VecEpdSideShiftWest, v_mQ1VecEpdSideShiftFull;

    int mFlagEpdGrp0Ep; 
    TVector2 v_mQ1VecEpdGrp0ShiftEast, v_mQ1VecEpdGrp0ShiftWest, v_mQ1VecEpdGrp0ShiftFull;
    int mFlagEpdGrp1Ep; 
    TVector2 v_mQ1VecEpdGrp1ShiftEast, v_mQ1VecEpdGrp1ShiftWest, v_mQ1VecEpdGrp1ShiftFull;

    int mFlagTpcEp;
    TVector2 v_mQ1VecTpcReCtrEast, v_mQ1VecTpcReCtrWest; 
    TVector2 v_mQ2VecTpcReCtrEast, v_mQ2VecTpcReCtrWest; 
    TVector2 v_mQ3VecTpcReCtrEast, v_mQ3VecTpcReCtrWest;
    int mNumTrkReCtrEast, mNumTrkReCtrWest;

    int mRunIdx, mCent9, mCent16;
    double mRefWgt, mVz, mPsiShiftFull;

    const int mType;

  ClassDef(StPhiMesonTree,1)
};
#endif
