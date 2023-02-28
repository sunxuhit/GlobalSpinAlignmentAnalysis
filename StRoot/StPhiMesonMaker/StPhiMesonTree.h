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
class StPhiMesonEvent;
class StPhiMesonTrack;

class StPhiMesonTree : public TObject
{
  public:
    StPhiMesonTree(int beamType);
    virtual ~StPhiMesonTree();

    void initPhiTree();
    void clearPhiMixBuffer(int cent9, int vzBin, int PsiBin);

    int getPhiMixKey(int cent9, int vzBin, int PsiBin, int evtBin); // return 1000*cent9 + 100*vzBin + 10*PsiBin + evtBin
    int getVzMixBin(double vz);
    int getPsi2MixBin(double Psi2);
    void getPhiEvtSize(int,int,int);

    void recoPhi(int cent9, int vzBin, int Psi2Bin); // reconstruct phi meson in the same event
    void mixPhi(int cent9, int vzBin, int Psi2Bin); // reconstruct phi meson in the mixed event
    void fillPhiEvent(StPicoDst* picoDst, int flagME);

    void writePhiTree();

    void clearEvtInfo();
    void setEvtInfo(int cent9, int cent16, double refwgt, double vz, double Psi2ShiftFull);
    void setZdcQ1Vec(int flagEp, TVector2 Q1VecZdcShiftEast,TVector2 Q1VecZdcShiftWest, TVector2 Q1VecZdcShiftFull); // Shift Corrected ZDC Q1Vector: East & West & Full
    void setEpdQ1Vec(int flagEp, TVector2 Q1VecEpdShiftEast,TVector2 Q1VecEpdShiftWest, TVector2 Q1VecEpdShiftFull); // Shift Corrected EPD Q1Vector: East & West & Full
    void setTpcQVec(int flagEp, TVector2 Q2VecTpcReCtrEast,TVector2 Q2VecTpcReCtrWest, TVector2 Q3VecTpcReCtrEast,TVector2 Q3VecTpcReCtrWest); // ReCenter Corrected TPC Q2Vector & Q3Vector: East & West
    void setNumTrks(int,int); // Number of Tracks used in ReCenter Correction: East & West

  private:
    static const int mNumCentrality = 9;  // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    static const int mNumMixVzBin   = 10; // 10 vz bins for event mixing
    static const int mNumMixPsiBin  = 5;  // 5 TPC Psi2 bins for event mixing
    static const int mNumMixBuffer  = 5;  // 5 events for event mixing

    TTree *t_mPhiMesonTree;
    StPhiMesonEvent *mPhiMesonEvent;
    StPhiMesonTrack *mPhiMesonTrack;

    // event information container 
    // 0 = centrality bin, 1 = vertexZ bin, 2 = Psi bin || push_back->event
    int mEventCounter[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // General Event Info
    std::vector<int> vec_mRunId[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
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
    std::vector<TVector2> vec_mQ1ZdcShiftFull[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];

    std::vector<int> vec_mFlagEpdEp[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // EPD EP Info
    std::vector<TVector2> vec_mQ1EpdShiftEast[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ1EpdShiftWest[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];
    std::vector<TVector2> vec_mQ1EpdShiftFull[mNumCentrality][mNumMixVzBin][mNumMixPsiBin];

    std::vector<int> vec_mFlagTpcEp[mNumCentrality][mNumMixVzBin][mNumMixPsiBin]; // TPC EP Info
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
    std::map<int, std::vector<double> > map_mNSigKp;
    std::map<int, std::vector<double> > map_mNSigKm;
    std::map<int, std::vector<double> > map_mDcaKp;
    std::map<int, std::vector<double> > map_mDcaKm;
    std::map<int, std::vector<double> > map_mChargeKp;
    std::map<int, std::vector<double> > map_mChargeKm;
    std::map<int, std::vector<double> > map_mNHitsFitKp;
    std::map<int, std::vector<double> > map_mNHitsFitKm;

    TH2F *h_mPhiMass2[mNumCentrality];

    // set QVector
    int mFlagZdcEp, mFlagEpdEp, mFlagTpcEp;
    TVector2 v_mQ1VecZdcShiftEast, v_mQ1VecZdcShiftWest, v_mQ1VecZdcShiftFull;
    TVector2 v_mQ1VecEpdShiftEast, v_mQ1VecEpdShiftWest, v_mQ1VecEpdShiftFull;
    TVector2 v_mQ2VecTpcReCtrEast, v_mQ2VecTpcReCtrWest, v_mQ3VecTpcReCtrEast, v_mQ3VecTpcReCtrWest;
    int mNumTrkReCtrEast, mNumTrkReCtrWest;

    int mCent9, mCent16;
    double mRefWgt, mVz, mPsi2ShiftFull;
    const int mType;

  ClassDef(StPhiMesonTree,1)
};
#endif
