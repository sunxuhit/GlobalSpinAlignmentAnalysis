#ifndef StTpcEpManager_h
#define StTpcEpManager_h

#include <string>
#include "TObject.h"
#include "TVector2.h"
#include "TVector3.h"

class StPicoTrack;
class TFile;
class TProfile2D;
class TProfile;
class TH2F;

class StTpcEpManager : public TObject
{
  public:
    StTpcEpManager(int beamType);
    virtual ~StTpcEpManager();
    void clearTpcEpManager();
    void initTpcEpManager(int cent9, int runIndex, int vzBin);

    // Utilities
    TVector2 calq1Vector(StPicoTrack* picoTrack);
    TVector2 calq2Vector(StPicoTrack* picoTrack);
    TVector2 calq3Vector(StPicoTrack* picoTrack);
    double getWeight(StPicoTrack* picoTrack);
    void addTrackRawEast(StPicoTrack* picoTrack);
    void addTrackRawWest(StPicoTrack* picoTrack);
    void addTrackRawFull(StPicoTrack* picoTrack);
    void addTrackReCtrEast(StPicoTrack* picoTrack);
    void addTrackReCtrWest(StPicoTrack* picoTrack);
    void addTrackReCtrFull(StPicoTrack* picoTrack);
    TVector2 calq1Vector(TVector3 primMom); // used in StPhiMesonAnalyzer
    TVector2 calq2Vector(TVector3 primMom);
    TVector2 calq3Vector(TVector3 primMom);
    double getWeight(TVector3 primMom);

    // ReCenter Correction
    void initTpcReCtr();
    void fillTpcReCtrEast(StPicoTrack* picoTrack);
    void fillTpcReCtrWest(StPicoTrack* picoTrack);
    void fillTpcReCtrFull(StPicoTrack* picoTrack);
    void writeTpcReCtr();
    void readTpcReCtr();
    TVector2 getq1VecCtrEast(); // 1st ReCtr Parameter
    TVector2 getq1VecCtrWest();
    TVector2 getq1VecCtrFull();
    TVector2 getq2VecCtrEast(); // 2nd ReCtr Parameter
    TVector2 getq2VecCtrWest();
    TVector2 getq2VecCtrFull();
    TVector2 getq3VecCtrEast(); // 3rd ReCtr Parameter
    TVector2 getq3VecCtrWest();
    TVector2 getq3VecCtrFull();

    // Shift Correction
    void initTpcShift();
    void fillTpcShiftEast();
    void fillTpcShiftWest();
    void fillTpcShiftFull();
    void writeTpcShift();
    void readTpcShift();
    double getPsi1ShiftEast(TVector2 Q1Vector); // 1st shift angle
    double getPsi1ShiftWest(TVector2 Q1Vector);
    double getPsi1ShiftFull(TVector2 Q1Vector);
    double transPsi1(double Psi1);
    bool isPsi1InRange(double Psi1);
    double getPsi2ShiftEast(TVector2 Q2Vector); // 2nd shift angle
    double getPsi2ShiftWest(TVector2 Q2Vector);
    double getPsi2ShiftFull(TVector2 Q2Vector);
    double transPsi2(double Psi2);
    bool isPsi2InRange(double Psi2);
    double getPsi3ShiftEast(TVector2 Q3Vector); // 3rd shift angle
    double getPsi3ShiftWest(TVector2 Q3Vector);
    double getPsi3ShiftFull(TVector2 Q3Vector);
    double transPsi3(double Psi3);
    bool isPsi3InRange(double Psi3);

    // Event Plane Resolution
    void initTpcResolution();
    void fillTpcResolution(double Psi1East, double Psi1West, double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeTpcResolution();
    void readTpcResolution();
    double getTpcSubEp1ResVal(int cent9);
    double getTpcSubEp1ResErr(int cent9);
    double getTpcSubEp2ResVal(int cent9);
    double getTpcSubEp2ResErr(int cent9);
    double getTpcSubEp3ResVal(int cent9);
    double getTpcSubEp3ResErr(int cent9);

    void initTpcSubEpFlow(); // Sub EP
    void fillTpcSubEpV1(double pt, double v1, double reweight);
    void fillTpcSubEpV2(double pt, double v2, double reweight);
    void fillTpcSubEpV3(double pt, double v3, double reweight);
    void writeTpcSubEpFlow();

    // QVector
    TVector2 getQ1VecRawEast(); // Q1Vector
    TVector2 getQ1VecRawWest();
    TVector2 getQ1VecRawFull();
    TVector2 getQ1VecReCtrEast();
    TVector2 getQ1VecReCtrWest();
    TVector2 getQ1VecReCtrFull();
    double getPsi1RawEast();
    double getPsi1RawWest();
    double getPsi1RawFull();
    double getPsi1ReCtrEast();
    double getPsi1ReCtrWest();
    double getPsi1ReCtrFull();

    TVector2 getQ2VecRawEast(); // Q2Vector
    TVector2 getQ2VecRawWest();
    TVector2 getQ2VecRawFull();
    TVector2 getQ2VecReCtrEast();
    TVector2 getQ2VecReCtrWest();
    TVector2 getQ2VecReCtrFull();
    double getPsi2RawEast();
    double getPsi2RawWest();
    double getPsi2RawFull();
    double getPsi2ReCtrEast();
    double getPsi2ReCtrWest();
    double getPsi2ReCtrFull();

    TVector2 getQ3VecRawEast(); // Q3Vector
    TVector2 getQ3VecRawWest();
    TVector2 getQ3VecRawFull();
    TVector2 getQ3VecReCtrEast();
    TVector2 getQ3VecReCtrWest();
    TVector2 getQ3VecReCtrFull();
    double getPsi3RawEast();
    double getPsi3RawWest();
    double getPsi3RawFull();
    double getPsi3ReCtrEast();
    double getPsi3ReCtrWest();
    double getPsi3ReCtrFull();

    int getNumTrkRawEast(); // num of Tracks used in EP reconstruction
    int getNumTrkRawWest();
    int getNumTrkRawFull();
    int getNumTrkReCtrEast();
    int getNumTrkReCtrWest();
    int getNumTrkReCtrFull();

    // Event Plane Distribution
    void initTpcSubEpRaw(); // raw Sub EP
    void fillTpcSubEpRaw(double Psi1East, double Psi1West, double Psi1Full, double Psi2East, double Psi2West, double Psi2Full, double Psi3East, double Psi3West, double Psi3Full);
    void writeTpcSubEpRaw();

    void initTpcSubEpReCtr(); // recenter Sub EP
    void fillTpcSubEpReCtr(double Psi1East, double Psi1West, double Psi1Full, double Psi2East, double Psi2West, double Psi2Full, double Psi3East, double Psi3West, double Psi3Full);
    void writeTpcSubEpReCtr();

    void initTpcSubEpShift(); // shift Sub EP
    void fillTpcSubEpShift(double Psi1East, double Psi1West, double Psi1Full, double Psi2East, double Psi2West, double Psi2Full, double Psi3East, double Psi3West, double Psi3Full);
    void writeTpcSubEpShift();

  private:
    static const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
    static const int mNumShiftCorr = 20;
    static const int mNumCentrality = 9; // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%

    int mCent9;
    int mRunIndex;
    int mVzBin;

    int mQCouRawEast, mQCouRawWest, mQCouRawFull;
    int mQCouReCtrEast, mQCouReCtrWest, mQCouReCtrFull;

    TVector2 v_mQ1RawEast, v_mQ1RawWest, v_mQ1RawFull; // 1st EP
    TVector2 v_mQ1ReCtrEast, v_mQ1ReCtrWest, v_mQ1ReCtrFull; 

    TVector2 v_mQ2RawEast, v_mQ2RawWest, v_mQ2RawFull; // 2nd EP
    TVector2 v_mQ2ReCtrEast, v_mQ2ReCtrWest, v_mQ2ReCtrFull; 

    TVector2 v_mQ3RawEast, v_mQ3RawWest, v_mQ3RawFull; // 3rd EP
    TVector2 v_mQ3ReCtrEast, v_mQ3ReCtrWest, v_mQ3ReCtrFull; 

    // ReCenter Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mTpcQ1ReCtrXEast[mNumVzBin]; // 1st EP
    TProfile2D *p_mTpcQ1ReCtrYEast[mNumVzBin];
    TProfile2D *p_mTpcQ1ReCtrXWest[mNumVzBin];
    TProfile2D *p_mTpcQ1ReCtrYWest[mNumVzBin];
    TProfile2D *p_mTpcQ1ReCtrXFull[mNumVzBin];
    TProfile2D *p_mTpcQ1ReCtrYFull[mNumVzBin];

    TProfile2D *p_mTpcQ2ReCtrXEast[mNumVzBin]; // 2nd EP
    TProfile2D *p_mTpcQ2ReCtrYEast[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCtrXWest[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCtrYWest[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCtrXFull[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCtrYFull[mNumVzBin];

    TProfile2D *p_mTpcQ3ReCtrXEast[mNumVzBin]; // 3rd EP
    TProfile2D *p_mTpcQ3ReCtrYEast[mNumVzBin];
    TProfile2D *p_mTpcQ3ReCtrXWest[mNumVzBin];
    TProfile2D *p_mTpcQ3ReCtrYWest[mNumVzBin];
    TProfile2D *p_mTpcQ3ReCtrXFull[mNumVzBin];
    TProfile2D *p_mTpcQ3ReCtrYFull[mNumVzBin];

    // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mTpcQ1ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 1st EP
    TProfile2D *p_mTpcQ1ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ1ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ1ShiftSinWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ1ShiftCosFull[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ1ShiftSinFull[mNumVzBin][mNumShiftCorr];

    TProfile2D *p_mTpcQ2ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 2nd EP
    TProfile2D *p_mTpcQ2ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ2ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ2ShiftSinWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ2ShiftCosFull[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ2ShiftSinFull[mNumVzBin][mNumShiftCorr];

    TProfile2D *p_mTpcQ3ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 3rd EP
    TProfile2D *p_mTpcQ3ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ3ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ3ShiftSinWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ3ShiftCosFull[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ3ShiftSinFull[mNumVzBin][mNumShiftCorr];

    // TPC EP Resolution
    TProfile *p_mTpcSubEp1Res; // 1st EP
    double mTpcSubEp1ResVal[mNumCentrality];
    double mTpcSubEp1ResErr[mNumCentrality];

    TProfile *p_mTpcSubEp2Res; // 2nd EP
    double mTpcSubEp2ResVal[mNumCentrality];
    double mTpcSubEp2ResErr[mNumCentrality];

    TProfile *p_mTpcSubEp3Res; // 3rd EP
    double mTpcSubEp3ResVal[mNumCentrality];
    double mTpcSubEp3ResErr[mNumCentrality];

    // Charged Hadron Elliptic and Triangular Flow
    TProfile *p_mTpcSubEpV1[mNumCentrality]; // v1 vs. pT
    TProfile *p_mTpcSubEpV2[mNumCentrality]; // v2 vs. pT
    TProfile *p_mTpcSubEpV3[mNumCentrality]; // v3 vs. pT

    // Event Plane Distribution
    TH2F *h_mTpcEp1RawEast[mNumCentrality]; // 1st raw EP
    TH2F *h_mTpcEp1RawWest[mNumCentrality];
    TH2F *h_mTpcEp1RawFull[mNumCentrality];
    TH2F *h_mTpcEp1RawCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mTpcEp2RawEast[mNumCentrality]; // 2nd raw EP
    TH2F *h_mTpcEp2RawWest[mNumCentrality];
    TH2F *h_mTpcEp2RawFull[mNumCentrality];
    TH2F *h_mTpcEp2RawCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp3RawEast[mNumCentrality]; // 3rd raw EP
    TH2F *h_mTpcEp3RawWest[mNumCentrality];
    TH2F *h_mTpcEp3RawFull[mNumCentrality];
    TH2F *h_mTpcEp3RawCorr[mNumCentrality]; // Psi3East vs Psi3West

    TH2F *h_mTpcEp1ReCtrEast[mNumCentrality]; // 1st recenter EP
    TH2F *h_mTpcEp1ReCtrWest[mNumCentrality];
    TH2F *h_mTpcEp1ReCtrFull[mNumCentrality];
    TH2F *h_mTpcEp1ReCtrCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mTpcEp2ReCtrEast[mNumCentrality]; // 2nd recenter EP
    TH2F *h_mTpcEp2ReCtrWest[mNumCentrality];
    TH2F *h_mTpcEp2ReCtrFull[mNumCentrality];
    TH2F *h_mTpcEp2ReCtrCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp3ReCtrEast[mNumCentrality]; // 3rd recenter EP
    TH2F *h_mTpcEp3ReCtrWest[mNumCentrality];
    TH2F *h_mTpcEp3ReCtrFull[mNumCentrality];
    TH2F *h_mTpcEp3ReCtrCorr[mNumCentrality]; // Psi3East vs Psi3West

    TH2F *h_mTpcEp1ShiftEast[mNumCentrality]; // 1st shift EP
    TH2F *h_mTpcEp1ShiftWest[mNumCentrality];
    TH2F *h_mTpcEp1ShiftFull[mNumCentrality];
    TH2F *h_mTpcEp1ShiftCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mTpcEp2ShiftEast[mNumCentrality]; // 2nd shift EP
    TH2F *h_mTpcEp2ShiftWest[mNumCentrality];
    TH2F *h_mTpcEp2ShiftFull[mNumCentrality];
    TH2F *h_mTpcEp2ShiftCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp3ShiftEast[mNumCentrality]; // 3rd shift EP
    TH2F *h_mTpcEp3ShiftWest[mNumCentrality];
    TH2F *h_mTpcEp3ShiftFull[mNumCentrality];
    TH2F *h_mTpcEp3ShiftCorr[mNumCentrality]; // Psi3East vs Psi3West

    TFile *file_mReCtrPar;
    TFile *file_mShiftPar;
    TFile *file_mResolution;

    const int mType;

  ClassDef(StTpcEpManager,1)
};

#endif
