#ifndef StTpcEpManager_h
#define StTpcEpManager_h

#include <string>
#include "TObject.h"
#include "TVector2.h"

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
    TVector2 calq2Vector(StPicoTrack* picoTrack);
    TVector2 calq3Vector(StPicoTrack* picoTrack);
    double getWeight(StPicoTrack* picoTrack);
    void addTrackRawEast(StPicoTrack* picoTrack);
    void addTrackRawWest(StPicoTrack* picoTrack);
    void addTrackReCtrEast(StPicoTrack* picoTrack);
    void addTrackReCtrWest(StPicoTrack* picoTrack);

    // ReCenter Correction
    void initTpcReCtr();
    void fillTpcReCtrEast(StPicoTrack* picoTrack);
    void fillTpcReCtrWest(StPicoTrack* picoTrack);
    void writeTpcReCtr();
    void readTpcReCtr();
    TVector2 getq2VecCtrEast(); // 2nd ReCenter Parameter
    TVector2 getq2VecCtrWest();
    TVector2 getq3VecCtrEast(); // 3rd ReCenter Parameter
    TVector2 getq3VecCtrWest();

    // Shift Correction
    void initTpcShift();
    void fillTpcShiftEast();
    void fillTpcShiftWest();
    void writeTpcShift();
    void readTpcShift();
    double getPsi2ShiftEast(TVector2 Q2Vector); // 2nd shift Psi2
    double getPsi2ShiftWest(TVector2 Q2Vector);
    double transPsi2(double Psi2);
    double getPsi3ShiftEast(TVector2 Q3Vector); // 3rd shift Psi3
    double getPsi3ShiftWest(TVector2 Q3Vector);
    double transPsi3(double Psi3);

    // Event Plane Resolution
    void initTpcResolution();
    void fillTpcResolution(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeTpcResolution();
    void readTpcResolution();
    double getTpcSubEp2ResVal(int cent9);
    double getTpcSubEp2ResErr(int cent9);
    double getTpcSubEp3ResVal(int cent9);
    double getTpcSubEp3ResErr(int cent9);

    void initTpcSubEpFlow(); // Sub EP
    void fillTpcSubEpEFlow(double pt, double v2, double reweight);
    void fillTpcSubEpTFlow(double pt, double v3, double reweight);
    void writeTpcSubEpFlow();

    // Q2Vector
    TVector2 getQ2VecRawEast(); // Q2Vector
    TVector2 getQ2VecRawWest();
    TVector2 getQ2VecReCtrEast();
    TVector2 getQ2VecReCtrWest();

    double getPsi2RawEast();
    double getPsi2RawWest();
    double getPsi2ReCtrEast();
    double getPsi2ReCtrWest();

    TVector2 getQ3VecRawEast(); // Q3Vector
    TVector2 getQ3VecRawWest();
    TVector2 getQ3VecReCtrEast();
    TVector2 getQ3VecReCtrWest();

    double getPsi3RawEast();
    double getPsi3RawWest();
    double getPsi3ReCtrEast();
    double getPsi3ReCtrWest();

    int getNumTrkRawEast(); // num of Tracks used in EP reconstruction
    int getNumTrkRawWest();
    int getNumTrkReCtrEast();
    int getNumTrkReCtrWest();

    // Event Plane Distribution
    void initTpcSubEpRaw(); // raw Sub EP
    void fillTpcSubEpRaw(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeTpcSubEpRaw();

    void initTpcSubEpReCtr(); // recenter Sub EP
    void fillTpcSubEpReCtr(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeTpcSubEpReCtr();

    void initTpcSubEpShift(); // shift Sub EP
    void fillTpcSubEpShift(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeTpcSubEpShift();

  private:
    static const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
    static const int mNumShiftCorr = 20;
    static const int mNumCentrality = 9; // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%

    int mCent9;
    int mRunIndex;
    int mVzBin;

    int mQCouRawEast, mQCouRawWest;
    int mQCouReCtrEast, mQCouReCtrWest;

    TVector2 v_mQ2RawEast, v_mQ2RawWest; // 2nd EP
    TVector2 v_mQ2ReCtrEast, v_mQ2ReCtrWest; 

    TVector2 v_mQ3RawEast, v_mQ3RawWest; // 3rd EP
    TVector2 v_mQ3ReCtrEast, v_mQ3ReCtrWest; 

    // ReCenter Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mTpcQ2ReCtrXEast[mNumVzBin]; // 2nd EP
    TProfile2D *p_mTpcQ2ReCtrYEast[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCtrXWest[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCtrYWest[mNumVzBin];

    TProfile2D *p_mTpcQ3ReCtrXEast[mNumVzBin]; // 3rd EP
    TProfile2D *p_mTpcQ3ReCtrYEast[mNumVzBin];
    TProfile2D *p_mTpcQ3ReCtrXWest[mNumVzBin];
    TProfile2D *p_mTpcQ3ReCtrYWest[mNumVzBin];

    // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mTpcQ2ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 2nd EP
    TProfile2D *p_mTpcQ2ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ2ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ2ShiftSinWest[mNumVzBin][mNumShiftCorr];

    TProfile2D *p_mTpcQ3ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 3rd EP
    TProfile2D *p_mTpcQ3ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ3ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ3ShiftSinWest[mNumVzBin][mNumShiftCorr];

    // TPC EP Resolution
    TProfile *p_mTpcSubEp2Res; // 2nd EP
    double mTpcSubEp2ResVal[mNumCentrality];
    double mTpcSubEp2ResErr[mNumCentrality];

    TProfile *p_mTpcSubEp3Res; // 3rd EP
    double mTpcSubEp3ResVal[mNumCentrality];
    double mTpcSubEp3ResErr[mNumCentrality];

    // Charged Hadron Elliptic and Triangular Flow
    TProfile *p_mTpcSubEpEFlow[mNumCentrality]; // v2 vs. pT
    TProfile *p_mTpcSubEpTFlow[mNumCentrality]; // v3 vs. pT

    // Event Plane Distribution
    TH2F *h_mTpcEp2RawEast[mNumCentrality]; // 2nd raw EP
    TH2F *h_mTpcEp2RawWest[mNumCentrality];
    TH2F *h_mTpcEp2RawCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp3RawEast[mNumCentrality]; // 3rd raw EP
    TH2F *h_mTpcEp3RawWest[mNumCentrality];
    TH2F *h_mTpcEp3RawCorr[mNumCentrality]; // Psi3East vs Psi3West

    TH2F *h_mTpcEp2ReCtrEast[mNumCentrality]; // 2nd recenter EP
    TH2F *h_mTpcEp2ReCtrWest[mNumCentrality];
    TH2F *h_mTpcEp2ReCtrCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp3ReCtrEast[mNumCentrality]; // 3rd recenter EP
    TH2F *h_mTpcEp3ReCtrWest[mNumCentrality];
    TH2F *h_mTpcEp3ReCtrCorr[mNumCentrality]; // Psi3East vs Psi3West

    TH2F *h_mTpcEp2ShiftEast[mNumCentrality]; // 3rd shift EP
    TH2F *h_mTpcEp2ShiftWest[mNumCentrality];
    TH2F *h_mTpcEp2ShiftCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp3ShiftEast[mNumCentrality]; // 3rd shift EP
    TH2F *h_mTpcEp3ShiftWest[mNumCentrality];
    TH2F *h_mTpcEp3ShiftCorr[mNumCentrality]; // Psi3East vs Psi3West

    TFile *file_mReCtrPar;
    TFile *file_mShiftPar;
    TFile *file_mResolution;

    const int mType;

  ClassDef(StTpcEpManager,1)
};

#endif
