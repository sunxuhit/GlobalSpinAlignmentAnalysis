#ifndef StTpcEpManager_h
#define StTpcEpManager_h

#include <string>
#include "TObject.h"
#include "TVector2.h"

class StPicoTrack;
class TProfile2D;
class TProfile;
class TFile;

class StTpcEpManager : public TObject
{
  public:
    StTpcEpManager(int beamType);
    virtual ~StTpcEpManager();
    void clearTpcEp();
    void initTpcEpManager(int cent9, int runIndex, int vzBin);

    // Utilities
    TVector2 calq2Vector(StPicoTrack* picoTrack);
    TVector2 calq3Vector(StPicoTrack* picoTrack);
    double getWeight(StPicoTrack* picoTrack);
    void addTrackRawEast(StPicoTrack* picoTrack);
    void addTrackRawWest(StPicoTrack* picoTrack);
    void addTrackReCenterEast(StPicoTrack* picoTrack);
    void addTrackReCenterWest(StPicoTrack* picoTrack);

    // ReCenter Correction
    void initTpcReCenter();
    void fillTpcReCenterEast(StPicoTrack* picoTrack);
    void fillTpcReCenterWest(StPicoTrack* picoTrack);
    void writeTpcReCenter();
    void readTpcReCenterCorr();
    TVector2 getq2VecReCenterEast(); // 2nd ReCenter Parameter
    TVector2 getq2VecReCenterWest();
    TVector2 getq3VecReCenterEast(); // 3rd ReCenter Parameter
    TVector2 getq3VecReCenterWest();

    // Shift Correction
    void initTpcShift();
    void fillTpcShiftEast();
    void fillTpcShiftWest();
    void writeTpcShift();
    void readTpcShiftCorr();
    double getPsi2ShiftEast(); // 2nd shift Psi2
    double getPsi2ShiftWest();
    double transPsi2(double Psi2);
    double getPsi3ShiftEast(); // 3rd shift Psi3
    double getPsi3ShiftWest();
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

    // Q2Vector
    TVector2 getQ2VecRawEast(); // Q2Vector
    TVector2 getQ2VecRawWest();
    TVector2 getQ2VecReCenterEast();
    TVector2 getQ2VecReCenterWest();
    double getPsi2RawEast();
    double getPsi2RawWest();
    double getPsi2ReCenterEast();
    double getPsi2ReCenterWest();

    TVector2 getQ3VecRawEast(); // Q3Vector
    TVector2 getQ3VecRawWest();
    TVector2 getQ3VecReCenterEast();
    TVector2 getQ3VecReCenterWest();
    double getPsi3RawEast();
    double getPsi3RawWest();
    double getPsi3ReCenterEast();
    double getPsi3ReCenterWest();

    int getNumTrkRawEast();
    int getNumTrkRawWest();
    int getNumTrkReCenterEast();
    int getNumTrkReCenterWest();

    // Event Plane Distribution
    void initTpcSubEpRaw(); // raw Sub EP
    void fillTpcSubEpRaw(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeTpcSubEpRaw();

    void initTpcSubEpReCenter(); // recenter Sub EP
    void fillTpcSubEpReCenter(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeTpcSubEpReCenter();

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
    int mQCouReCenterEast, mQCouReCenterWest;

    TVector2 v_mQ2RawEast, v_mQ2RawWest; // 2nd EP
    TVector2 v_mQ2ReCenterEast, v_mQ2ReCenterWest; 

    TVector2 v_mQ3RawEast, v_mQ3RawWest; // 3rd EP
    TVector2 v_mQ3ReCenterEast, v_mQ3ReCenterWest; 

    // ReCenter Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mTpcQ2ReCenterXEast[mNumVzBin]; // 2nd EP
    TProfile2D *p_mTpcQ2ReCenterYEast[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCenterXWest[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCenterYWest[mNumVzBin];

    TProfile2D *p_mTpcQ3ReCenterXEast[mNumVzBin]; // 3rd EP
    TProfile2D *p_mTpcQ3ReCenterYEast[mNumVzBin];
    TProfile2D *p_mTpcQ3ReCenterXWest[mNumVzBin];
    TProfile2D *p_mTpcQ3ReCenterYWest[mNumVzBin];

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

    // Event Plane Distribution
    TH2F *h_mTpcEp2RawEast[mNumCentrality]; // 2nd raw EP
    TH2F *h_mTpcEp2RawWest[mNumCentrality];
    TH2F *h_mTpcEp2RawCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp3RawEast[mNumCentrality]; // 3rd raw EP
    TH2F *h_mTpcEp3RawWest[mNumCentrality];
    TH2F *h_mTpcEp3RawCorr[mNumCentrality]; // Psi3East vs Psi3West

    TH2F *h_mTpcEp2ReCenterEast[mNumCentrality]; // 2nd recenter EP
    TH2F *h_mTpcEp2ReCenterWest[mNumCentrality];
    TH2F *h_mTpcEp2ReCenterCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp3ReCenterEast[mNumCentrality]; // 3rd recenter EP
    TH2F *h_mTpcEp3ReCenterWest[mNumCentrality];
    TH2F *h_mTpcEp3ReCenterCorr[mNumCentrality]; // Psi3East vs Psi3West

    TH2F *h_mTpcEp2ShiftEast[mNumCentrality]; // 3rd shift EP
    TH2F *h_mTpcEp2ShiftWest[mNumCentrality];
    TH2F *h_mTpcEp2ShiftCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp3ShiftEast[mNumCentrality]; // 3rd shift EP
    TH2F *h_mTpcEp3ShiftWest[mNumCentrality];
    TH2F *h_mTpcEp3ShiftCorr[mNumCentrality]; // Psi3East vs Psi3West

    TFile *file_mReCenterPar;
    TFile *file_mShiftPar;
    TFile *file_mResolution;

    const int mType;

  ClassDef(StTpcEpManager,1)
};

#endif
