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
    TVector2 getReCenterParEast();
    TVector2 getReCenterParWest();

    void print(TVector2);

    // Shift Correction
    void initTpcShift();
    void fillTpcShiftEast();
    void fillTpcShiftWest();
    void writeTpcShift();
    void readTpcShiftCorr();
    double getPsi2ShiftEast();
    double getPsi2ShiftWest();
    double angleShift(double PsiRaw);

    // Event Plane Resolution
    void initTpcResolution();
    void fillTpcResolution();
    void writeTpcResolution();
    void readTpcResolution();
    double getTpcSubEp2ResVal(int cent9);
    double getTpcSubEp2ResErr(int cent9);

    // Q2Vector
    TVector2 getQ2VecRawEast();
    TVector2 getQ2VecRawWest();
    TVector2 getQ2VecReCenterEast();
    TVector2 getQ2VecReCenterWest();
    double getPsi2RawEast();
    double getPsi2RawWest();
    double getPsi2ReCenterEast();
    double getPsi2ReCenterWest();
    int getNumTrkRawEast();
    int getNumTrkRawWest();
    int getNumTrkReCenterEast();
    int getNumTrkReCenterWest();

    // Event Plane Distribution
    void initTpcSubEpRaw(); // raw Sub EP
    void fillTpcSubEpRaw(double Psi2East, double Psi2West);
    void writeTpcSubEpRaw();

    void initTpcSubEpReCenter(); // recenter Sub EP
    void fillTpcSubEpReCenter(double Psi2East, double Psi2West);
    void writeTpcSubEpReCenter();

    void initTpcSubEpShift(); // shift Sub EP
    void fillTpcSubEpShift(double Psi2East, double Psi2West);
    void writeTpcSubEpShift();

  private:
    static const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
    static const int mNumShiftCorr = 20;

    int mCent9;
    int mRunIndex;
    int mVzBin;

    TVector2 v_mQ2RawEast, v_mQ2RawWest;
    TVector2 v_mQ2ReCenterEast, v_mQ2ReCenterWest; 

    int mQCouRawEast, mQCouRawWest;
    int mQCouReCenterEast, mQCouReCenterWest;

    // ReCenter Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mTpcQ2ReCenterXEast[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCenterYEast[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCenterXWest[mNumVzBin];
    TProfile2D *p_mTpcQ2ReCenterYWest[mNumVzBin];

    // Shift Correction for East/West
    TProfile2D *p_mTpcQ2ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mTpcQ2ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ2ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQ2ShiftSinWest[mNumVzBin][mNumShiftCorr];

    // TPC EP Resolution
    TProfile *p_mTpcSubEp2Res;
    double mTpcSubEp2ResVal[9];
    double mTpcSubEp2ResErr[9];

    // Event Plane Distribution
    TH2F *h_mTpcEp2RawEast[9]; // raw EP
    TH2F *h_mTpcEp2RawWest[9];
    TH2F *h_mTpcEp2RawCorr[9]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp2ReCenterEast[9]; // recenter EP
    TH2F *h_mTpcEp2ReCenterWest[9];
    TH2F *h_mTpcEp2ReCenterCorr[9]; // Psi2East vs Psi2West

    TH2F *h_mTpcEp2ShiftEast[9]; // shift EP
    TH2F *h_mTpcEp2ShiftWest[9];
    TH2F *h_mTpcEp2ShiftCorr[9]; // Psi2East vs Psi2West

    TFile *file_mReCenterPar;
    TFile *file_mShiftPar;
    TFile *file_mResolution;

    const int mType;

  ClassDef(StTpcEpManager,1)
};

#endif
