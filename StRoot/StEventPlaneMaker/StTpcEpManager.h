#ifndef StTpcEpManager_h
#define StTpcEpManager_h

#include <string>
#include "TObject.h"
#include "TVector2.h"
// #include "TString.h"

class StPicoDst;
class StPicoTrack;
class TProfile2D;
class TFile;
// class TNtuple;

class StTpcEpManager : public TObject
{
  public:
    StTpcEpManager(int beamType);
    virtual ~StTpcEpManager();
    void clearTpcEp();
    void initTpcEp(int cent9, int runIndex, int vzBin);

    TVector2 calq2Vector(StPicoTrack* picoTrack);
    double getWeight(StPicoTrack* picoTrack);

    void readReCenterCorr();

    void addTrackEastRaw(StPicoTrack* picoTrack);
    void addTrackEast(StPicoTrack* picoTrack);
    TVector2 getReCenterParEast();

    void addTrackWestRaw(StPicoTrack* picoTrack);
    void addTrackWest(StPicoTrack* picoTrack);
    TVector2 getReCenterParWest();

    void print(TVector2);


    // Shift Correction
    // TVector2 calPsi2_East_EP(int); // 0 = ShiftOrder: 2, 4, 6, 8, 10
    // TVector2 calPsi2_West_EP(int);

    void readShiftCorr();
    double angleShift(double PsiRaw);

    // Event Plane method
    double calShiftAngle2East();
    double calShiftAngle2West();

    void readResolution();
    double getResolutionVal(int cent9);
    double getResolutionErr(int cent9);

    TVector2 getQVector(int epMode); // east/west
    TVector2 getQVectorRaw(int epMode);
    int getNumTrack(int epMode);

  private:
    static const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
    static const int mNumShiftCorr = 20;

    int mCent9;
    int mRunIndex;
    int mVzBin;

    TVector2 v_mQ2EastRaw, v_mQ2WestRaw;
    TVector2 v_mQ2EastReCenter, v_mQ2WestReCenter; 

    int mQCouEastRaw, mQCouWestRaw;
    int mQCouEastReCenter, mQCouWestReCenter;

    // ReCenter Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mTpcQReCenterXEast[mNumVzBin];
    TProfile2D *p_mTpcQReCenterYEast[mNumVzBin];
    TProfile2D *p_mTpcQReCenterXWest[mNumVzBin];
    TProfile2D *p_mTpcQReCenterYWest[mNumVzBin];

    // Shift Correction for East/West
    TProfile2D *p_mTpcQShiftCosEast[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mTpcQShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mZdcQShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mZdcQShiftSinWest[mNumVzBin][mNumShiftCorr];

    // TPC EP Resolution
    TProfile *p_mTpcEpResolution;
    double mResolutionVal[9];
    double mResolutionErr[9];

    TFile *file_mReCenterPar;
    TFile *file_mShiftPar;
    TFile *file_mResolution;

    const int mType;

    std::string mVStr[2] = {"pos","neg"};
    std::string mOrder = "2nd";

  ClassDef(StTpcEpManager,1)
};

#endif
