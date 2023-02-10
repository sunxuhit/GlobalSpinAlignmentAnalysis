#ifndef StTpcEpManager_h
#define StTpcEpManager_h

#include <string>
#include "TObject.h"
#include "TVector2.h"

class StPicoTrack;
class TProfile2D;
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
    void addTrackEastRaw(StPicoTrack* picoTrack);
    void addTrackWestRaw(StPicoTrack* picoTrack);
    void addTrackEastReCenter(StPicoTrack* picoTrack);
    void addTrackWestReCenter(StPicoTrack* picoTrack);

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
    double angleShift(double PsiRaw);

    // Event Plane method
    double calShiftAngle2East();
    double calShiftAngle2West();

    void readTpcResolution();
    double getTpcResSubVal(int cent9);
    double getTpcResSubErr(int cent9);

    TVector2 getQVectorEastRaw(); // east/west
    TVector2 getQVectorWestRaw(); // east/west
    TVector2 getQVectorEastReCenter();
    TVector2 getQVectorWestReCenter();
    int getNumTrackEastRaw();
    int getNumTrackWestRaw();
    int getNumTrackEastReCenter();
    int getNumTrackWestReCenter();

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
    TProfile2D *p_mTpcQShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mTpcQShiftSinWest[mNumVzBin][mNumShiftCorr];

    // TPC EP Resolution
    TProfile *p_mTpcEpResolutionSub;
    double mTpcResSubVal[9];
    double mTpcResSubErr[9];

    TFile *file_mReCenterPar;
    TFile *file_mShiftPar;
    TFile *file_mResolution;

    const int mType;

  ClassDef(StTpcEpManager,1)
};

#endif
