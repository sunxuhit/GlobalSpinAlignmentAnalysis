#ifndef StEpdEpManager_h
#define StEpdEpManager_h

#include <string>
#include "TObject.h"
#include "TVector2.h"
#include "TVector3.h"

class StPicoEpdHit
class StEpdGeom;
class TFile;
class TProfile2D;
class TProfile;
class TH2F;

class StEpdEpManager : public TObject
{
  public:
    StEpdEpManager(int beamType);
    virtual ~StEpdEpManager();
    void clearEpdEpManager();
    void initEpdEpManager(int cent9, int runIndex, int vzBin);

    // Utilities
    TVector2 calq1Vector(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    double getTileWeight(StPicoEpdHit* picoEpdHit);
    double getPhiWeight(StPicoEpdHit* picoEpdHit);
    double getEtaWeight(StPicoEpdHit* picoEpdHit);

    // void addHitRawEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    // void addHitRawWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    // void addHitReCenterEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    // void addHitReCenterWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);

    // phi Weight Correction
    void initEpdPhiWeight();
    void fillEpdPhiWeightEast(StPicoEpdHit* picoEpdHit);
    void fillEpdPhiWeightWest(StPicoEpdHit* picoEpdHit);
    void writeEpdPhiWeight();
    void readEpdPhiWeight();

#if 0
    // ReCenter Correction
    void initEpdReCenter();
    void fillEpdReCenterEast(StPicoEpdHit* picoEpdHit);
    void fillEpdReCenterWest(StPicoEpdHit* picoEpdHit);
    void writeEpdReCenter();
    void readEpdReCenter();
    TVector2 getq2VecReCenterEast(); // 2nd ReCenter Parameter
    TVector2 getq2VecReCenterWest();
    TVector2 getq3VecReCenterEast(); // 3rd ReCenter Parameter
    TVector2 getq3VecReCenterWest();

    // Shift Correction
    void initEpdShift();
    void fillEpdShiftEast();
    void fillEpdShiftWest();
    void writeEpdShift();
    void readEpdShift();
    double getPsi2ShiftEast(); // 2nd shift Psi2
    double getPsi2ShiftWest();
    double transPsi2(double Psi2);
    double getPsi3ShiftEast(); // 3rd shift Psi3
    double getPsi3ShiftWest();
    double transPsi3(double Psi3);

    // Event Plane Resolution
    void initEpdResolution();
    void fillEpdResolution(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeEpdResolution();
    void readEpdResolution();
    double getEpdSubEp2ResVal(int cent9);
    double getEpdSubEp2ResErr(int cent9);
    double getEpdSubEp3ResVal(int cent9);
    double getEpdSubEp3ResErr(int cent9);

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
    void initEpdSubEpRaw(); // raw Sub EP
    void fillEpdSubEpRaw(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeEpdSubEpRaw();

    void initEpdSubEpReCenter(); // recenter Sub EP
    void fillEpdSubEpReCenter(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeEpdSubEpReCenter();

    void initEpdSubEpShift(); // shift Sub EP
    void fillEpdSubEpShift(double Psi2East, double Psi2West, double Psi3East, double Psi3West);
    void writeEpdSubEpShift();
#endif

  private:
    int mCent9;
    int mRunIndex;
    int mVzBin;

    int mQCouRawEast, mQCouRawWest;
    int mQCouReCenterEast, mQCouReCenterWest;

    TVector2 v_mQ2RawEast, v_mQ2RawWest; // 2nd EP
    TVector2 v_mQ2ReCenterEast, v_mQ2ReCenterWest; 

    TVector2 v_mQ3RawEast, v_mQ3RawWest; // 3rd EP
    TVector2 v_mQ3ReCenterEast, v_mQ3ReCenterWest; 

#if 0
    static const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
    static const int mNumShiftCorr = 20;
    static const int mNumCentrality = 9; // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    // ReCenter Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mEpdQ2ReCenterXEast[mNumVzBin]; // 2nd EP
    TProfile2D *p_mEpdQ2ReCenterYEast[mNumVzBin];
    TProfile2D *p_mEpdQ2ReCenterXWest[mNumVzBin];
    TProfile2D *p_mEpdQ2ReCenterYWest[mNumVzBin];

    TProfile2D *p_mEpdQ3ReCenterXEast[mNumVzBin]; // 3rd EP
    TProfile2D *p_mEpdQ3ReCenterYEast[mNumVzBin];
    TProfile2D *p_mEpdQ3ReCenterXWest[mNumVzBin];
    TProfile2D *p_mEpdQ3ReCenterYWest[mNumVzBin];

    // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mEpdQ2ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 2nd EP
    TProfile2D *p_mEpdQ2ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ2ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ2ShiftSinWest[mNumVzBin][mNumShiftCorr];

    TProfile2D *p_mEpdQ3ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 3rd EP
    TProfile2D *p_mEpdQ3ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ3ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ3ShiftSinWest[mNumVzBin][mNumShiftCorr];

    // TPC EP Resolution
    TProfile *p_mEpdSubEp2Res; // 2nd EP
    double mEpdSubEp2ResVal[mNumCentrality];
    double mEpdSubEp2ResErr[mNumCentrality];

    TProfile *p_mEpdSubEp3Res; // 3rd EP
    double mEpdSubEp3ResVal[mNumCentrality];
    double mEpdSubEp3ResErr[mNumCentrality];

    // Event Plane Distribution
    TH2F *h_mEpdEp2RawEast[mNumCentrality]; // 2nd raw EP
    TH2F *h_mEpdEp2RawWest[mNumCentrality];
    TH2F *h_mEpdEp2RawCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mEpdEp3RawEast[mNumCentrality]; // 3rd raw EP
    TH2F *h_mEpdEp3RawWest[mNumCentrality];
    TH2F *h_mEpdEp3RawCorr[mNumCentrality]; // Psi3East vs Psi3West

    TH2F *h_mEpdEp2ReCenterEast[mNumCentrality]; // 2nd recenter EP
    TH2F *h_mEpdEp2ReCenterWest[mNumCentrality];
    TH2F *h_mEpdEp2ReCenterCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mEpdEp3ReCenterEast[mNumCentrality]; // 3rd recenter EP
    TH2F *h_mEpdEp3ReCenterWest[mNumCentrality];
    TH2F *h_mEpdEp3ReCenterCorr[mNumCentrality]; // Psi3East vs Psi3West

    TH2F *h_mEpdEp2ShiftEast[mNumCentrality]; // 3rd shift EP
    TH2F *h_mEpdEp2ShiftWest[mNumCentrality];
    TH2F *h_mEpdEp2ShiftCorr[mNumCentrality]; // Psi2East vs Psi2West

    TH2F *h_mEpdEp3ShiftEast[mNumCentrality]; // 3rd shift EP
    TH2F *h_mEpdEp3ShiftWest[mNumCentrality];
    TH2F *h_mEpdEp3ShiftCorr[mNumCentrality]; // Psi3East vs Psi3West

    TFile *file_mReCenterPar;
    TFile *file_mShiftPar;
    TFile *file_mResolution;
#endif

    const int mType;
    StEpdGeom* mEpdGeom;

  ClassDef(StEpdEpManager,1)
};

#endif
