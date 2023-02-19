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

    void addHitRawEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt
    void addHitRawWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void addHitWgtEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt * phiWgt * etaWgt
    void addHitWgtWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void addHitReCtrEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt * phiWgt * etaWgt
    void addHitReCtrWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);

    // phi Weight Correction
    void initEpdPhiWgt();
    void fillEpdPhiWgtEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void fillEpdPhiWgtWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void writeEpdPhiWgt();
    void readEpdPhiWgt();

    // ReCenter Correction
    void initEpdReCtr();
    void fillEpdReCtrEast();
    void fillEpdReCtrWest();
    void writeEpdReCtr();
    void readEpdReCtr();
    TVector2 getq1VecCtrEast(); // 1st ReCtr Parameter
    TVector2 getq1VecCtrWest();

    // Shift Correction
    void initEpdShift();
    void fillEpdShiftEast();
    void fillEpdShiftWest();
    void writeEpdShift();
    void readEpdShift();
    double getPsi1ShiftEast(); // 1st shift Psi1
    double getPsi1ShiftWest();
    double getPsi1ShiftFull();
    double transPsi1(double Psi1);

    void initEpdShiftFull(); // Full
    void fillEpdShiftFull();
    void writeEpdShiftFull();
    void readEpdShiftFull();
    TVector2 getQ1VecShiftFullCorr();
    double getPsi1ShiftFullCorr();

    // Event Plane Resolution
    void initEpdResolution();
    void fillEpdResolution(double Psi1East, double Psi1West);
    void writeEpdResolution();
    void readEpdResolution();
    double getEpdSubEp1ResVal(int cent9);
    double getEpdSubEp1ResErr(int cent9);
    double getEpdFullEp1ResVal(int cent9);
    double getEpdFullEp1ResErr(int cent9);

    // Q1Vector
    TVector2 getQ1VecRawEast(); // Q1Vector
    TVector2 getQ1VecRawWest();
    TVector2 getQ1VecRawFull();
    TVector2 getQ1VecWgtEast();
    TVector2 getQ1VecWgtWest();
    TVector2 getQ1VecWgtFull();
    TVector2 getQ1VecReCtrEast();
    TVector2 getQ1VecReCtrWest();
    TVector2 getQ1VecReCtrFull();
    TVector2 getQ1VecShiftEast();
    TVector2 getQ1VecShiftWest();
    TVector2 getQ1VecShiftFull();
    double getPsi1RawEast();
    double getPsi1RawWest();
    double getPsi1RawFull();
    double getPsi1WgtEast();
    double getPsi1WgtWest();
    double getPsi1WgtFull();
    double getPsi1ReCtrEast();
    double getPsi1ReCtrWest();
    double getPsi1ReCtrFull();

    // Event Plane Distribution
    void initEpdSubEpRaw(); // raw Sub EP
    void fillEpdSubEpRaw(double Psi1East, double Psi1West, double Psi1Full);
    void writeEpdSubEpRaw();

    void initEpdSubEpWgt(); // phi weighted Sub EP
    void fillEpdSubEpWgt(double Psi1East, double Psi1West, double Psi1Full);
    void writeEpdSubEpWgt();

    void initEpdSubEpReCtr(); // recenter Sub EP
    void fillEpdSubEpReCtr(double Psi1East, double Psi1West, double Psi1Full);
    void writeEpdSubEpReCtr();

    void initEpdSubEpShift(); // shift Sub EP
    void fillEpdSubEpShift(double Psi1East, double Psi1West, double Psi1Full);
    void writeEpdSubEpShift();

    void initEpdFullEpShift(); // shift Sub EP
    void fillEpdFullEpShift(double Psi1FullCorr);
    void writeEpdFullEpShift();

  private:
    static const int mNumVzBin      = 2;  // 0: vz < 0 | 1: vz >= 0
    static const int mNumShiftCorr  = 20;
    static const int mNumCentrality = 9;  // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    static const int mNumSectors    = 12; // picoEpdHit->position(): return 1-12 | use sec-1
    static const int mNumTiles      = 31; // picoEpdHit->tile(): return 1-31 | use tile-1
    static const int mNumRings      = 16; // picoEpdHit->row(): return 1-16 | use ring-1 = > 0: most inner ring (12 tiles), 1-15: 24 tiles

    double mQ1WgtSideRawEast, mQ1WgtSideRawWest; // tileWgt only
    double mQ1WgtSideWgtEast, mQ1WgtSideWgtWest; // tileWgt * phiWgt(if avaliable) * etaWgt(if avaliable)
    double mQ1WgtSideReCtrEast, mQ1WgtSideReCtrWest; // tileWgt * phiWgt(if avaliable) * etaWgt(if avaliable)
    double mQ1WgtRingRawEast[mNumRings], mQ1WgtRingRawWest[mNumRings]; // tileWgt only for each ring
    double mQ1WgtRingWgtEast[mNumRings], mQ1WgtRingWgtWest[mNumRings]; // tileWgt * phiWgt(if avaliable) * etaWgt(if avaliable) for each ring
    // double mQ1WgtRingReCtrEast[mNumRings], mQ1WgtRingReCtrWest[mNumRings]; // tileWgt * phiWgt(if avaliable) * etaWgt(if avaliable) for each ring
    
    TVector2 v_mQ1SideRawEast, v_mQ1SideRawWest; // raw Q1 Vector
    TVector2 v_mQ1SideWgtEast, v_mQ1SideWgtWest; // phi&eta weighted Q1 Vector
    TVector2 v_mQ1SideReCtrEast, v_mQ1SideReCtrWest; // phi&eta weighted Q1 Vector
    TVector2 v_mQ1RingRawEast[mNumRings], v_mQ1RingRawWest[mNumRings]; // raw Q1 Vector for each ring
    TVector2 v_mQ1RingWgtEast[mNumRings], v_mQ1RingWgtWest[mNumRings]; // phi&eta weighted Q1 Vector for each ring
    // TVector2 v_mQ1RingReCtrEast[mNumRings], v_mQ1RingReCtrWest[mNumRings]; // phi&eta weighted Q1 Vector for each ring


    // phi Weight Correction | x-axis is super sector Id, y-axis is tile Id
    TH2F *h_mEpdPhiWgtEast[mNumCentrality];
    TH2F *h_mEpdPhiAveEast[mNumCentrality];
    TH2F *h_mEpdPhiWgtWest[mNumCentrality];
    TH2F *h_mEpdPhiAveWest[mNumCentrality];

    // ReCtr Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mEpdQ1ReCtrXEast[mNumVzBin]; // 1st EP
    TProfile2D *p_mEpdQ1ReCtrYEast[mNumVzBin];
    TProfile2D *p_mEpdQ1ReCtrXWest[mNumVzBin];
    TProfile2D *p_mEpdQ1ReCtrYWest[mNumVzBin];

    // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mEpdQ1ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 1st EP
    TProfile2D *p_mEpdQ1ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1ShiftSinWest[mNumVzBin][mNumShiftCorr];

    // Shift Correction for Full EP
    TProfile2D *p_mEpdQ1ShiftCosFull[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1ShiftSinFull[mNumVzBin][mNumShiftCorr];

    // EPD EP Resolution
    TProfile *p_mEpdSubEp1Res; // 1st EP
    double mEpdSubEp1ResVal[mNumCentrality];
    double mEpdSubEp1ResErr[mNumCentrality];
    double mEpdFullEp1ResVal[mNumCentrality];
    double mEpdFullEp1ResErr[mNumCentrality];

    // Event Plane Distribution
    TH2F *h_mEpdEp1RawEast[mNumCentrality]; // 1st raw EP
    TH2F *h_mEpdEp1RawWest[mNumCentrality];
    TH2F *h_mEpdEp1RawFull[mNumCentrality];
    TH2F *h_mEpdEp1RawCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mEpdEp1WgtEast[mNumCentrality]; // 1st weighted EP
    TH2F *h_mEpdEp1WgtWest[mNumCentrality];
    TH2F *h_mEpdEp1WgtFull[mNumCentrality];
    TH2F *h_mEpdEp1WgtCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mEpdEp1ReCtrEast[mNumCentrality]; // 1st recenter EP
    TH2F *h_mEpdEp1ReCtrWest[mNumCentrality];
    TH2F *h_mEpdEp1ReCtrFull[mNumCentrality];
    TH2F *h_mEpdEp1ReCtrCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mEpdEp1ShiftEast[mNumCentrality]; // 1st shift EP
    TH2F *h_mEpdEp1ShiftWest[mNumCentrality];
    TH2F *h_mEpdEp1ShiftFull[mNumCentrality];
    TH2F *h_mEpdEp1ShiftCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mEpdEp1ShiftFullCorr[mNumCentrality];

    TFile *file_mPhiWgtPar;
    TFile *file_mReCtrPar;
    TFile *file_mShiftPar;
    // TFile *file_mResolution;

    int mCent9;
    int mRunIndex;
    int mVzBin;
    bool mUsePhiWgt;
    bool mUseEtaWgt;
    const int mType;
    StEpdGeom* mEpdGeom;

  ClassDef(StEpdEpManager,1)
};

#endif
