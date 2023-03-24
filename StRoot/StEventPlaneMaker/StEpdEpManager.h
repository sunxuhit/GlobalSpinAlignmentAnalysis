#ifndef StEpdEpManager_h
#define StEpdEpManager_h

#include <string>
#include "TObject.h"
#include "TVector2.h"
#include "TVector3.h"

class TFile;
class TProfile2D;
class TProfile;
class TH2F;

class StPicoEpdHit;
class StEpdGeom;

class StEpdEpManager : public TObject
{
  public:
    StEpdEpManager(int beamType);
    virtual ~StEpdEpManager();
    void clearEpdEpManager();
    void initEpdEpManager(int cent9, int runIndex, int vzBin);

    // Utilities
    TVector3 getEpdRanVec(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    TVector2 calq1Vector(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    double getTileWeight(StPicoEpdHit* picoEpdHit);
    double getPhiWeight(StPicoEpdHit* picoEpdHit);
    double getEtaWeight(StPicoEpdHit* picoEpdHit);
    int getEpdEpGrp(StPicoEpdHit* picoEpdHit);

    void addHitRawEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt
    void addHitRawWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void addHitWgtEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt * phiWgt * etaWgt
    void addHitWgtWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void addHitReCtrEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt * phiWgt * etaWgt
    void addHitReCtrWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);

    // phi Weight Correction
    void initEpdPhiWgt();
    void fillEpdPhiWgtEast(StPicoEpdHit* picoEpdHit);
    void fillEpdPhiWgtWest(StPicoEpdHit* picoEpdHit);
    void writeEpdPhiWgt();
    void readEpdPhiWgt();

    // ReCenter Correction
    void initEpdReCtr();
    void fillEpdSideReCtrEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx);
    void fillEpdSideReCtrWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx);
    void fillEpdGrpReCtrEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx);
    void fillEpdGrpReCtrWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx);
    void writeEpdReCtr();
    void readEpdReCtr();
    TVector2 getq1VecSideCtrEast(); // 1st ReCtr Parameter
    TVector2 getq1VecSideCtrWest();
    TVector2 getq1VecGrpCtrEast(int grpId);
    TVector2 getq1VecGrpCtrWest(int grpId);

    // Shift Correction
    void initEpdShift();
    void fillEpdSideShiftEast();
    void fillEpdSideShiftWest();
    void fillEpdGrpShiftEast(int grpId);
    void fillEpdGrpShiftWest(int grpId);
    void writeEpdShift();
    void readEpdShift();
    double getPsi1SideShiftEast(); // 1st shift Psi1
    double getPsi1SideShiftWest();
    double getPsi1SideShiftFull();
    double getPsi1GrpShiftEast(int grpId);
    double getPsi1GrpShiftWest(int grpId);
    double getPsi1GrpShiftFull(int grpId);
    double transPsi1(double Psi1);
    bool isPsi1InRange(double Psi1);
    TVector2 getQ1VecSideShiftEast();
    TVector2 getQ1VecSideShiftWest();
    TVector2 getQ1VecSideShiftFull();
    TVector2 getQ1VecGrpShiftEast(int grpId);
    TVector2 getQ1VecGrpShiftWest(int grpId);
    TVector2 getQ1VecGrpShiftFull(int grpId);

    void initEpdShiftFull(); // Full
    void fillEpdSideShiftFull();
    void fillEpdGrpShiftFull(int grpId);
    void writeEpdShiftFull();
    void readEpdShiftFull();
    double getPsi1SideShiftFullCorr();
    double getPsi1GrpShiftFullCorr(int grpId);
    TVector2 getQ1VecSideShiftFullCorr();
    TVector2 getQ1VecGrpShiftFullCorr(int grpId);

    // Event Plane Resolution
    void initEpdResolution();
    void fillEpdSideResolution(double Psi1East, double Psi1West);
    void fillEpdGrpResolution(double Psi1East, double Psi1West, int grpId);
    void writeEpdResolution();
    void readEpdResolution();
    double getEpdSubEp1SideResVal(int cent9);
    double getEpdSubEp1SideResErr(int cent9);
    double getEpdFullEp1SideResVal(int cent9);
    double getEpdFullEp1SideResErr(int cent9);
    double getEpdSubEp1GrpResVal(int cent9, int grpId);
    double getEpdSubEp1GrpResErr(int cent9, int grpId);
    double getEpdFullEp1GrpResVal(int cent9, int grpId);
    double getEpdFullEp1GrpResErr(int cent9, int grpId);

    // Charged Hadron Directed Flow
    void initEpdSubEpFlow(); // Sub EP
    void fillEpdSubEpV1(double eta, double v1, double reweight);
    void writeEpdSubEpFlow();

    // Q1Vector
    TVector2 getQ1VecSideRawEast(); // Q1Vector for full EPD
    TVector2 getQ1VecSideRawWest();
    TVector2 getQ1VecSideRawFull();
    TVector2 getQ1VecSideWgtEast();
    TVector2 getQ1VecSideWgtWest();
    TVector2 getQ1VecSideWgtFull();
    TVector2 getQ1VecSideReCtrEast();
    TVector2 getQ1VecSideReCtrWest();
    TVector2 getQ1VecSideReCtrFull();
    double getPsi1SideRawEast();
    double getPsi1SideRawWest();
    double getPsi1SideRawFull();
    double getPsi1SideWgtEast();
    double getPsi1SideWgtWest();
    double getPsi1SideWgtFull();
    double getPsi1SideReCtrEast();
    double getPsi1SideReCtrWest();
    double getPsi1SideReCtrFull();

    TVector2 getQ1VecGrpRawEast(int grpId); // Q1Vector for each group
    TVector2 getQ1VecGrpRawWest(int grpId);
    TVector2 getQ1VecGrpRawFull(int grpId);
    TVector2 getQ1VecGrpWgtEast(int grpId);
    TVector2 getQ1VecGrpWgtWest(int grpId);
    TVector2 getQ1VecGrpWgtFull(int grpId);
    TVector2 getQ1VecGrpReCtrEast(int grpId);
    TVector2 getQ1VecGrpReCtrWest(int grpId);
    TVector2 getQ1VecGrpReCtrFull(int grpId);
    double getPsi1GrpRawEast(int grpId);
    double getPsi1GrpRawWest(int grpId);
    double getPsi1GrpRawFull(int grpId);
    double getPsi1GrpWgtEast(int grpId);
    double getPsi1GrpWgtWest(int grpId);
    double getPsi1GrpWgtFull(int grpId);
    double getPsi1GrpReCtrEast(int grpId);
    double getPsi1GrpReCtrWest(int grpId);
    double getPsi1GrpReCtrFull(int grpId);

    // Event Plane Distribution
    void initEpdSubEpRaw(); // raw Sub EP
    void fillEpdSubEpSideRaw(double Psi1East, double Psi1West, double Psi1Full);
    void fillEpdSubEpGrpRaw(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void writeEpdSubEpRaw();

    void initEpdSubEpWgt(); // phi weighted Sub EP
    void fillEpdSubEpSideWgt(double Psi1East, double Psi1West, double Psi1Full);
    void fillEpdSubEpGrpWgt(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void writeEpdSubEpWgt();

    void initEpdSubEpReCtr(); // recenter Sub EP
    void fillEpdSubEpSideReCtr(double Psi1East, double Psi1West, double Psi1Full);
    void fillEpdSubEpGrpReCtr(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void writeEpdSubEpReCtr();

    void initEpdSubEpShift(); // shift Sub EP
    void fillEpdSubEpSideShift(double Psi1East, double Psi1West, double Psi1Full);
    void fillEpdSubEpGrpShift(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void writeEpdSubEpShift();

    void initEpdFullEpShift(); // shift Sub EP
    void fillEpdFullEpSideShift(double Psi1FullCorr);
    void fillEpdFullEpGrpShift(double Psi1FullCorr, int grpId);
    void writeEpdFullEpShift();

  private:
    static const int mNumVzBin      = 2;  // 0: vz < 0 | 1: vz >= 0
    static const int mNumShiftCorr  = 20;
    static const int mNumCentrality = 9;  // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    static const int mNumSectors    = 12; // picoEpdHit->position(): return 1-12 | use sec-1
    static const int mNumTiles      = 31; // picoEpdHit->tile(): return 1-31 | use tile-1
    static const int mNumRings      = 16; // picoEpdHit->row(): return 1-16 | use ring-1 = > 0: most inner ring (12 tiles), 1-15: 24 tiles
    static const int mNumRingsUsed  = 4;  // rings used in Q1Vector calculation: [0, mNumRingsUsed) | set to 4 to get eta weight then to 16 for EPD EP
    static const int mNumRingsGrps  = 2;  // Group 0: 0-7 rings | Group 1: 8-15 rings

    // Q1Vector of all EPD rings
    double mQ1WgtSideRawEast, mQ1WgtSideRawWest; // tileWgt only
    double mQ1WgtSideWgtEast, mQ1WgtSideWgtWest; // tileWgt * phiWgt(if avaliable) * etaWgt(if avaliable)
    double mQ1WgtSideReCtrEast, mQ1WgtSideReCtrWest; // tileWgt * phiWgt(if avaliable) * etaWgt(if avaliable)
    TVector2 v_mQ1SideRawEast, v_mQ1SideRawWest; // raw Q1 Vector
    TVector2 v_mQ1SideWgtEast, v_mQ1SideWgtWest; // phi&eta weighted Q1 Vector
    TVector2 v_mQ1SideReCtrEast, v_mQ1SideReCtrWest; // phi&eta weighted Q1 Vector

    // Q1Vector of groups of EPD rings
    double mQ1WgtGrpRawEast[mNumRingsGrps], mQ1WgtGrpRawWest[mNumRingsGrps]; // tileWgt only for each ring
    double mQ1WgtGrpWgtEast[mNumRingsGrps], mQ1WgtGrpWgtWest[mNumRingsGrps]; // tileWgt * phiWgt(if avaliable) * etaWgt(if avaliable) for each ring
    double mQ1WgtGrpReCtrEast[mNumRingsGrps], mQ1WgtGrpReCtrWest[mNumRingsGrps]; // tileWgt * phiWgt(if avaliable) * etaWgt(if avaliable) for each ring
    TVector2 v_mQ1GrpRawEast[mNumRingsGrps], v_mQ1GrpRawWest[mNumRingsGrps]; // raw Q1 Vector for each ring
    TVector2 v_mQ1GrpWgtEast[mNumRingsGrps], v_mQ1GrpWgtWest[mNumRingsGrps]; // phi&eta weighted Q1 Vector for each ring
    TVector2 v_mQ1GrpReCtrEast[mNumRingsGrps], v_mQ1GrpReCtrWest[mNumRingsGrps]; // phi&eta weighted Q1 Vector for each ring

    // phi Weight Correction | x-axis is super sector Id, y-axis is tile Id
    TH2F *h_mEpdPhiWgtEast[mNumCentrality];
    TH2F *h_mEpdPhiAveEast[mNumCentrality];
    TH2F *h_mEpdPhiWgtWest[mNumCentrality];
    TH2F *h_mEpdPhiAveWest[mNumCentrality];

    // ReCtr Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mEpdQ1SideReCtrXEast[mNumVzBin]; // 1st EP
    TProfile2D *p_mEpdQ1SideReCtrYEast[mNumVzBin];
    TProfile2D *p_mEpdQ1SideReCtrXWest[mNumVzBin];
    TProfile2D *p_mEpdQ1SideReCtrYWest[mNumVzBin];
    TProfile2D *p_mEpdQ1GrpReCtrXEast[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrYEast[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrXWest[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrYWest[mNumVzBin][mNumRingsGrps];

    // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mEpdQ1SideShiftCosEast[mNumVzBin][mNumShiftCorr]; // 1st EP
    TProfile2D *p_mEpdQ1SideShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1SideShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1SideShiftSinWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1GrpShiftCosEast[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftSinEast[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftCosWest[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftSinWest[mNumVzBin][mNumShiftCorr][mNumRingsGrps];

    // Shift Correction for Full EP
    TProfile2D *p_mEpdQ1SideShiftCosFull[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1SideShiftSinFull[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1GrpShiftCosFull[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftSinFull[mNumVzBin][mNumShiftCorr][mNumRingsGrps];

    // EPD EP Resolution
    TProfile *p_mEpdSubEp1SideRes; // 1st EP
    double mEpdSubEp1SideResVal[mNumCentrality];
    double mEpdSubEp1SideResErr[mNumCentrality];
    double mEpdFullEp1SideResVal[mNumCentrality];
    double mEpdFullEp1SideResErr[mNumCentrality];
    TProfile *p_mEpdSubEp1GrpRes[mNumRingsGrps]; // resolution of same group
    double mEpdSubEp1GrpResVal[mNumCentrality][mNumRingsGrps];
    double mEpdSubEp1GrpResErr[mNumCentrality][mNumRingsGrps];
    double mEpdFullEp1GrpResVal[mNumCentrality][mNumRingsGrps];
    double mEpdFullEp1GrpResErr[mNumCentrality][mNumRingsGrps];

    // Charged Hadron Directed Flow
    TProfile *p_mEpdSubEpV1[mNumCentrality]; // v1 vs. eta

    // Event Plane Distribution
    TH2F *h_mEpdEp1SideRawEast[mNumCentrality]; // 1st raw EP
    TH2F *h_mEpdEp1SideRawWest[mNumCentrality];
    TH2F *h_mEpdEp1SideRawFull[mNumCentrality];
    TH2F *h_mEpdEp1SideRawCorr[mNumCentrality]; // Psi1East vs Psi1West
    TH2F *h_mEpdEp1GrpRawEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpRawWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpRawFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpRawCorr[mNumCentrality][mNumRingsGrps];
    // TH2F *h_mEpdEp1GrpXRawCorr[mNumCentrality][mNumRingsGrps]; // 0: Psi1Grp0East vs. Psi1Grp1East | 1: Psi1Grp0West vs. Psi1Grp1West

    TH2F *h_mEpdEp1SideWgtEast[mNumCentrality]; // 1st weighted EP
    TH2F *h_mEpdEp1SideWgtWest[mNumCentrality];
    TH2F *h_mEpdEp1SideWgtFull[mNumCentrality];
    TH2F *h_mEpdEp1SideWgtCorr[mNumCentrality]; // Psi1East vs Psi1West
    TH2F *h_mEpdEp1GrpWgtEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpWgtWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpWgtFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpWgtCorr[mNumCentrality][mNumRingsGrps];
    // TH2F *h_mEpdEp1GrpXWgtCorr[mNumCentrality][mNumRingsGrps]; // 0: Psi1Grp0East vs. Psi1Grp1East | 1: Psi1Grp0West vs. Psi1Grp1West

    TH2F *h_mEpdEp1SideReCtrEast[mNumCentrality]; // 1st recenter EP
    TH2F *h_mEpdEp1SideReCtrWest[mNumCentrality];
    TH2F *h_mEpdEp1SideReCtrFull[mNumCentrality];
    TH2F *h_mEpdEp1SideReCtrCorr[mNumCentrality]; // Psi1East vs Psi1West
    TH2F *h_mEpdEp1GrpReCtrEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrCorr[mNumCentrality][mNumRingsGrps];
    // TH2F *h_mEpdEp1GrpXReCtrCorr[mNumCentrality][mNumRingsGrps]; // 0: Psi1Grp0East vs. Psi1Grp1East | 1: Psi1Grp0West vs. Psi1Grp1West

    TH2F *h_mEpdEp1SideShiftEast[mNumCentrality]; // 1st shift EP
    TH2F *h_mEpdEp1SideShiftWest[mNumCentrality];
    TH2F *h_mEpdEp1SideShiftFull[mNumCentrality];
    TH2F *h_mEpdEp1SideShiftCorr[mNumCentrality]; // Psi1East vs Psi1West
    TH2F *h_mEpdEp1GrpShiftEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftCorr[mNumCentrality][mNumRingsGrps];
    // TH2F *h_mEpdEp1GrpXShiftCorr[mNumCentrality][mNumRingsGrps]; // 0: Psi1Grp0East vs. Psi1Grp1East | 1: Psi1Grp0West vs. Psi1Grp1West

    TH2F *h_mEpdEp1SideShiftFullCorr[mNumCentrality]; // 1st shift full EP
    TH2F *h_mEpdEp1GrpShiftFullCorr[mNumCentrality][mNumRingsGrps];

    TFile *file_mPhiWgtPar;
    TFile *file_mReCtrPar;
    TFile *file_mShiftPar;
    TFile *file_mResolution;

    int mCent9;
    int mRunIndex;
    int mVzBin;
    bool mUsePhiWgt;
    bool mUseEtaWgt;
    const int mType;
    StEpdGeom *mEpdGeom;

  ClassDef(StEpdEpManager,1)
};

#endif
