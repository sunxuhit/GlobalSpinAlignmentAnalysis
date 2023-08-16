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
    TVector3 getEpdCtrVec(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    TVector3 getEpdRanVec(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    TVector2 calq1Vector(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    double getTileWgt(StPicoEpdHit* picoEpdHit);
    double getPhiWgt(StPicoEpdHit* picoEpdHit);
    double getEtaWgt(StPicoEpdHit* picoEpdHit);
    double getEpdWgt(StPicoEpdHit* picoEpdHit);
    int getEpdEpGrp(StPicoEpdHit* picoEpdHit); // used in FXT
    double transPsi1(double Psi1);
    bool isPsi1InRange(double Psi1);

    // Calculate Q1Vector
    // QVector of each Side & used in IsoBar
    void addHitSideRawEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt
    void addHitSideRawWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void addHitSideWgtEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt * phiWgt * etaWgt
    void addHitSideWgtWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void addHitSideReCtrEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt * phiWgt * etaWgt
    void addHitSideReCtrWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    // QVector of each Group on each Side & used in FXT
    void addHitGrpRawEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt
    void addHitGrpRawWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void addHitGrpWgtEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt * phiWgt * etaWgt
    void addHitGrpWgtWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);
    void addHitGrpReCtrTrkAveEast(StPicoEpdHit* picoEpdHit, TVector3 primVtx); // tileWgt * phiWgt * etaWgt
    void addHitGrpReCtrTrkAveWest(StPicoEpdHit* picoEpdHit, TVector3 primVtx);

    // phi Weight Correction
    void initEpdPhiWgt();
    void fillEpdPhiWgtEast(StPicoEpdHit* picoEpdHit);
    void fillEpdPhiWgtWest(StPicoEpdHit* picoEpdHit);
    void writeEpdPhiWgt();
    void readEpdPhiWgt();

    // ReCenter Correction
    // each Side & used in IsoBar
    void initEpdSideReCtr();
    void fillEpdSideReCtrEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx);
    void fillEpdSideReCtrWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx);
    void writeEpdSideReCtr();
    void readEpdSideReCtr();
    TVector2 getq1VecSideCtrEast();
    TVector2 getq1VecSideCtrWest();
    // each Group on each Side & used in FXT
    void initEpdGrpReCtr();
    void fillEpdGrpReCtrTrkAveEast(StPicoEpdHit *picoEpdHit, TVector3 primVtx);
    void fillEpdGrpReCtrTrkAveWest(StPicoEpdHit *picoEpdHit, TVector3 primVtx);
    void fillEpdGrpReCtrEvtAveEast(TVector2 Q1VecGrp, int grpId); // need addHitGrpRawEast in StEventPlaneMaker
    void fillEpdGrpReCtrEvtAveWest(TVector2 Q1VecGrp, int grpId); // need addHitGrpRawWest in StEventPlaneMaker
    void writeEpdGrpReCtr();
    void readEpdGrpReCtr();
    TVector2 getq1VecGrpCtrTrkAveEast(int grpId); // FXT
    TVector2 getq1VecGrpCtrTrkAveWest(int grpId);
    TVector2 getQ1VecGrpCtrEvtAveEast(int grpId); // FXT
    TVector2 getQ1VecGrpCtrEvtAveWest(int grpId);

    // Shift Correction for sub EP
    // each Side & used in IsoBar
    void initEpdSideShift();
    void fillEpdSideShiftEast();
    void fillEpdSideShiftWest();
    void writeEpdSideShift();
    void readEpdSideShift();
    double getPsi1SideShiftEast(); // 1st shift Psi1
    double getPsi1SideShiftWest();
    double getPsi1SideShiftFull();
    TVector2 getQ1VecSideShiftEast();
    TVector2 getQ1VecSideShiftWest();
    TVector2 getQ1VecSideShiftFull();
    // each Group on each Side & used in FXT
    void initEpdGrpShift();
    void fillEpdGrpShiftTrkAveEast(int grpId);
    void fillEpdGrpShiftTrkAveWest(int grpId);
    void fillEpdGrpShiftEvtAveEast(int grpId);
    void fillEpdGrpShiftEvtAveWest(int grpId);
    void writeEpdGrpShift();
    void readEpdGrpShift();
    double getPsi1GrpShiftTrkAveEast(int grpId);
    double getPsi1GrpShiftTrkAveWest(int grpId);
    double getPsi1GrpShiftTrkAveFull(int grpId);
    double getPsi1GrpShiftEvtAveEast(int grpId);
    double getPsi1GrpShiftEvtAveWest(int grpId);
    double getPsi1GrpShiftEvtAveFull(int grpId);
    TVector2 getQ1VecGrpShiftTrkAveEast(int grpId);
    TVector2 getQ1VecGrpShiftTrkAveWest(int grpId);
    TVector2 getQ1VecGrpShiftTrkAveFull(int grpId);
    TVector2 getQ1VecGrpShiftEvtAveEast(int grpId);
    TVector2 getQ1VecGrpShiftEvtAveWest(int grpId);
    TVector2 getQ1VecGrpShiftEvtAveFull(int grpId);

    // Shift Correction for full EP
    // each Side & used in IsoBar
    void initEpdSideShiftFull();
    void fillEpdSideShiftFull();
    void writeEpdSideShiftFull();
    void readEpdSideShiftFull();
    double getPsi1SideShiftFullCorr();
    TVector2 getQ1VecSideShiftFullCorr();

    // Event Plane Resolution
    // each Side & used in IsoBar
    void initEpdSideResolution();
    void fillEpdSideResolution(double Psi1East, double Psi1West);
    void writeEpdSideResolution();
    void readEpdSideResolution();
    double getEpdSubEp1SideResVal(int cent9);
    double getEpdSubEp1SideResErr(int cent9);
    double getEpdFullEp1SideResVal(int cent9);
    double getEpdFullEp1SideResErr(int cent9);
    // each Group on each Side & used in FXT
    void initEpdGrpResolution();
    void fillEpdGrpResolution(double Psi1East, double Psi1West, int grpId);
    void writeEpdGrpResolution();
    void readEpdGrpResolution();
    double getEpdSubEp1GrpResVal(int cent9, int grpId);
    double getEpdSubEp1GrpResErr(int cent9, int grpId);
    double getEpdFullEp1GrpResVal(int cent9, int grpId);
    double getEpdFullEp1GrpResErr(int cent9, int grpId);

    // Charged Hadron Directed Flow w.r.t. EPD Side
    void initEpdSubEpSideFlow(); // Sub EP
    void fillEpdSubEpSideV1(double eta, double v1, double reweight);
    void writeEpdSubEpSideFlow();

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
    TVector2 getQ1VecGrpReCtrTrkAveEast(int grpId);
    TVector2 getQ1VecGrpReCtrTrkAveWest(int grpId);
    TVector2 getQ1VecGrpReCtrTrkAveFull(int grpId);
    TVector2 getQ1VecGrpReCtrEvtAveEast(int grpId);
    TVector2 getQ1VecGrpReCtrEvtAveWest(int grpId);
    TVector2 getQ1VecGrpReCtrEvtAveFull(int grpId);
    double getPsi1GrpRawEast(int grpId);
    double getPsi1GrpRawWest(int grpId);
    double getPsi1GrpRawFull(int grpId);
    double getPsi1GrpWgtEast(int grpId);
    double getPsi1GrpWgtWest(int grpId);
    double getPsi1GrpWgtFull(int grpId);
    double getPsi1GrpReCtrTrkAveEast(int grpId);
    double getPsi1GrpReCtrTrkAveWest(int grpId);
    double getPsi1GrpReCtrTrkAveFull(int grpId);
    double getPsi1GrpReCtrEvtAveEast(int grpId);
    double getPsi1GrpReCtrEvtAveWest(int grpId);
    double getPsi1GrpReCtrEvtAveFull(int grpId);

    // Event Plane Distribution
    void initEpdSubEpSideRaw(); // raw Sub EP
    void fillEpdSubEpSideRaw(double Psi1East, double Psi1West, double Psi1Full);
    void writeEpdSubEpSideRaw();
    void initEpdSubEpGrpRaw();
    void fillEpdSubEpGrpRaw(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void writeEpdSubEpGrpRaw();

    void initEpdSubEpSideWgt(); // phi weighted Sub EP
    void fillEpdSubEpSideWgt(double Psi1East, double Psi1West, double Psi1Full);
    void writeEpdSubEpSideWgt();
    void initEpdSubEpGrpWgt();
    void fillEpdSubEpGrpWgt(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void writeEpdSubEpGrpWgt();

    void initEpdSubEpSideReCtr(); // recenter Sub EP
    void fillEpdSubEpSideReCtr(double Psi1East, double Psi1West, double Psi1Full);
    void writeEpdSubEpSideReCtr();
    void initEpdSubEpGrpReCtr();
    void fillEpdSubEpGrpReCtrTrkAve(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void fillEpdSubEpGrpReCtrEvtAve(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void writeEpdSubEpGrpReCtr();

    void initEpdSubEpSideShift(); // shift Sub EP
    void fillEpdSubEpSideShift(double Psi1East, double Psi1West, double Psi1Full);
    void writeEpdSubEpSideShift();
    void initEpdSubEpGrpShift();
    void fillEpdSubEpGrpShiftTrkAve(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void fillEpdSubEpGrpShiftEvtAve(double Psi1East, double Psi1West, double Psi1Full, int grpId);
    void writeEpdSubEpGrpShift();

    void initEpdFullEpSideShift(); // shift Full EP
    void fillEpdFullEpSideShift(double Psi1FullCorr);
    void writeEpdFullEpSideShift();

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
    double mQ1WgtGrpReCtrTrkAveEast[mNumRingsGrps], mQ1WgtGrpReCtrTrkAveWest[mNumRingsGrps]; // tileWgt * phiWgt(if avaliable) * etaWgt(if avaliable) for each ring
    TVector2 v_mQ1GrpRawEast[mNumRingsGrps], v_mQ1GrpRawWest[mNumRingsGrps]; // raw Q1 Vector for each ring
    TVector2 v_mQ1GrpWgtEast[mNumRingsGrps], v_mQ1GrpWgtWest[mNumRingsGrps]; // phi&eta weighted Q1 Vector for each ring
    TVector2 v_mQ1GrpReCtrTrkAveEast[mNumRingsGrps], v_mQ1GrpReCtrTrkAveWest[mNumRingsGrps]; // phi&eta weighted Q1 Vector for each ring
    // TVector2 v_mQ1GrpReCtrEvtAveEast[mNumRingsGrps], v_mQ1GrpReCtrEvtAveWest[mNumRingsGrps]; // phi&eta weighted Q1 Vector for each ring

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
    TProfile2D *p_mEpdQ1GrpReCtrTrkAveXEast[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrTrkAveYEast[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrTrkAveXWest[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrTrkAveYWest[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrEvtAveXEast[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrEvtAveYEast[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrEvtAveXWest[mNumVzBin][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpReCtrEvtAveYWest[mNumVzBin][mNumRingsGrps];

    // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mEpdQ1SideShiftCosEast[mNumVzBin][mNumShiftCorr]; // 1st EP
    TProfile2D *p_mEpdQ1SideShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1SideShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1SideShiftSinWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1GrpShiftCosTrkAveEast[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftSinTrkAveEast[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftCosTrkAveWest[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftSinTrkAveWest[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftCosEvtAveEast[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftSinEvtAveEast[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftCosEvtAveWest[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
    TProfile2D *p_mEpdQ1GrpShiftSinEvtAveWest[mNumVzBin][mNumShiftCorr][mNumRingsGrps];

    // Shift Correction for Full EP
    TProfile2D *p_mEpdQ1SideShiftCosFull[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mEpdQ1SideShiftSinFull[mNumVzBin][mNumShiftCorr];

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
    TProfile *p_mEpdSubEpSideV1[mNumCentrality]; // v1 vs. eta

    // Event Plane Distribution
    TH2F *h_mEpdEp1SideRawEast[mNumCentrality]; // 1st raw EP
    TH2F *h_mEpdEp1SideRawWest[mNumCentrality];
    TH2F *h_mEpdEp1SideRawFull[mNumCentrality];
    TH2F *h_mEpdEp1SideRawCorr[mNumCentrality]; // Psi1East vs Psi1West
    TH2F *h_mEpdEp1GrpRawEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpRawWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpRawFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpRawCorr[mNumCentrality][mNumRingsGrps];

    TH2F *h_mEpdEp1SideWgtEast[mNumCentrality]; // 1st weighted EP
    TH2F *h_mEpdEp1SideWgtWest[mNumCentrality];
    TH2F *h_mEpdEp1SideWgtFull[mNumCentrality];
    TH2F *h_mEpdEp1SideWgtCorr[mNumCentrality]; // Psi1East vs Psi1West
    TH2F *h_mEpdEp1GrpWgtEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpWgtWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpWgtFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpWgtCorr[mNumCentrality][mNumRingsGrps];

    TH2F *h_mEpdEp1SideReCtrEast[mNumCentrality]; // 1st recenter EP
    TH2F *h_mEpdEp1SideReCtrWest[mNumCentrality];
    TH2F *h_mEpdEp1SideReCtrFull[mNumCentrality];
    TH2F *h_mEpdEp1SideReCtrCorr[mNumCentrality]; // Psi1East vs Psi1West
    TH2F *h_mEpdEp1GrpReCtrTrkAveEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrTrkAveWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrTrkAveFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrTrkAveCorr[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrEvtAveEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrEvtAveWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrEvtAveFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpReCtrEvtAveCorr[mNumCentrality][mNumRingsGrps];

    TH2F *h_mEpdEp1SideShiftEast[mNumCentrality]; // 1st shift EP
    TH2F *h_mEpdEp1SideShiftWest[mNumCentrality];
    TH2F *h_mEpdEp1SideShiftFull[mNumCentrality];
    TH2F *h_mEpdEp1SideShiftCorr[mNumCentrality]; // Psi1East vs Psi1West
    TH2F *h_mEpdEp1GrpShiftTrkAveEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftTrkAveWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftTrkAveFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftTrkAveCorr[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftEvtAveEast[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftEvtAveWest[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftEvtAveFull[mNumCentrality][mNumRingsGrps];
    TH2F *h_mEpdEp1GrpShiftEvtAveCorr[mNumCentrality][mNumRingsGrps];

    TH2F *h_mEpdEp1SideShiftFullCorr[mNumCentrality]; // 1st shift full EP

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
