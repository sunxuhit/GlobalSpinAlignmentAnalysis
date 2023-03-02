#ifndef StZdcEpManager_h
#define StZdcEpManager_h

#include <string>
#include "TObject.h"
#include "TVector2.h"

class TFile;
class TProfile2D;
class TProfile;
class TH2F;

class StZdcEpManager : public TObject
{
  public:
    StZdcEpManager(int beamType);
    virtual ~StZdcEpManager();
    void clearZdcEpManager();
    void initZdcEpManager(int cent9, int runIndex, int vzBin);

    void setZdcSmd(int eastwest,int verthori,int strip,const double zdcsmd);
    double getZdcSmd(int eastwest,int verthori,int strip);

    // Gain Correction
    void initZdcGain();
    void fillZdcGain(int i_eastwest, int i_verthori, int i_slat, double zdcsmd);
    void writeZdcGain();
    void readZdcGain();
    void setZdcSmdGainCorr(int eastwest,int verthori,int strip,const double zdcsmd);
    double getZdcSmdGainCorr(int eastwest,int verthori,int strip);

    // ReCenter Correction
    void initZdcReCtr();
    void fillZdcReCtrEast(TVector2 QVector);
    void fillZdcReCtrWest(TVector2 QVector);
    void writeZdcReCtr();
    void readZdcReCtr();
    void setZdcSmdCenter();

    // Shift Correction
    void initZdcShift(); // East/West
    void fillZdcShiftEast(TVector2 QVector);
    void fillZdcShiftWest(TVector2 QVector);
    void writeZdcShift();
    void readZdcShift();
    TVector2 applyZdcSmdShiftCorrEast(TVector2 QVector);
    TVector2 applyZdcSmdShiftCorrWest(TVector2 QVector);
    double transPsi1(double Psi1);
    bool isPsi1InRange(double Psi1);

    void initZdcShiftFull(); // Full
    void fillZdcShiftFull(TVector2 QVector);
    void writeZdcShiftFull();
    void readZdcShiftFull();
    TVector2 applyZdcSmdShiftCorrFull(TVector2 QVector);

    // Event Plane Resolution
    void initZdcResolution(); // Full
    void fillZdcResolution(TVector2 QEast, TVector2 QWest);
    void writeZdcResolution();
    void readZdcResolution();
    double getZdcSubEpResVal(int cent9);
    double getZdcSubEpResErr(int cent9);
    double getZdcFullEpResVal(int cent9);
    double getZdcFullEpResErr(int cent9);

    // Charged Hadron Directed Flow
    void initZdcFullEpFlow(); // Full EP
    void fillZdcFullEpDFlow(double eta, double pt, double v1, double reweight);
    void writeZdcFullEpFlow();

    // QVector
    double getPosition(int eastwest,int verthori,int strip, int mode);
    TVector2 getQ1VecEast(int mode);
    TVector2 getQ1VecWest(int mode);
    TVector2 getQ1VecFull(TVector2 QEast, TVector2 QWest, int mode);

    // Event Plane Distribution
    void initZdcSubEpRaw(); // raw Sub EP
    void fillZdcSubEpRaw(TVector2 QEast, TVector2 QWest, TVector2 QFull);
    void writeZdcSubEpRaw();

    void initZdcSubEpReCtr(); // recenter Sub EP
    void fillZdcSubEpReCtr(TVector2 QEast, TVector2 QWest, TVector2 QFull);
    void writeZdcSubEpReCtr();

    void initZdcSubEpShift(); // shift Sub EP
    void fillZdcSubEpShift(TVector2 QEast, TVector2 QWest, TVector2 QFull);
    void writeZdcSubEpShift();

    void initZdcFullEpShift(); // shift Full EP
    void fillZdcFullEpShift(TVector2 QFull);
    void writeZdcFullEpShift();

  private:
    static const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
    static const int mNumShiftCorr = 20;
    static const int mNumCentrality = 9; // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%

    int mCent9;
    int mRunIndex;
    int mVzBin;
    double mZdcSmd[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y);
    double mGainFactor[2][2][8];

    // Gain Correction
    TH2F *h_mZdcGainCorr[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y); | x-axis: runIndex | y-axis: ADC

    // ReCenter Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mZdcQ1ReCtrVertEast[mNumVzBin];
    TProfile2D *p_mZdcQ1ReCtrHoriEast[mNumVzBin];
    TProfile2D *p_mZdcQ1ReCtrVertWest[mNumVzBin];
    TProfile2D *p_mZdcQ1ReCtrHoriWest[mNumVzBin];
    double mCenterVertEast, mCenterHoriEast, mCenterVertWest, mCenterHoriWest;

    // Shift Correction for East/West EP
    TProfile2D *p_mZdcQ1ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mZdcQ1ShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mZdcQ1ShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mZdcQ1ShiftSinWest[mNumVzBin][mNumShiftCorr];

    // Shift Correction for Full EP
    TProfile2D *p_mZdcQ1ShiftCosFull[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mZdcQ1ShiftSinFull[mNumVzBin][mNumShiftCorr];

    // Event Plane Resolution
    TProfile *p_mZdcSubEp1Res;
    double mZdcSubEp1ResVal[mNumCentrality];
    double mZdcSubEp1ResErr[mNumCentrality];
    double mZdcFullEp1ResVal[mNumCentrality];
    double mZdcFullEp1ResErr[mNumCentrality];

    // Charged Hadron Directed Flow
    TProfile *p_mZdcFullEpDFlow[mNumCentrality]; // v1 vs. eta

    // Event Plane Distribution | x axis is runIndex, y axis is EP angle
    TH2F *h_mZdcEp1RawEast[mNumCentrality]; // raw EP
    TH2F *h_mZdcEp1RawWest[mNumCentrality];
    TH2F *h_mZdcEp1RawFull[mNumCentrality]; // Qwest-QEast
    TH2F *h_mZdcEp1RawCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mZdcEp1ReCtrEast[mNumCentrality]; // recenter EP
    TH2F *h_mZdcEp1ReCtrWest[mNumCentrality];
    TH2F *h_mZdcEp1ReCtrFull[mNumCentrality]; // Qwest-QEast
    TH2F *h_mZdcEp1ReCtrCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mZdcEp1ShiftEast[mNumCentrality]; // shift EP
    TH2F *h_mZdcEp1ShiftWest[mNumCentrality];
    TH2F *h_mZdcEp1ShiftFull[mNumCentrality]; // Qwest-QEast
    TH2F *h_mZdcEp1ShiftCorr[mNumCentrality]; // Psi1East vs Psi1West

    TH2F *h_mZdcEp1ShiftFullCorr[mNumCentrality]; // Qwest-QEast

    TFile *file_mGainCorrPar;
    TFile *file_mReCtrPar;
    TFile *file_mShiftPar;
    TFile *file_mShiftParFull;
    TFile *file_mResolution;

    std::string str_mEastWest[2] = {"East","West"};
    std::string str_mVertHori[2] = {"Vert","Hori"};

    const int mType;

  ClassDef(StZdcEpManager,1)
};

#endif
