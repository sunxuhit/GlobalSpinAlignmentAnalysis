#ifndef StZdcEpManager_h
#define StZdcEpManager_h

#include <string>
#include "TObject.h"
#include "TVector2.h"

class TProfile2D;
class TProfile;
class TH2F;
class TFile;

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
    void initZdcGainCorr();
    void fillZdcGainCorr(int i_eastwest, int i_verthori, int i_slat, double zdcsmd);
    void writeZdcGainCorr();
    void readZdcGainCorr();
    void setZdcSmdGainCorr(int eastwest,int verthori,int strip,const double zdcsmd);
    double getZdcSmdGainCorr(int eastwest,int verthori,int strip);

    // ReCenter Correction
    void initZdcReCenter();
    void fillZdcReCenterEast(TVector2 QVector);
    void fillZdcReCenterWest(TVector2 QVector);
    void writeZdcReCenter();
    void readZdcReCenterCorr();
    void setZdcSmdCenter();

    // Shift Correction
    void initZdcShift(); // East/West
    void fillZdcShiftEast(TVector2 QVector);
    void fillZdcShiftWest(TVector2 QVector);
    void writeZdcShift();
    void readZdcShiftCorr();
    TVector2 applyZdcSmdShiftCorrEast(TVector2 QVector);
    TVector2 applyZdcSmdShiftCorrWest(TVector2 QVector);
    double angleShift(double Psi_shifted);

    void initZdcShiftFull(); // Full
    void fillZdcShiftFull(TVector2 QVector);
    void writeZdcShiftFull();
    void readZdcShiftCorrFull();
    TVector2 applyZdcSmdShiftCorrFull(TVector2 QVector);

    // QVector
    double getPosition(int eastwest,int verthori,int strip, int mode);
    TVector2 getQEast(int mode);
    TVector2 getQWest(int mode);
    TVector2 getQFull(TVector2 QEast, TVector2 QWest, int mode);

    // Event Plane Resolution
    void readZdcResolution();
    void calZdcResolution();
    double getZdcResFullVal(int cent9);
    double getZdcResFullErr(int cent9);

    // Event Plane Distribution
    void initZdcRawEP(); // raw EP
    void fillZdcRawEP(TVector2 QEast, TVector2 QWest, TVector2 QFull);
    void writeZdcRawEP();

  private:
    static const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
    static const int mNumShiftCorr = 20;

    int mCent9;
    int mRunIndex;
    int mVzBin;
    double mZdcSmd[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y);
    double mGainFactor[2][2][8];

    // Gain Correction
    TH2F *h_mZdcGainCorr[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y); | x-axis: runIndex | y-axis: ADC

    // ReCenter Correction | x axis is runIndex, y axis is Centrality
    TProfile2D *p_mZdcQReCenterVertEast[mNumVzBin];
    TProfile2D *p_mZdcQReCenterHoriEast[mNumVzBin];
    TProfile2D *p_mZdcQReCenterVertWest[mNumVzBin];
    TProfile2D *p_mZdcQReCenterHoriWest[mNumVzBin];
    double mCenterVertEast, mCenterHoriEast, mCenterVertWest, mCenterHoriWest;

    // Shift Correction for East/West EP
    TProfile2D *p_mZdcQShiftCosEast[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mZdcQShiftSinEast[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mZdcQShiftCosWest[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mZdcQShiftSinWest[mNumVzBin][mNumShiftCorr];

    // Shift Correction for Full EP
    TProfile2D *p_mZdcQShiftCosFull[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mZdcQShiftSinFull[mNumVzBin][mNumShiftCorr];

    // charged hadron v1 calculation
    TProfile *p_mZdcEpResolution;
    double mZdcResFullVal[9];
    double mZdcResFullErr[9];

    // Event Plane Distribution
    TH2F *h_mZdcRawEpEast[9]; // raw EP
    TH2F *h_mZdcRawEpWest[9];
    TH2F *h_mZdcRawEpFull[9]; // Qwest-QEast

    TFile *file_mGainCorrPar;
    TFile *file_mReCenterPar;
    TFile *file_mShiftPar;
    TFile *file_mShiftParFull;
    TFile *file_mResolution;

    std::string str_mEastWest[2] = {"East","West"};
    std::string str_mVertHori[2] = {"Vert","Hori"};

    const int mType;

  ClassDef(StZdcEpManager,1)
};

#endif
