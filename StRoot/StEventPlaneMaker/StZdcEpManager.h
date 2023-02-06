#ifndef StZdcEpManager_h
#define StZdcEpManager_h

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
    void clearZdcEp();
    void initZdcEp(int Cent9, int RunIndex, int vzBin);

    void setZdcSmd(int eastwest,int verthori,int strip,const double zdcsmd);
    double getZdcSmd(int eastwest,int verthori,int strip);

    // Gain Correction
    void initZdcGainCorr();
    void fillZdcGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, double zdcsmd);
    void writeZdcGainCorr();

    void readGainCorr();
    void setZdcSmdGainCorr(int eastwest,int verthori,int strip,const double zdcsmd);
    double getZdcSmdGainCorr(int eastwest,int verthori,int strip);
    double getPosition(int eastwest,int verthori,int strip, int mode);
    TVector2 getQEast(int mode);
    TVector2 getQWest(int mode);

    void initZdcRawEP(); // raw EP
    void fillZdcRawEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex);
    void writeZdcRawEP();

    // ReCenter Correction
    void initZdcReCenter();
    void fillZdcReCenterEast(TVector2 qVector, int Cent9, int RunIndex, int vzBin);
    void fillZdcReCenterWest(TVector2 qVector, int Cent9, int RunIndex, int vzBin);
    void writeZdcReCenter();

    void readReCenterCorr();
    void setZdcSmdCenter();

    void readShiftCorr();
    TVector2 applyZdcSmdShiftCorrEast(TVector2 qVector);
    TVector2 applyZdcSmdShiftCorrWest(TVector2 qVector);
    double angleShift(double Psi_shifted);

    // Shift Correction
    void readShiftCorrFull();
    TVector2 applyZdcSmdShiftCorrFull(TVector2 qVector);
    TVector2 getQFull(TVector2 QEast, TVector2 QWest);

    void readResolution();
    void calResolution();
    double getResolution(int Cent9);

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
    TH2F *h_mZdcRawEpEast[9]; // raw EP
    TH2F *h_mZdcRawEpWest[9];
    TH2F *h_mZdcRawEpFull[9]; // Qwest-QEast

    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mZdcQEastVertical[mNumVzBin];
    TProfile2D *p_mZdcQEastHorizontal[mNumVzBin];
    TProfile2D *p_mZdcQWestVertical[mNumVzBin];
    TProfile2D *p_mZdcQWestHorizontal[mNumVzBin];
    double mCenterEastVertical, mCenterEastHorizontal, mCenterWestVertical, mCenterWestHorizontal;

    // Shift Correction for East/West
    TProfile2D *p_mZdcQEastCos[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mZdcQEastSin[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mZdcQWestCos[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mZdcQWestSin[mNumVzBin][mNumShiftCorr];

    // Shift Correction for East/West
    TProfile2D *p_mZdcQFullCos[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mZdcQFullSin[mNumVzBin][mNumShiftCorr];

    // charged hadron v1 calculation
    TProfile *p_mZdcEpResolution;
    double mResolution[9];

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
