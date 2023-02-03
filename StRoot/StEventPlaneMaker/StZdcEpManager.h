#ifndef StZdcEpManager_h
#define StZdcEpManager_h

#include "TObject.h"
#include "TVector2.h"
// #include "TString.h"

class TProfile2D;
class TProfile;
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

    void readGainCorr();
    void setZdcSmdGainCorr(int eastwest,int verthori,int strip,const double zdcsmd);
    double getZdcSmdGainCorr(int eastwest,int verthori,int strip);
    double getPosition(int eastwest,int verthori,int strip, int mode);
    TVector2 getQEast(int mode);
    TVector2 getQWest(int mode);

    void readReCenterCorr();
    void setZdcSmdCenter();

#if 0
    void readShiftCorr();
    TVector2 ApplyZdcSmdShiftCorrEast(TVector2 qVector);
    TVector2 ApplyZdcSmdShiftCorrWest(TVector2 qVector);
    double AngleShift(double Psi_shifted);

    void readShiftCorrFull();
    TVector2 ApplyZdcSmdShiftCorrFull(TVector2 qVector);
    TVector2 getQFull(TVector2 QEast, TVector2 QWest);

    void readResolution();
    void calResolution();
    double getResolution(int Cent9);
#endif

  private:
    int mCent9;
    int mRunIndex;
    int mVzBin;
    double mZdcSmd[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y);
    double mGainCorrFactor[2][2][8];

    static const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mZdcQEastVertical[mNumVzBin];
    TProfile2D *p_mZdcQEastHorizontal[mNumVzBin];
    TProfile2D *p_mZdcQWestVertical[mNumVzBin];
    TProfile2D *p_mZdcQWestHorizontal[mNumVzBin];
    double mCenterEastVertical, mCenterEastHorizontal, mCenterWestVertical, mCenterWestHorizontal;

#if 0
    static const int mNumShiftCorr = 20;
    // Shift Correction for East/West
    TProfile2D *p_mQEastCos[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mQEastSin[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mQWestCos[mNumVzBin][mNumShiftCorr];
    TProfile2D *p_mQWestSin[mNumVzBin][mNumShiftCorr];

    // Shift Correction for East/West
    TProfile2D *p_mQFullCos[[mNumVzBin]mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
    TProfile2D *p_mQFullSin[mNumVzBin][mNumShiftCorr];

    // charged hadron v1 calculation
    TProfile *p_mResolution;
    double mResolution[9];
#endif

    TFile *file_mGainCorrPar;
    TFile *file_mReCenterPar;
    // TFile *file_mShiftPar;
    // TFile *file_mShiftParFull;
    // TFile *file_mResolution;

    std::string str_mEastWest[2] = {"East","West"};
    std::string str_mVertHori[2] = {"Vert","Hori"};

    const int mType;

  ClassDef(StZdcEpManager,1)
};

#endif
