#ifndef StZdcEpManager_h
#define StZdcEpManager_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"

class TProfile2D;
class TProfile;
class TFile;

class StZdcEpManager : public TObject
{
  public:
    StZdcEpManager(int energy);
    virtual ~StZdcEpManager();
    void clearZdcEp();
    void initZdcEp(int Cent9, int RunIndex);

    void setZdcSmd(int eastwest,int verthori,int strip,const float zdcsmd);
    float getZdcSmd(int eastwest,int verthori,int strip);

    void readGainCorr();
    void setZdcSmdGainCorr(int eastwest,int verthori,int strip,const float zdcsmd);
    float getZdcSmdGainCorr(int eastwest,int verthori,int strip);
    float getPosition(int eastwest,int verthori,int strip, int mode);
    TVector2 getQEast(int mode);
    TVector2 getQWest(int mode);

    void readReCenterCorr();
    void setZdcSmdCenter();

    void readShiftCorr();
    TVector2 ApplyZdcSmdShiftCorrEast(TVector2 qVector);
    TVector2 ApplyZdcSmdShiftCorrWest(TVector2 qVector);
    float AngleShift(float Psi_shifted);

    void readShiftCorrFull();
    TVector2 ApplyZdcSmdShiftCorrFull(TVector2 qVector);
    TVector2 getQFull(TVector2 QEast, TVector2 QWest);

    void readResolution();
    void calResolution();
    float getResolution(int Cent9);

  private:

    int mEnergy;
    int mCent9;
    int mRunIndex;
    int mVz_sign;
    float mZdcSmd[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y);
    float mGainCorrFactor[2][2][8];

    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mQEastVertical; // vz_sign
    TProfile2D *p_mQEastHorizontal;
    TProfile2D *p_mQWestVertical;
    TProfile2D *p_mQWestHorizontal;
    float mCenterEastVertical, mCenterEastHorizontal, mCenterWestVertical, mCenterWestHorizontal;

    // Shift Correction for East/West
    TProfile2D *p_mQEastCos[20]; // 0 = vertex pos/neg | 1 = shift correction harmonics
    TProfile2D *p_mQEastSin[20];
    TProfile2D *p_mQWestCos[20];
    TProfile2D *p_mQWestSin[20];

    // Shift Correction for East/West
    TProfile2D *p_mQFullCos[20]; // 0 = vertex pos/neg | 1 = shift correction harmonics
    TProfile2D *p_mQFullSin[20];

    // charged hadron v1 calculation
    TProfile *p_mResolution;
    float mResolution[9];

    TFile *mFile_GainCorrPar;
    TFile *mFile_ReCenterPar;
    TFile *mFile_ShiftPar;
    TFile *mFile_ShiftParFull;
    TFile *mFile_Resolution;

  ClassDef(StZdcEpManager,1)
};

#endif
