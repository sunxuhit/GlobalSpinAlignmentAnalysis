#ifndef StTpcEpManager_h
#define StTpcEpManager_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class TProfile2D;
class TFile;
class TNtuple;

class StTpcEpManager : public TObject
{
  public:
    StTpcEpManager(Int_t energy);
    virtual ~StTpcEpManager();


    // ReCenter Correction
    bool passTrackEtaEast(StPicoTrack*);
    bool passTrackEtaWest(StPicoTrack*);
    bool passTrackFull(StPicoTrack*);

    TVector2 calq2Vector(StPicoTrack*);
    Float_t getWeight(StPicoTrack*);

    void InitReCenterCorrection();
    void addTrack_EastRaw(StPicoTrack* track, Int_t Cent9, Int_t RunIndex);
    void addTrack_WestRaw(StPicoTrack* track, Int_t Cent9, Int_t RunIndex);
    void addTrack_FullRaw(StPicoTrack* track, Int_t Cent9, Int_t RunIndex);

    void addTrack_East(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i); // i = vz_sign
    void addTrack_West(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i);
    void addTrack_Full(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i);

    void addTrack_A(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i); // i = vz_sign || random sub A
    void addTrack_B(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i); // i = vz_sign || random sub B
    void Randomization();

    TVector2 getReCenterPar_East(Int_t Cent9, Int_t RunIndex, Int_t vz_sign);
    TVector2 getReCenterPar_West(Int_t Cent9, Int_t RunIndex, Int_t vz_sign);
    TVector2 getReCenterPar_Full(Int_t Cent9, Int_t RunIndex, Int_t vz_sign);

    void print(TVector2);
    void clear();

    // Shift Correction
    bool passTrackEtaNumCut();
    bool passTrackFullNumCut();
    bool passTrackEtaNumRawCut();
    bool passTrackFullNumRawCut();

    // Event Plane method
    TVector2 calPsi2_East_EP(Int_t); // 0 = ShiftOrder: 2, 4, 6, 8, 10
    TVector2 calPsi2_West_EP(Int_t);
    TVector2 calPsi2_Full_EP(Int_t);

    void InitShiftCorrection();
    Float_t AngleShift(Float_t Psi_raw);

    // Event Plane method
    Float_t calShiftAngle2East_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign);
    Float_t calShiftAngle2West_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign);
    Float_t calShiftAngle2A_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign);
    Float_t calShiftAngle2B_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign);
    Float_t calShiftAngle2Full_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign);
    Float_t calShiftAngle2Full_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track); // subtract self-correlation

    void InitResolutionCorr();
    Float_t getResolution2_EP(Int_t Cent9);
    Float_t getResolution2_Full_EP(Int_t Cent9);

    TVector2 getQVector(Int_t l); // east/west
    TVector2 getQVectorRaw(Int_t l);
    Int_t getNumTrack(Int_t);

  private:
    //Event Plane method
    TVector2 mQ2Vector_EastRaw_EP, mQ2Vector_WestRaw_EP, mQ2Vector_FullRaw_EP;
    TVector2 mQ2Vector_East_EP, mQ2Vector_West_EP, mQ2Vector_Full_EP, mQ2Vector_A_EP, mQ2Vector_B_EP;

    Int_t    mQCounter_RawEast, mQCounter_RawWest, mQCounter_RawFull;
    Int_t    mQCounter_East, mQCounter_West, mQCounter_Full;
    Int_t    mQCounter_Full_East, mQCounter_Full_West;
    Int_t    mQCounter_A, mQCounter_B;
    Int_t    mEnergy;

    TFile *mInPutFile;
    TFile *mInPutFile_Shift;
    TFile *mInPutFile_Res;

    static TString mVStr[2];
    static TString mOrder;

  ClassDef(StTpcEpManager,1)
};

#endif
