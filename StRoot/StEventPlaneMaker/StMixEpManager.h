#ifndef StMixEpManager_h
#define StMixEpManager_h

#include "TObject.h"

class TFile;
class TProfile;

class StMixEpManager : public TObject
{
  public:
    StMixEpManager(int beamType);
    virtual ~StMixEpManager();
    void clearMixEpManager();
    void initMixEpManager(int cent9, int runIndex, int vzBin);

    // Event Plane Resolution
    void initMixEpRes();
    void fillMixEpRes(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest);
    void writeMixEpRes();
    void readMixEpRes();
    double propMixEpResErr(double valA, double sigA, double valB, double sigB, double valC, double sigC); // return the error of valSubResA*valSubResB/valSubResC
    double getMixSubEp1ResVal(int cent9, int grpId);
    double getMixSubEp1ResErr(int cent9, int grpId);

    // Deuteron Directed Flow
    void initMixSubEpFlow(); // Sub EP
    void fillMixSubEpDeuV1(double rap, double v1, double reweight);
    void writeMixSubEpFlow();

    // Event Plane Distribution
    void initMixSubEpRaw();
    void fillMixSubEpRaw(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest);
    void writeMixSubEpRaw();

    void initMixSubEpReCtr();
    void fillMixSubEpReCtr(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest);
    void writeMixSubEpReCtr();

    void initMixSubEpShift();
    void fillMixSubEpShift(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest);
    void writeMixSubEpShift();

  private:
    static const int mNumCentrality = 9; // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
    static const int mNumEpGroup    = 6;  

    // EPD EP Resolution
    // 0: EpdEpGrp0 vs. TpcEpEast | 1: EpdEpGrp0 vs. TpcEpWest | 2: EpdEpGrp1 vs. TpcEpEast
    // 3: EpdEpGrp1 vs. TpcEpWest | 4: EpdEpGrp0 vs. EpdEpGrp1 | 5: TpcEpEast vs. TpcEpWest
    TProfile *p_mMixSubEp1Res[mNumEpGroup]; // 1st EP
    // 0: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpWest (default) | 1: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpEast | 2: EpdEpGrp0 vs. TpcEpEast && TpcEpWest
    // 3: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpWest (mainSys) | 4: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpEast | 5: EpdEpGrp1 vs. TpcEpEast && TpcEpWest
    double mMixSubEp1ResVal[mNumCentrality][mNumEpGroup];
    double mMixSubEp1ResErr[mNumCentrality][mNumEpGroup];

    // deutron Directed Flow
    TProfile *p_mMixSubEpDeuV1; // deutron v1(EpdEpGrp0) vs. rap

    // Event Plane Distribution
    // 0: EpdEpGrp0 vs. TpcEpEast | 1: EpdEpGrp0 vs. TpcEpWest | 2: EpdEpGrp1 vs. TpcEpEast
    // 3: EpdEpGrp1 vs. TpcEpWest | 4: EpdEpGrp0 vs. EpdEpGrp1 | 5: TpcEpEast vs. TpcEpWest
    TH2F *h_mMixEp1RawCorr[mNumCentrality][mNumEpGroup]; 
    TH2F *h_mMixEp1ReCtrCorr[mNumCentrality][mNumEpGroup]; 
    TH2F *h_mMixEp1ShiftCorr[mNumCentrality][mNumEpGroup]; 

    TFile *file_mResolution;

    int mCent9, mRunIndex, mVzBin;
    const int mType;

  ClassDef(StMixEpManager,1)
};

#endif
