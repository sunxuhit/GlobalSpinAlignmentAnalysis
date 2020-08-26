#ifndef StVecMesonMaker_h
#define StVecMesonMaker_h

#include "StMaker.h"
// #include "TString.h"
#include <string>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StVecMesonCut;
class StVecMesonHistoManager;
class StVecMesonUtility;
class StVecMesonProManager;
class StVecMesonZdcEpManager;
// class StVecMesonTree;

class StVecMesonMaker : public StMaker {
  public:
    StVecMesonMaker(const char *name, StPicoDstMaker *picoMaker, const string jobId, const int Mode, const int Energy, const int Flag_ME);
    virtual ~StVecMesonMaker();
    
    virtual int Init();
    virtual int Make();
    virtual void  Clear(Option_t *opt="");
    virtual int Finish();
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent    *mPicoEvent;
    StRefMultCorr  *mRefMultCorr;
    StVecMesonCut  *mVecMesonCut;
    StVecMesonHistoManager *mVecMesonHistoManager;
    StVecMesonUtility *mVecMesonUtility;
    StVecMesonProManager *mVecMesonProManager;
    StVecMesonZdcEpManager *mVecMesonZdcEpManager;
    // StVecMesonBbcEpManager *mVecMesonBbcEpManager;
    // StVecMesonEpdEpManager *mVecMesonEpdEpManager;
    // StVecMesonTpcEpManager *mVecMesonTpcEpManager;
    // StVecMesonTree *mVecMesonTree;
    
    int mMode;
    int mEnergy;
    int mFlag_ME;

    string mInPut_Corr_ReCenter;

    string mOutPut_QA;
    string mOutPut_ZdcGainCorr;
    // string mOutPut_BbcGainCorr;
    string mOutPut_ReCenterPar;
    string mOutPut_ShiftPar;
    string mOutPut_Resolution;
    string mOutPut_Phi;

    TFile *mFile_QA;
    TFile *mFile_ZdcGainCorr;
    // TFile *mFile_BbcGainCorr;
    TFile *mFile_ReCenterPar;
    TFile *mFile_ShiftPar;
    TFile *mFile_Resolution;
    TFile *mFile_Phi;

    int mUsedTrackCounter;

    ClassDef(StVecMesonMaker, 1)
};

#endif
