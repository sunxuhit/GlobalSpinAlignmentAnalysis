#ifndef StEventPlaneMaker_h
#define StEventPlaneMaker_h

#include "StMaker.h"
#include <string>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;

class StAnalysisUtils;
class StAnalysisCut;
class StZdcEpManager;

class StEventPlaneMaker : public StMaker {
  public:
    StEventPlaneMaker(const char *name, StPicoDstMaker *picoMaker, const string jobId, const int mode, const int beamType);
    virtual ~StEventPlaneMaker();
    
    virtual int Init();
    virtual int Make();
    virtual void Clear(Option_t *opt="");
    virtual int Finish();
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent    *mPicoEvent;
    StRefMultCorr  *mRefMultCorr;

    StAnalysisUtils *mAnaUtils; // Analysis Utilities
    StAnalysisCut   *mAnaCut; // Analysis Cuts
    StZdcEpManager  *mZdcEpManager;
    // StEventPlaneTpcEpManager *mTpcEpManager;
    // StEventPlaneBbcEpManager *mBbcEpManager;
    // StEventPlaneEpdEpManager *mEpdEpManager;
    
    const int mMode;
    const int mType;

    string str_mOutPutGainCorr;
    string str_mOutPutReCenterPar;
    // string str_mOutPutShiftPar;
    // string str_mOutPutResolution;

    TFile *file_mOutPutGainCorr;
    TFile *file_mOutPutReCenterPar;
    // TFile *file_mOutPutShiftPar;
    // TFile *file_mOutPutResolution;

    int mUsedTrackCounter;

    ClassDef(StEventPlaneMaker, 1)
};

#endif
