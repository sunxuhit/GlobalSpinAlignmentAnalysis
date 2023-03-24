#ifndef StEventPlaneMaker_h
#define StEventPlaneMaker_h

#include "StMaker.h"
#include <string>

class TProfile;

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StPileupUtil;

class StAnalysisUtils;
class StAnalysisCut;
class StZdcEpManager;
class StEpdEpManager;
class StTpcEpManager;
class StMixEpManager;

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
    StPileupUtil   *mPileupUtilFxt3p85; // from Xionghong He (Developed by Sooraj Radhakrishnan & only for Fxt3p85GeV_2018)

    StAnalysisUtils *mAnaUtils; // Analysis Utilities
    StAnalysisCut   *mAnaCut;   // Analysis Cuts
    StZdcEpManager  *mZdcEpManager;
    StEpdEpManager  *mEpdEpManager;
    StTpcEpManager  *mTpcEpManager;
    StMixEpManager  *mMixEpManager; // for FXT ONLY
    
    string str_mOutPutGainCorr;
    string str_mOutPutReCenterPar;
    string str_mOutPutShiftPar;
    string str_mOutPutResolution;
    string str_mOutPutFlow;

    TFile *file_mOutPutGainCorr;
    TFile *file_mOutPutReCenterPar;
    TFile *file_mOutPutShiftPar;
    TFile *file_mOutPutResolution;
    TFile *file_mOutPutFlow;

    static const int mNumRingsGrps = 2;  // EPD Group 0: 0-7 rings | Group 1: 8-15 rings
    const int mMode;
    const int mType;

    ClassDef(StEventPlaneMaker, 1)
};

#endif
