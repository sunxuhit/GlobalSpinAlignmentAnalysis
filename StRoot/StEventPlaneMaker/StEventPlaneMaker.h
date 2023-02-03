#ifndef StEventPlaneMaker_h
#define StEventPlaneMaker_h

#include "StMaker.h"
// #include "TString.h"
#include <string>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;

class StEventPlaneUtility;
class StEventPlaneCut;
class StEventPlaneHistoManager;
class StEventPlaneProManager;
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

    StEventPlaneUtility *mEventPlaneUtility;
    StEventPlaneCut *mEventPlaneCut;
    StEventPlaneHistoManager *mEventPlaneHistoManager;
    StEventPlaneProManager *mEventPlaneProManager;
    StZdcEpManager *mZdcEpManager;
    // StEventPlaneBbcEpManager *mEventPlaneBbcEpManager;
    // StEventPlaneEpdEpManager *mEventPlaneEpdEpManager;
    // StEventPlaneTpcEpManager *mEventPlaneTpcEpManager;
    
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
