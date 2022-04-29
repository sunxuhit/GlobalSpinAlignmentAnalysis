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
    StEventPlaneMaker(const char *name, StPicoDstMaker *picoMaker, const string jobId, const int Mode, const int Energy);
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
    StEventPlaneCut  *mEventPlaneCut;
    StEventPlaneHistoManager *mEventPlaneHistoManager;
    StEventPlaneProManager *mEventPlaneProManager;
    StZdcEpManager *mZdcEpManager;
    // StEventPlaneBbcEpManager *mEventPlaneBbcEpManager;
    // StEventPlaneEpdEpManager *mEventPlaneEpdEpManager;
    // StEventPlaneTpcEpManager *mEventPlaneTpcEpManager;
    
    int mMode;
    int mEnergy;

    string mOutPut_GainCorr;
    string mOutPut_ReCenterPar;
    // string mOutPut_ShiftPar;
    // string mOutPut_Resolution;

    TFile *mFile_GainCorr;
    TFile *mFile_ReCenterPar;
    // TFile *mFile_ShiftPar;
    // TFile *mFile_Resolution;

    int mUsedTrackCounter;

    ClassDef(StEventPlaneMaker, 1)
};

#endif
