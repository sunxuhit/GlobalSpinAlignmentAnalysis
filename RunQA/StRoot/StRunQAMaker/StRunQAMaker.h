#ifndef StRunQAMaker_h
#define StRunQAMaker_h

#include "StMaker.h"
// #include "TString.h"
#include <string>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StRunQACut;
class StRunQAHistoManager;
class StRunQAUtility;
class StRunQAProManager;

class StRunQAMaker : public StMaker {
  public:
    StRunQAMaker(const char *name, StPicoDstMaker *picoMaker, const string jobId, const int Mode, const int Energy);
    virtual ~StRunQAMaker();
    
    virtual int Init();
    virtual int Make();
    virtual void  Clear(Option_t *opt="");
    virtual int Finish();
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent    *mPicoEvent;
    StRefMultCorr  *mRefMultCorr;
    StRunQACut  *mRunQACut;
    StRunQAHistoManager *mRunQAHistoManager;
    StRunQAUtility *mRunQAUtility;
    StRunQAProManager *mRunQAProManager;
    
    int mMode;
    int mEnergy;

    string mInPut_Corr_ReCenter;

    string mOutPut_RunQA;

    TFile *mFile_QA;

    int mUsedTrackCounter;

    ClassDef(StRunQAMaker, 1)
};

#endif
