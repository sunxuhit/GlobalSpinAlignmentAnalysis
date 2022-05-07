#ifndef StRunQAMaker_h
#define StRunQAMaker_h

#include "StMaker.h"
#include <string>

class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StRefMultCorr;

class StAnalysisUtils;
class StAnalysisCut;
class StRunQAHistoManager;
class StRunQAProManager;

class StRunQAMaker : public StMaker {
  public:
    StRunQAMaker(const char *name, StPicoDstMaker *picoMaker, const string jobId, const int beamType);
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

    StAnalysisUtils     *mAnaUtils;
    StAnalysisCut       *mAnaCut;
    StRunQAHistoManager *mRunQAHistoManager;
    StRunQAProManager   *mRunQAProManager;
    
    int mType;

    string mOutPut_RunQA;
    TFile *mFile_QA;

    ClassDef(StRunQAMaker, 1)
};

#endif
