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
    StRunQAMaker(const char *name, StPicoDstMaker *picoMaker, string jobId, int beamType);
    virtual ~StRunQAMaker();
    
    virtual int Init();
    virtual int Make();
    virtual void  Clear(Option_t *opt="");
    virtual int Finish();
    
  private:
    StPicoDstMaker *mPicoDstMaker; // STAR Lib
    StPicoDst      *mPicoDst;
    StPicoEvent    *mPicoEvent;
    StRefMultCorr  *mRefMultCorr;
    StPileupUtil   *mPileupUtilFxt3p85; // from Xionghong He (Developed by Sooraj Radhakrishnan & only for Fxt3p85GeV_2018)

    StAnalysisUtils     *mAnaUtils; // Analysis Module
    StAnalysisCut       *mAnaCut;
    StRunQAHistoManager *mRunQAHistoManager;
    StRunQAProManager   *mRunQAProManager;
    
    const int mType;

    std::string str_mOutPutRunQA;
    TFile *file_mOutPutRunQA;

    ClassDef(StRunQAMaker, 1)
};

#endif
