#ifndef StPhiMesonMaker_h
#define StPhiMesonMaker_h

#include "StMaker.h"
#include <string>

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
// class StMixEpManager;
class StPhiMesonTree;

class StPhiMesonMaker : public StMaker {
  public:
    StPhiMesonMaker(const char *name, StPicoDstMaker *picoMaker, const string jobId, const int mode, const int beamType, const int flagME);
    virtual ~StPhiMesonMaker();
    
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
    // StMixEpManager  *mMixEpManager;
    StPhiMesonTree  *mPhiMesonTree;

    const int mMode;
    const int mType;
    const int mFlagME;

    string str_mOutPutRecoPhi;
    string str_mOutPutPhiTree;

    TFile *file_mOutPutRecoPhi;
    TFile *file_mOutPutPhiTree;

    std::string str_mFlagME[2] = {"SE","ME"};

    ClassDef(StPhiMesonMaker, 1)
};

#endif
