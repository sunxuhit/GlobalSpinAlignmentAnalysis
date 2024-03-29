#ifndef StPhiMesonAnalyzer_h
#define StPhiMesonAnalyzer_h

#include "TObject.h"
#include <string>

class TChain;
class TFile;

class StAnalysisUtils;
class StAnalysisCut;
class StTpcEpManager;
class StMixEpManager;
class StPhiMesonEvent;
class StPhiMesonTrack;
class StPhiMesonHistoManger;

class StPhiMesonAnalyzer : public TObject
{
  public:
    StPhiMesonAnalyzer(const string inputList, const string jobId, const int beamType, const int mode, const int flagME, const long startEvt, const long stopEvt); // mode: 0 for QA, 1 for phi invMass | flagME: 0 for Same Event, 1 for Mixed Event | listId: list ID
    virtual ~StPhiMesonAnalyzer();

    void initChain();

    void Init();
    void Make();
    void Finish();

  private:
    // static const int mNumList = 20; // number of files per list

    StAnalysisUtils *mAnaUtils; // Analysis Utilities
    StAnalysisCut   *mAnaCut;   // Analysis Cuts
    StTpcEpManager  *mTpcEpManager;
    StMixEpManager  *mMixEpManager; // for FXT ONLY
    StPhiMesonEvent *mPhiEvt;
    StPhiMesonTrack *mPhiTrk;
    StPhiMesonHistoManger *mHistManager;

    int mFlagInPut; // 0 for problematic input list, 1 for good input list
    const int mType;
    const int mMode; // 0 for QA | 1 for phi flow | 2 for phi spin alignment
    const int mFlagME; // 0 for Same Event, 1 for Mixed Event

    string str_mInPutList;
    string str_mOutPutQA;
    string str_mOutPutFlow;
    string str_mOutPutSpin;

    long mStopEvt;
    long mStartEvt;

    TChain *c_mInPutPhiEvt;
    TFile *file_mOutPutQA;
    TFile *file_mOutPutFlow;
    TFile *file_mOutPutSpin;

    std::string str_mMixEvt[2] = {"SE","ME"};

  ClassDef(StPhiMesonAnalyzer,1)
};
#endif
