#include <algorithm>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TVector3.h>
#include <TMath.h>

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
#include "StRoot/StEventPlaneMaker/StEpdEpManager.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"
// #include "StRoot/StPhiMesonMaker/StPhiMesonTree.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonMaker.h"

ClassImp(StPhiMesonMaker)

//-----------------------------------------------------------------------------
StPhiMesonMaker::StPhiMesonMaker(const char* name, StPicoDstMaker *picoMaker, const string jobId, const int mode, const int beamType, const int flagME) : StMaker(name), mMode(mode), mType(beamType), mFlagME(flagME)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = NULL;
  mRefMultCorr = NULL;

  if(mMode == 0)
  {
    // mOutPut_QA = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/QA/file_%s_QA_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),jobId.c_str());
    mOutPut_QA = Form("./file_%s_QA_%s.root",vmsa::mBeamEnergy[energy].c_str(),jobId.c_str());
  }
}

//----------------------------------------------------------------------------- 
StPhiMesonMaker::~StPhiMesonMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
int StPhiMesonMaker::Init() 
{
  mZdcEpManager = new StZdcEpManager(mType); // initialize ZDC EP Manager
  mEpdEpManager = new StEpdEpManager(mType); // initialize EPD EP Manager
  mTpcEpManager = new StTpcEpManager(mType); // initialize TPC EP Manager
  mAnaCut       = new StAnalysisCut(mType);
  mAnaUtils     = new StAnalysisUtils(mType);
  mAnaUtils->initRunIndex(); // initialize std::map for run index

  if(!mRefMultCorr)
  {
    if( mAnaCut->isIsobar() ) mRefMultCorr = CentralityMaker::instance()->getRefMultCorr_Isobar();
  }

  if(mMode == 0)
  { // QA
    mFile_QA= new TFile(mOutPut_QA.c_str(),"RECREATE");
    mFile_QA->cd();
    mVecMesonHistoManager->initEventQA();
    mVecMesonHistoManager->initTrackQA();
    mVecMesonProManager->initRunQA();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StPhiMesonMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_QA != "")
    {
      mFile_QA->cd();
      mVecMesonHistoManager->writeEventQA();
      mVecMesonHistoManager->writeTrackQA();
      mVecMesonProManager->writeRunQA();
      mFile_QA->Close();
    }
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StPhiMesonMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
int StPhiMesonMaker::Make() 
{

  return kStOK;
}

