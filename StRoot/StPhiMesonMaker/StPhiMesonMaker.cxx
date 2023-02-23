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
  { // Test phi Meson Reconstruction
    str_mOutPutRecoPhi = Form("./file_RecoPhi_%s_%s.root",globCons::str_mBeamType[beamType].c_str(),jobId.c_str());
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
  { // Test phi Meson Reconstruction
    file_mOutPutRecoPhi = new TFile(str_mOutPutRecoPhi.c_str(),"RECREATE");

    mZdcEpManager->readZdcGain(); // ZDC
    mZdcEpManager->readZdcReCtr();
    mZdcEpManager->readZdcShift();
    mZdcEpManager->readZdcShiftFull();
    mZdcEpManager->readZdcResolution();

    mEpdEpManager->readEpdPhiWgt(); // EPD
    mEpdEpManager->readEpdReCtr();
    mEpdEpManager->readEpdShift();
    mEpdEpManager->readEpdShiftFull();
    mEpdEpManager->readEpdResolution();

    mTpcEpManager->readTpcReCtr(); // TPC
    mTpcEpManager->readTpcShift();
    mTpcEpManager->readTpcResolution();
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
int StPhiMesonMaker::Finish() 
{
  if(mMode == 0)
  {
    if(str_mOutPutRecoPhi != "")
    {
      file_mOutPutRecoPhi->cd();

      file_mOutPutRecoPhi->Close();
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

