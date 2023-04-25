#include <fstream>

#include "TFile.h"
#include "TChain.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TString.h"

#include "StMessMgr.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
// #include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
// #include "StRoot/StEventPlaneMaker/StEpdEpManager.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonEvent.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonAnalyzer.h"

ClassImp(StPhiMesonAnalyzer)

//----------------------------------------------------
StPhiMesonAnalyzer::StPhiMesonAnalyzer(const int beamType, const int mode, const int flagME, const int listId, const long startEvt, const long stopEvt) : mType(beamType), mMode(mode), mFlagME(flagME), mListId(listId)
{
  mFlagInPut = 1;
  setStartEvent(startEvt);
  setStopEvent(stopEvt);

  string inputdir = Form("/star/data01/pwg/sunxuhit/%s/SpinAlignment/PhiMesonMaker/Forest/",globCons::str_mBeamType[mType].c_str());
  setInPutDir(inputdir);

  const int startListId = mNumList*mListId + 1; // start list
  const int stopListId  = mNumList*(mListId+1); // stop list

  string inputList = Form("/star/data01/pwg/sunxuhit/%s/SpinAlignment/PhiMesonMaker/List/List_RecoPhi%s_%s_%d_%d.list",globCons::str_mBeamType[mType].c_str(),str_mMixEvt[mFlagME].c_str(),globCons::str_mBeamType[mType].c_str(),startListId,stopListId);
  setInPutList(inputList);

  if(mMode == 0)
  {
    string outputfile = Form("./file_QaPhi%s_%s_%d.root",str_mMixEvt[mFlagME].c_str(),globCons::str_mBeamType[mType].c_str(),mListId);
    setOutPutfile(outputfile);
  }
}

StPhiMesonAnalyzer::~StPhiMesonAnalyzer()
{
}
//----------------------------------------------------
// set Input/Output
void StPhiMesonAnalyzer::setInPutDir(const string inputDir)
{
  str_mInPutDir = inputDir;
  cout << "Input directory was set to: " << str_mInPutDir.c_str() << endl;
}
void StPhiMesonAnalyzer::setInPutList(const string inputList)
{
  str_mInPutList = inputList;
  string InFo_InPutList = Form("InPut %s list was set to: %s",str_mMixEvt[mFlagME].c_str(),str_mInPutList.c_str());
  cout << InFo_InPutList.c_str() << endl;
}
void StPhiMesonAnalyzer::setOutPutFile(const string outputFile)
{
  str_mOutPutFile = outputFile;
  cout << "Output file was set to: " << str_mOutPutFile.c_str() << endl;
}
void StPhiMesonAnalyzer::setStartEvent(const long startEvt)
{
    mStartEvt = startEvt;
    cout << "nStartEvent = " << mStartEvt << endl;
}
void StPhiMesonAnalyzer::setStopEvent(const long stopEvt)
{
    mStopEvt = stopEvt;
    cout << "nStopEvent = " << mStopEvt << endl;
}
void StPhiMesonAnalyzer::initChain()
{ // initialize TChain
  if (!str_mInPutList.empty())   // if input file is ok
  {
    string infoList = Form("Open %s file list ",str_mMixEvt[mFlagME].c_str());
    cout << infoList.c_str() << endl;
    ifstream in(str_mInPutList);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      c_mInPut = new TChain("PhiMesonEvent", "PhiMesonEvent");
      char str[255];       // char array for each file name
      long entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  string addfile;
	  addfile = str;
	  addfile = str_mInPutDir+addfile;
	  c_mInPut->AddFile(addfile.c_str(),-1,"PhiMesonEvent");
	  long file_entries = c_mInPut->GetEntries();
	  cout << "File added to data chain: " << addfile.c_str() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      string infoWarning = Form("WARNING: %s file input is problemtic",str_mMixEvt[mFlagME].c_str());
      cout << infoWarning.c_str() << endl;
      mFlagInPut = 0;
    }
  }

  // Set the input tree
  if (mFlagInPut == 1 && !c_mInPut->GetBranch( "phiSpinAlignmentBranch" ))
  {
    cerr << "ERROR: Could not find branch 'phiSpinAlignmentBranch' in tree!" << endl;
  }

  if(mFlagInPut == 1)
  {
    mPhiEvt = new StPhiMesonEvent();
    c_mInPut->SetBranchAddress("phiSpinAlignmentBranch",&mPhiEvt);

    int nEvts = c_mInPut->GetEntriesFast();
    cout << "Number of events in file(s) = " << nEvts << endl;
    if(mStartEvt > nEvts) mStartEvt = nEvts;
    if(mStopEvt > nEvts) mStopEvt   = nEvts;

    cout << "New nStartEvent = " << mStartEvt << ", new nStopEvent = " << mStopEvt << endl;
  }
}
//----------------------------------------------------
// initial functions
void StPhiMesonAnalyzer::Init()
{
  mZdcEpManager = new StZdcEpManager(mType); // initialize ZDC EP Manager
  mEpdEpManager = new StEpdEpManager(mType); // initialize EPD EP Manager
  mTpcEpManager = new StTpcEpManager(mType); // initialize TPC EP Manager
  mHistManager  = new StPhiMesonHistoManger(mType); // initialize histogram manager
  mAnaCut       = new StAnalysisCut(mType);
  mAnaUtils     = new StAnalysisUtils(mType);
  // mAnaUtils->initRunIndex(); // initialize std::map for run index

  initChain();

  mTpcEpManager->readTpcReCtr(); // TPC
  mTpcEpManager->readTpcShift();

  if(mMode == 0)
  {
    file_mOutPutQA = new TFile(str_mOutPutFile.c_str(),"RECREATE");
    mHistManager->initPhiQA();
  }
  if(mMode == 1)
  {
    file_mOutPutFlow = new TFile(str_mOutPutFile.c_str(),"RECREATE");
  }
  if(mMode == 1)
  {
    file_mOutPutSpin = new TFile(str_mOutPutFile.c_str(),"RECREATE");
  }
}
//-------------------------------------------------------------------
void StPhiMesonAnalyzer::Finish()
{
  if(mMode == 0)
  {
    if(str_mOutPutFile != "")
    {
      file_mOutPutQA->cd();
      mHistManager->writePhiQA();
      file_mOutPutQA->Close();
    }
  }
  if(mMode == 1)
  {
    if(str_mOutPutFlow != "")
    {
      file_mOutPutFlow->cd();
      file_mOutPutFlow->Close();
    }
  }
  if(mMode == 2)
  {
    if(str_mOutPutSpin != "")
    {
      file_mOutPutSpin->cd();
      file_mOutPutSpin->Close();
    }
  }
}
//-------------------------------------------------------------------

// loop phi meson events
void StPhiMesonAnalyzer::Make()
{
  long startEvtUsed;
  long stopEvtUsed;

  startEvtUsed = mStartEvt;
  stopEvtUsed  = mStopEvt;
  c_mInPut->SetBranchAddress("phiSpinAlignmentBranch",&mPhiEvt);
  c_mInPut->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  // Initialise Event Head
  int runId       = -1;
  int runIdx      = -1;
  int evtId       = -1;
  int refMult     = -1;
  int numTofMatch = -1;
  int cent9       = -1;
  int cent16      = -1;
  double refWgt   = -1.0;
  double zDCx     = -1.0;
  double bBCx     = -1.0;
  double vzVpd    = -1.0;
  TVector3 mPrimVtx(-999.9,-999.9,-999.9);

  int flagZdcEp = -1; // ZDC EP
  TVector2 vQ1ZdcShiftEast(0.0,0.0);
  TVector2 vQ1ZdcShiftWest(0.0,0.0);
  TVector2 vQ1ZdcShiftFull(0.0,0.0);

  int flagEpdSideEp = -1; // EPD Side EP
  TVector2 vQ1EpdSideShiftEast(0.0,0.0);
  TVector2 vQ1EpdSideShiftWest(0.0,0.0);
  TVector2 vQ1EpdSideShiftFull(0.0,0.0);

  int flagEpdGrp0Ep = -1; // EPD Grp0 EP
  TVector2 vQ1EpdGrp0ShiftEast(0.0,0.0);
  TVector2 vQ1EpdGrp0ShiftWest(0.0,0.0);
  TVector2 vQ1EpdGrp0ShiftFull(0.0,0.0);
  int flagEpdGrp1Ep = -1; // EPD Grp1 EP
  TVector2 vQ1EpdGrp1ShiftEast(0.0,0.0);
  TVector2 vQ1EpdGrp1ShiftWest(0.0,0.0);
  TVector2 vQ1EpdGrp1ShiftFull(0.0,0.0);

  int flagTpcEp = -1; // TPC EP
  int vQ1TpcReCtrEast(0.0,0.0);
  int vQ1TpcReCtrWest(0.0,0.0);
  int vQ2TpcReCtrEast(0.0,0.0);
  int vQ2TpcReCtrWest(0.0,0.0);
  int vQ3TpcReCtrEast(0.0,0.0);
  int vQ3TpcReCtrWest(0.0,0.0);
  int numTrkReCtrEast = -1;
  int numTrkReCtrWest = -1;

  unsigned short numTrkUsed = -1;

  for(long counter = startEvtUsed; counter < stopEvtUsed; counter++)
  {
    if (!c_mInPut->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    // display event process
    if(counter != 0  &&  counter % 1000 == 0)
    {
      cout << "." << flush;
    }
    if(counter != 0  &&  counter % 10000 == 0)
    {
      if((stopEvtUsed-startEvtUsed) > 0)
      {
	double evtPercent = 100.0*((double)(counter-startEvtUsed))/((double)(stopEvtUsed-startEvtUsed));
	cout << " " << counter-startEvtUsed << " (" << evtPercent << "%) " << "\n" << "==> Processing data (VecMesonSpinAlignment) " << flush;
      }
    }

    // get Event Header
    runId       = mPhiEvt->getRunId(); // event header
    runIdx      = mPhiEvt->getRunIdx();
    evtId       = mPhiEvt->getEvtId();
    refMult     = mPhiEvt->getRefMult();
    numTofMatch = mPhiEvt->getNumTofMatch();
    cent9       = mPhiEvt->getCentrality9();
    cent16      = mPhiEvt->getCentrality16();
    refWgt      = mPhiEvt->getRefWgt();
    zDCx        = mPhiEvt->getZDCx();
    bBCx        = mPhiEvt->getBBCx();
    vzVpd       = mPhiEvt->getVzVpd();
    mPrimVtx    = mPhiEvt->getPrimVtx();

    flagZdcEp       = mPhiEvt->getFlagZdcEp(); // ZDC EP Info
    vQ1ZdcShiftEast = mPhiEvt->getQ1VecZdcEast();
    vQ1ZdcShiftWest = mPhiEvt->getQ1VecZdcWest();
    vQ1ZdcShiftFull = mPhiEvt->getQ1VecZdcFull();

    flagEpdSideEp       = mPhiEvt->getFlagEpdSideEp(); // EPD EP Side Info
    vQ1EpdSideShiftEast = mPhiEvt->getQ1VecEpdSideEast();
    vQ1EpdSideShiftWest = mPhiEvt->getQ1VecEpdSideWest();
    vQ1EpdSideShiftFull = mPhiEvt->getQ1VecEpdSideFull();

    flagEpdGrp0Ep       = mPhiEvt->getFlagEpdGrp0Ep(); // EPD EP Grp0 Info
    vQ1EpdGrp0ShiftEast = mPhiEvt->getQ1VecEpdGrp0East();
    vQ1EpdGrp0ShiftWest = mPhiEvt->getQ1VecEpdGrp0West();
    vQ1EpdGrp0ShiftFull = mPhiEvt->getQ1VecEpdGrp0Full();
    flagEpdGrp1Ep       = mPhiEvt->getFlagEpdGrp1Ep(); // EPD EP Grp1 Info
    vQ1EpdGrp1ShiftEast = mPhiEvt->getQ1VecEpdGrp1East();
    vQ1EpdGrp1ShiftWest = mPhiEvt->getQ1VecEpdGrp1West();
    vQ1EpdGrp1ShiftFull = mPhiEvt->getQ1VecEpdGrp1Full();

    flagTpcEp       = mPhiEvt->getFlagTpcEp(); // TPC EP Info
    vQ1TpcReCtrEast = mPhiEvt->getQ1VecTpcEast();
    vQ1TpcReCtrWest = mPhiEvt->getQ1VecTpcWest();
    vQ2TpcReCtrEast = mPhiEvt->getQ2VecTpcEast();
    vQ2TpcReCtrWest = mPhiEvt->getQ2VecTpcWest();
    vQ3TpcReCtrEast = mPhiEvt->getQ3VecTpcEast();
    vQ3TpcReCtrWest = mPhiEvt->getQ3VecTpcWest();
    numTrkReCtrEast = mPhiEvt->getNumTrkReCtrEast();
    numTrkReCtrWest = mPhiEvt->getNumTrkReCtrWest();
    numTrkUsed      = mPhiEvt->getNumTracks();

    const int vzBin  = mAnaUtils->getVzBin(vz); // 0 for -vz || 1 for +vz

    // Initialise Track 
    TVector3 vTrkMomKp(0.0,0.0,0.0);
    TVector3 vTrkMomKm(0.0,0.0,0.0);
    double mass2Kp = -999.9;
    double mass2Km = -999.9;
    double nSigKp  = -999.9;
    double nSigKm  = -999.9;
    double dcaKp   = -999.9;
    double dcaKm   = -999.9;
    int chargeKp   = -999;
    int chargeKm   = -999;
    int nHitsFitKp = -999;
    int nHitsFitKm = -999;
    int flagKp     = -1;
    int flagKm     = -1;

    if(mMode == 0)
    {
      for(unsigned short iTrk = 0; iTrk < numTrkUsed; ++iTrk) // loop over all tracks of the actual event
      {
	mPhiTrk    = mPhiEvt->getTrack(iTrk); // get Track Information
	vTrkMomKp = mPhiTrk->getTrkMomKp(); // K+
	vTrkMomKm = mPhiTrk->getTrkMomKm(); // K-
	mass2Kp    = mPhiTrk->getMass2Kp(); // K+
	mass2Km    = mPhiTrk->getMass2Km(); // K-
	nSigKp     = mPhiTrk->getNSigKp(); // K+
	nSigKm     = mPhiTrk->getNSigKm(); // K-
	dcaKp      = mPhiTrk->getDcaKp(); // K+
	dcaKm      = mPhiTrk->getDcaKm(); // K-
	chargeKp   = mPhiTrk->getChargeKp(); // K+
	chargeKm   = mPhiTrk->getChargeKm(); // K-
	nHitsFitKp = mPhiTrk->getNHitsFitKp(); // K+
	nHitsFitKm = mPhiTrk->getNHitsFitKm(); // K-
	flagKp     = mPhiTrk->getFlagKp(); // K+
	flagKm     = mPhiTrk->getFlagKm(); // K-

	TLorentzVector lTrkKp, lTrkKm, lTrkPhi;
	lTrkKp.SetXYZM(vTrkMomKp.X(),vTrkMomKp.Y(),vTrkMomKp.Z(),anaUtils::mMassKaon);
	lTrkKm.SetXYZM(vTrkMomKm.X(),vTrkMomKm.Y(),vTrkMomKm.Z(),anaUtils::mMassKaon);
	lTrkPhi = lTrkKp+lTrkKm;

	if( (mass2Kp > 0.16 && mass2Kp < 0.36) && (mass2Km > 0.16 && mass2Km < 0.36) )
	{ // always require ToF Info
	  float ptPhi = lTrkPhi.Perp();
	  float yPhiLab = lTrkPhi.Rapidity();
	  float yPhiCms = mAnaUtils->getRapidityCMS(yPhiLab);
	  float invMassPhi = lTrkPhi.M();
	  mHistManager->fillPhiQA(cent9,ptPhi,yPhiLab,yPhiCms,invMassPhi,refWgt);
	}
      }
    }
    if(mMode == 1)
    {
      if(mAnaCut->isIsobar() && flagEpdSideEp == 1 && flagTpcEp == 1)
      { // Isobar with valid EPD side EP and TPC EP
      }
      if(mAnaCut->isFxt3p85GeV_2018() && flagEpdGrp0Ep == 1 && flagEpdGrp1Ep == 1 && flagTpcEp == 1)
      { // Fxt3p85GeV_2018 with valid EPD Grp0 & Grp1 EP and TPC EP
      }
    }
  }

  cout << "." << flush;
  cout << " " << stopEvtUsed-startEvtUsed << "(" << 100 << "%)";
  cout << endl;
}
