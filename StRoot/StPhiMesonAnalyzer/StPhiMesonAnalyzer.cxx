#include <fstream>

#include "TChain.h"
#include "TFile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TString.h"

#include "StMessMgr.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"
#include "StRoot/StAnalysisUtils/StAnalysisCut.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"
#include "StRoot/StEventPlaneMaker/StMixEpManager.h"
#include "StRoot/StPhiMesonMaker/StPhiMesonEvent.h"
#include "StRoot/StPhiMesonAnalyzer/StPhiMesonHistoManger.h"
#include "StRoot/StPhiMesonAnalyzer/StPhiMesonAnalyzer.h"

ClassImp(StPhiMesonAnalyzer)

//----------------------------------------------------
StPhiMesonAnalyzer::StPhiMesonAnalyzer(const string inputList, const string jobId, const int beamType, const int mode, const int flagME, const long startEvt, const long stopEvt) : mType(beamType), mMode(mode), mFlagME(flagME)
{
  mFlagInPut = 1;
  mStartEvt = startEvt;
  mStopEvt = stopEvt;
  cout << "nStartEvent = " << mStartEvt << ", nStopEvent = " << mStopEvt << endl;

  str_mInPutList = inputList;
  string infoInPutList = Form("InPut %s list was set to: %s",str_mMixEvt[mFlagME].c_str(),str_mInPutList.c_str());
  cout << infoInPutList.c_str() << endl;

  if(mMode == 0)
  {
    str_mOutPutQA = Form("./file_QaPhi%s_%s_%s.root",str_mMixEvt[mFlagME].c_str(),globCons::str_mBeamType[mType].c_str(),jobId.c_str());
    cout << "Output file was set to: " << str_mOutPutQA.c_str() << endl;
  }
  if(mMode == 1)
  {
    str_mOutPutFlow = Form("./file_FlowPhi%s_%s_%s.root",str_mMixEvt[mFlagME].c_str(),globCons::str_mBeamType[mType].c_str(),jobId.c_str());
    cout << "Output file was set to: " << str_mOutPutFlow.c_str() << endl;
  }
  if(mMode == 2)
  {
    str_mOutPutSpin = Form("./file_SpinPhi%s_%s_%s.root",str_mMixEvt[mFlagME].c_str(),globCons::str_mBeamType[mType].c_str(),jobId.c_str());
    cout << "Output file was set to: " << str_mOutPutSpin.c_str() << endl;
  }
}

StPhiMesonAnalyzer::~StPhiMesonAnalyzer()
{
}
//----------------------------------------------------
void StPhiMesonAnalyzer::initChain()
{ // initialize TChain
  if (!str_mInPutList.empty())   // if input file is ok
  {
    string infoList = Form("Open %s file list ",str_mMixEvt[mFlagME].c_str());
    cout << infoList.c_str() << endl;
    ifstream file_InPutList(str_mInPutList);  // input stream
    if(file_InPutList)
    {
      cout << "input file list is ok" << endl;
      c_mInPutPhiEvt = new TChain("PhiMesonEvent", "PhiMesonEvent");
      char str[255];       // char array for each file name
      long evtSave = 0;
      while(file_InPutList)
      {
	file_InPutList.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  string addfile = str;
	  c_mInPutPhiEvt->AddFile(addfile.c_str(),-1,"PhiMesonEvent");
	  long evtInPut = c_mInPutPhiEvt->GetEntries();
	  cout << "File added to data chain: " << addfile.c_str() << " with " << (evtInPut-evtSave) << " entries" << endl;
	  evtSave = evtInPut;
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
  if (mFlagInPut == 1 && !c_mInPutPhiEvt->GetBranch( "phiSpinAlignmentBranch" ))
  {
    cerr << "ERROR: Could not find branch 'phiSpinAlignmentBranch' in tree!" << endl;
  }

  if(mFlagInPut == 1)
  {
    mPhiEvt = new StPhiMesonEvent();
    c_mInPutPhiEvt->SetBranchAddress("phiSpinAlignmentBranch",&mPhiEvt);

    int nEvts = c_mInPutPhiEvt->GetEntriesFast();
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
  mTpcEpManager = new StTpcEpManager(mType); // initialize TPC EP Manager
  mMixEpManager = new StMixEpManager(mType); // initialize Mix EP Manager
  mHistManager  = new StPhiMesonHistoManger(mType,mFlagME); // initialize histogram manager
  mAnaCut       = new StAnalysisCut(mType);
  mAnaUtils     = new StAnalysisUtils(mType);

  initChain();

  if(mMode == 0)
  {
    file_mOutPutQA = new TFile(str_mOutPutQA.c_str(),"RECREATE");
    mHistManager->initPhiQA();
  }
  if(mMode == 1)
  {
    file_mOutPutFlow = new TFile(str_mOutPutFlow.c_str(),"RECREATE");
    if(mAnaCut->isIsobar()) 
    {
      mTpcEpManager->readTpcReCtr(); // TPC
      mTpcEpManager->readTpcShift();
      mTpcEpManager->readTpcResolution();
      mHistManager->initIsoPhiFlow();
    }
    if(mAnaCut->isFxt3p85GeV_2018())
    {
      mMixEpManager->readMixEpRes(); // Mix
      mHistManager->initFxtPhiFlow();
    }
  }
  if(mMode == 2)
  {
    file_mOutPutSpin = new TFile(str_mOutPutSpin.c_str(),"RECREATE");
  }
}
//-------------------------------------------------------------------
void StPhiMesonAnalyzer::Finish()
{
  if(mMode == 0)
  {
    if(str_mOutPutQA != "")
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
      if(mAnaCut->isIsobar()) mHistManager->writeIsoPhiFlow();
      if(mAnaCut->isFxt3p85GeV_2018()) mHistManager->writeFxtPhiFlow();
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
  long startEvtUsed = mStartEvt;
  long stopEvtUsed  = mStopEvt;

  c_mInPutPhiEvt->SetBranchAddress("phiSpinAlignmentBranch",&mPhiEvt);
  c_mInPutPhiEvt->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

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
  TVector2 vQ1TpcReCtrEast(0.0,0.0);
  TVector2 vQ1TpcReCtrWest(0.0,0.0);
  TVector2 vQ2TpcReCtrEast(0.0,0.0);
  TVector2 vQ2TpcReCtrWest(0.0,0.0);
  TVector2 vQ3TpcReCtrEast(0.0,0.0);
  TVector2 vQ3TpcReCtrWest(0.0,0.0);
  int numTrkReCtrEast = -1;
  int numTrkReCtrWest = -1;

  unsigned short numTrkUsed = -1;

  for(long iEvt = startEvtUsed; iEvt < stopEvtUsed; iEvt++)
  {
    if( !c_mInPutPhiEvt->GetEntry(iEvt) ) // take the event -> information is stored in event
      break;  // end of data chunk

    // display event process
    if(iEvt != 0  &&  iEvt % 1000 == 0) { cout << "." << flush; }
    if(iEvt != 0  &&  iEvt % 10000 == 0)
    {
      if((stopEvtUsed-startEvtUsed) > 0)
      {
	double evtPercent = 100.0*((double)(iEvt-startEvtUsed))/((double)(stopEvtUsed-startEvtUsed));
	cout << " " << iEvt-startEvtUsed << " (" << evtPercent << "%) " << "\n" << "==> Processing data (VecMesonSpinAlignment) " << flush;
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

    const int vzBin  = mAnaUtils->getVzBin(mPrimVtx.Z()); // 0 for -vz || 1 for +vz
    mTpcEpManager->initTpcEpManager(cent9,runIdx,vzBin); // initialize TPC EP Manager
    mMixEpManager->initMixEpManager(cent9,runIdx,vzBin); // initialize Mix EP Manager

    // Initialise Track 
    TVector3 vTrkMomKp(0.0,0.0,0.0);
    TVector3 vTrkMomKm(0.0,0.0,0.0);
    double mass2Kp = -999.9;
    double mass2Km = -999.9;
    double betaKp  = -999.9;
    double betaKm  = -999.9;
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

    for(unsigned short iTrk = 0; iTrk < numTrkUsed; ++iTrk) // loop over all tracks of the actual event
    {
      mPhiTrk    = mPhiEvt->getTrack(iTrk); // get Track Information
      vTrkMomKp  = mPhiTrk->getTrkMomKp(); // K+
      vTrkMomKm  = mPhiTrk->getTrkMomKm(); // K-
      mass2Kp    = mPhiTrk->getMass2Kp(); // K+
      mass2Km    = mPhiTrk->getMass2Km(); // K-
      betaKp     = mPhiTrk->getBetaKp(); // K+
      betaKm     = mPhiTrk->getBetaKm(); // K-
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

      TLorentzVector lTrkKp, lTrkKm;
      lTrkKp.SetXYZM(vTrkMomKp.X(),vTrkMomKp.Y(),vTrkMomKp.Z(),anaUtils::mMassKaon);
      lTrkKm.SetXYZM(vTrkMomKm.X(),vTrkMomKm.Y(),vTrkMomKm.Z(),anaUtils::mMassKaon);
      TLorentzVector lTrkPhi = lTrkKp+lTrkKm;

      double ptPhi      = lTrkPhi.Perp();
      double yPhiLab    = lTrkPhi.Rapidity();
      double yPhiCms    = mAnaUtils->getRapidityCMS(yPhiLab);
      double phiPhi     = lTrkPhi.Phi();
      double invMassPhi = lTrkPhi.M();

      if(mMode == 0)
      { // phi invMass QA
	if( mAnaCut->passTrkTofKaonBeta(mPhiTrk) )
	{ // always require ToF Info
	  mHistManager->fillPhiQA(cent9,ptPhi,yPhiLab,yPhiCms,invMassPhi,refWgt);
	}
      }
      if(mMode == 1)
      { // phi flow
	if(mAnaCut->isIsobar() && flagEpdSideEp == 1 && flagTpcEp == 1)
	{ // Isobar with valid EPD side EP and TPC EP
	  if( mAnaCut->passTrkTofKaonBeta(mPhiTrk) )
	  { // always require ToF Info
	    if(mAnaCut->passTrkPhiFlowEast(yPhiCms))
	    { // get EP from West
	      TVector2 Q2VecWest = vQ2TpcReCtrWest;
	      TVector2 Q3VecWest = vQ3TpcReCtrWest;
	      if(flagKp == 0 && mAnaCut->passTrkTpcEpWest(vTrkMomKp,dcaKp))
	      { // Kp used in West EP
		double wgtTpc = mTpcEpManager->getWeight(vTrkMomKp);

		TVector2 q2VecKp = mTpcEpManager->calq2Vector(vTrkMomKp);
		TVector2 q2VecCtrKp = mTpcEpManager->getq2VecCtrWest();
		Q2VecWest = Q2VecWest - wgtTpc*(q2VecKp-q2VecCtrKp);

		TVector2 q3VecKp = mTpcEpManager->calq3Vector(vTrkMomKp);
		TVector2 q3VecCtrKp = mTpcEpManager->getq3VecCtrWest();
		Q3VecWest = Q3VecWest - wgtTpc*(q3VecKp-q3VecCtrKp);
	      }
	      if(flagKm == 0 && mAnaCut->passTrkTpcEpWest(vTrkMomKm,dcaKm))
	      { // Km used in West EP
		double wgtTpc = mTpcEpManager->getWeight(vTrkMomKm);

		TVector2 q2VecKm = mTpcEpManager->calq2Vector(vTrkMomKm);
		TVector2 q2VecCtrKm = mTpcEpManager->getq2VecCtrWest();
		Q2VecWest = Q2VecWest - wgtTpc*(q2VecKm-q2VecCtrKm);

		TVector2 q3VecKm = mTpcEpManager->calq3Vector(vTrkMomKm);
		TVector2 q3VecCtrKm = mTpcEpManager->getq3VecCtrWest();
		Q3VecWest = Q3VecWest - wgtTpc*(q3VecKm-q3VecCtrKm);
	      }
	      double res2Sub = mTpcEpManager->getTpcSubEp2ResVal(cent9);
	      double Psi2West = mTpcEpManager->getPsi2ShiftWest(Q2VecWest);
	      double res3Sub = mTpcEpManager->getTpcSubEp3ResVal(cent9);
	      double Psi3West = mTpcEpManager->getPsi3ShiftWest(Q3VecWest);

	      mHistManager->fillIsoPhiV2(cent9, ptPhi, yPhiCms, phiPhi, Psi2West, invMassPhi, res2Sub, refWgt);
	      mHistManager->fillIsoPhiV3(cent9, ptPhi, yPhiCms, phiPhi, Psi3West, invMassPhi, res3Sub, refWgt);
	      mHistManager->fillIsoPhiYields(cent9, ptPhi, yPhiCms, invMassPhi, refWgt);
	    }
	    if(mAnaCut->passTrkPhiFlowWest(yPhiCms))
	    { // get EP from East
	      TVector2 Q2VecEast = vQ2TpcReCtrEast;
	      TVector2 Q3VecEast = vQ3TpcReCtrEast;
	      if(flagKp == 0 && mAnaCut->passTrkTpcEpEast(vTrkMomKp,dcaKp))
	      { // Kp used in East EP
		double wgtTpc = mTpcEpManager->getWeight(vTrkMomKp);

		TVector2 q2VecKp = mTpcEpManager->calq2Vector(vTrkMomKp);
		TVector2 q2VecCtrKp = mTpcEpManager->getq2VecCtrEast();
		Q2VecEast = Q2VecEast - wgtTpc*(q2VecKp-q2VecCtrKp);

		TVector2 q3VecKp = mTpcEpManager->calq3Vector(vTrkMomKp);
		TVector2 q3VecCtrKp = mTpcEpManager->getq3VecCtrEast();
		Q3VecEast = Q3VecEast - wgtTpc*(q3VecKp-q3VecCtrKp);
	      }
	      if(flagKm == 0 && mAnaCut->passTrkTpcEpEast(vTrkMomKm,dcaKm))
	      { // Km used in East EP
		double wgtTpc = mTpcEpManager->getWeight(vTrkMomKm);

		TVector2 q2VecKm = mTpcEpManager->calq2Vector(vTrkMomKm);
		TVector2 q2VecCtrKm = mTpcEpManager->getq2VecCtrEast();
		Q2VecEast = Q2VecEast - wgtTpc*(q2VecKm-q2VecCtrKm);

		TVector2 q3VecKm = mTpcEpManager->calq3Vector(vTrkMomKm);
		TVector2 q3VecCtrKm = mTpcEpManager->getq3VecCtrEast();
		Q3VecEast = Q3VecEast - wgtTpc*(q3VecKm-q3VecCtrKm);
	      }
	      double res2Sub = mTpcEpManager->getTpcSubEp2ResVal(cent9);
	      double Psi2East = mTpcEpManager->getPsi2ShiftEast(Q2VecEast);
	      double res3Sub = mTpcEpManager->getTpcSubEp3ResVal(cent9);
	      double Psi3East = mTpcEpManager->getPsi3ShiftEast(Q3VecEast);

	      mHistManager->fillIsoPhiV2(cent9, ptPhi, yPhiCms, phiPhi, Psi2East, invMassPhi, res2Sub, refWgt);
	      mHistManager->fillIsoPhiV3(cent9, ptPhi, yPhiCms, phiPhi, Psi3East, invMassPhi, res3Sub, refWgt);
	      mHistManager->fillIsoPhiYields(cent9, ptPhi, yPhiCms, invMassPhi, refWgt);
	    }
	  }
	}
	if(mAnaCut->isFxt3p85GeV_2018() && flagEpdGrp0Ep == 1 && flagEpdGrp1Ep == 1 && flagTpcEp == 1)
	{ // Fxt3p85GeV_2018 with valid EPD Grp0 & Grp1 EP and TPC EP
	  if( mAnaCut->passTrkTpcKaonFull(mPhiTrk) && mAnaCut->passTrkTofKaonBeta(mPhiTrk) )
	  { // always require ToF Info
	    double res12Sub = mMixEpManager->getMixSubEp1Res2Val(cent9,0);
	    double Psi1Grp0 = TMath::ATan2(vQ1EpdGrp0ShiftEast.Y(),vQ1EpdGrp0ShiftEast.X()); // Psi1 from EPD East Grp0

	    mHistManager->fillFxtPhiV2(cent9, ptPhi, yPhiCms, phiPhi, Psi1Grp0, invMassPhi, res12Sub, refWgt);
	    mHistManager->fillFxtPhiYields(cent9, ptPhi, yPhiCms, invMassPhi, refWgt);
	  }
	}
      }
    }
    mTpcEpManager->clearTpcEpManager();
    mMixEpManager->clearMixEpManager();
  }

  cout << "." << flush;
  cout << " " << stopEvtUsed-startEvtUsed << "(" << 100 << "%)";
  cout << endl;
}
