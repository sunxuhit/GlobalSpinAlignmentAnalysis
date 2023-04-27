#include <TSystem>
#include "TStopwatch.h"

void anaPhiMesonTree(const string inputList = "/star/u/sunxuhit/ZrZr200GeV_2018/SpinAlignment/PhiMesonMaker/List/RecoPhiSE_ZrZr200GeV_2018_1_20.list", const string jobId = "14", const int beamType = 0, const int mode = 0, const int flagME = 0, const long startEvt = 0, const long stopEvt = 10024)
{
  // mBeamType[NumBeamType] = {"ZrZr200GeV_2018","RuRu200GeV_2018","Fxt3p85GeV_2018"};
  // mode: 0 for QA, 1 for phi flow, 2 for phi spin alignment
  // flagME: 0 for Same Event, 1 for Mixed Event

  TStopwatch *stopWatch = new TStopwatch();
  stopWatch->Start();

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StAnalysisUtils");
  gSystem->Load("StEventPlaneMaker");
  gSystem->Load("StPhiMesonMaker");
  gSystem->Load("StPhiMesonAnalyzer");

  cout << "All libraries are loaded!!!!" << endl;

  cout << "Start to Read Trees!" << endl;

  StPhiMesonAnalyzer *phiMesonAna = new StPhiMesonAnalyzer(inputList,jobId,beamType,mode,flagME,startEvt,stopEvt);
  phiMesonAna->Init();
  phiMesonAna->Make();
  phiMesonAna->Finish();

  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;

  stopWatch->Stop();
  stopWatch->Print();

  delete phiMesonAna;
}
