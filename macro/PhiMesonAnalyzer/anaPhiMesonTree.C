#include <TSystem>
#include "TStopwatch.h"

void anaPhiMesonTree(const string inputList = "Utility/FileList/ZrZr200GeV_2018/forestRecoPhiSEtest_ZrZr200GeV_2018.list", const string jobId = "14", const int beamType = 0, const int mode = 0, const int flagME = 0, const long startEvt = 0, const long stopEvt = 100000024)
// void anaPhiMesonTree(const string inputList = "Utility/FileList/ZrZr200GeV_2018/forestRecoPhiMEtest_ZrZr200GeV_2018.list", const string jobId = "14", const int beamType = 0, const int mode = 0, const int flagME = 1, const long startEvt = 0, const long stopEvt = 100000024)
// void anaPhiMesonTree(const string inputList = "Utility/FileList/RuRu200GeV_2018/forestRecoPhiSEtest_RuRu200GeV_2018.list", const string jobId = "14", const int beamType = 1, const int mode = 0, const int flagME = 0, const long startEvt = 0, const long stopEvt = 100000024)
// void anaPhiMesonTree(const string inputList = "Utility/FileList/RuRu200GeV_2018/forestRecoPhiMEtest_RuRu200GeV_2018.list", const string jobId = "14", const int beamType = 1, const int mode = 0, const int flagME = 1, const long startEvt = 0, const long stopEvt = 100000024)
// void anaPhiMesonTree(const string inputList = "Utility/FileList/Fxt3p85GeV_2018/forestRecoPhiSEtest_Fxt3p85GeV_2018.list", const string jobId = "14", const int beamType = 2, const int mode = 0, const int flagME = 0, const long startEvt = 0, const long stopEvt = 100000024)
// void anaPhiMesonTree(const string inputList = "Utility/FileList/Fxt3p85GeV_2018/forestRecoPhiMEtest_Fxt3p85GeV_2018.list", const string jobId = "14", const int beamType = 2, const int mode = 0, const int flagME = 1, const long startEvt = 0, const long stopEvt = 100000024)
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
