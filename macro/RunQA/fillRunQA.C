#include <TSystem>
#include "TStopwatch.h"

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;

StChain *chain;

void fillRunQA(const char *inputFile="Utility/FileList/ZrZr200GeV_2018/pico_xrootd_local.list", const string jobId = "14", const int beamType = 0)
// void fillRunQA(const char *inputFile="Utility/FileList/RuRu200GeV_2018/pico_xrootd_local.list", const string jobId = "14", const int beamType = 1)
// void fillRunQA(const char *inputFile="Utility/FileList/Fxt3p85GeV_2018/pico_xrootd_local.list", const string jobId = "14", const int beamType = 2)
{
  // mBeamType[NumBeamType] = {"ZrZr200GeV_2018","RuRu200GeV_2018","Fxt3p85GeV_2018"};

  TStopwatch *stopWatch = new TStopwatch();
  stopWatch->Start();

  /*
  string SL_version = "pro";
  if(beamType == 0) SL_version = "SL20c"; // ZrZr200GeV_2018
  if(beamType == 1) SL_version = "SL20c"; // RuRu200GeV_2018
  if(beamType == 2) SL_version = "SL20d"; // Fxt3p85GeV_2018
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) 
  {
    cout<<"Environment Star Library does not match the requested library in RunQA.C. Exiting..."<<endl;
    exit(1);
  }
  */

  Int_t nEvents = 10000000000;
  // Int_t nEvents = 100000;

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StEpdUtil");
  gSystem->Load("StPileupUtil");
  gSystem->Load("StAnalysisUtils");
  gSystem->Load("StRunQAMaker");

  chain = new StChain();
  StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst");

  StRunQAMaker *RunQAMaker = new StRunQAMaker("RunQA",picoMaker,jobId,beamType);

  if( chain->Init()==kStErr ){ 
    cout<<"chain->Init();"<<endl;
    return;
  }

  int total = picoMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;

  for (Int_t i=0; i<nEvents; i++)
  {
    if(i != 0 && i%50 == 0)
    {
      cout << "." << flush;
    }
    if(i != 0 && i%500 == 0)
    {
      Float_t event_percent = 100.0*i/nEvents;
      cout << " " << i << "(" << event_percent << "%)" << "\n" << " ==> Processing data " << flush;
    }

    chain->Clear();
    int iret = chain->Make(i);

    if (iret)
    { 
      cout << "Bad return code!" << iret << endl;
      break;
    }

    total++;
  }

  cout << "." << flush;
  cout << " " << nEvents << "(" << 100 << "%)";
  cout << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;

  stopWatch->Stop();
  stopWatch->Print();

  delete RunQAMaker;
  delete picoMaker;
  delete chain;
}
