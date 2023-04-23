#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void plotTpcEpQA(int beamType = 0, string jobId = "CEBB7148080C11DF693CC3046509174F")
{
  TCanvas *c_TpcEpQA = new TCanvas("c_TpcEpQA","c_TpcEpQA",10,10,900,900);
  c_TpcEpQA->Divide(3,3);
  for(int iCent = 0; iCent < 9; ++iCent)
  {
    c_TpcEpQA->cd(iCent+1);
    c_TpcEpQA->cd(iCent+1)->SetLeftMargin(0.15);
    c_TpcEpQA->cd(iCent+1)->SetBottomMargin(0.15);
    c_TpcEpQA->cd(iCent+1)->SetTicks(1,1);
    c_TpcEpQA->cd(iCent+1)->SetGrid(0,0);
    c_TpcEpQA->cd(iCent+1)->SetLogy(1);
  }
  std::string figName = Form("../../figures/EventPlaneMaker/%s/TpcEpQA_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TpcEpQA->Print(figName.c_str());
  figName = Form("../../figures/EventPlaneMaker/%s/TpcEpQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());

  string inputList = Form("../../data/EventPlaneMaker/%s/%s.list",globCons::str_mBeamType[beamType].c_str(),jobId.c_str());
  if (!inputList.empty())   // if input file is ok
  {
    cout << "Open input probability file list" << endl;
    ifstream in(inputList.c_str());  // input stream
    if(in)
    {
      cout << "input file probability list is ok" << endl;
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
        in.getline(str,255);  // take the lines of the file list
        if(str[0] != 0)
        {
          string addfile = str;
	  string inputFile = Form("../../data/EventPlaneMaker/%s/%s",globCons::str_mBeamType[beamType].c_str(),addfile.c_str());
	  TFile *file_InPut = TFile::Open(inputFile.c_str());
	  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
	  cout << "inputFile sets to: " << inputFile.c_str() << endl;

	  TH2F *h_mTpcEp2[9];
	  for(int iCent = 0; iCent < 9; ++iCent)
	  {
	    string histName = Form("h_mTpcEp2ShiftWestCent%d",iCent);
	    h_mTpcEp2[iCent] = (TH2F*)file_InPut->Get(histName.c_str())->Clone();
	    h_mTpcEp2[iCent]->GetYaxis()->SetTitle(addfile.c_str());
	  }

	  for(int iCent = 0; iCent < 9; ++iCent)
	  {
	    c_TpcEpQA->cd(iCent+1);
	    h_mTpcEp2[iCent]->ProjectionY()->DrawCopy();
	  }
	  c_TpcEpQA->Update();
	  c_TpcEpQA->Print(figName.c_str());

	  file_InPut->Close();
        }
      }
    }
    else
    {
      cout << "WARNING: input probability file input is problemtic" << endl;
      return false;
    }
  }

  figName = Form("../../figures/EventPlaneMaker/%s/TpcEpQA_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TpcEpQA->Print(figName.c_str());
}
