#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void plotZdcGainCorr(int beamType = 0, string jobId = "test")
{
  TCanvas *c_ZdcGainCorr = new TCanvas("c_ZdcGainCorr","c_ZdcGainCorr",10,10,1600,800);
  c_ZdcGainCorr->Divide(4,2);
  for(int iSlat = 0; iSlat < 8; ++iSlat)
  {
    c_ZdcGainCorr->cd(iSlat+1);
    c_ZdcGainCorr->cd(iSlat+1)->SetLeftMargin(0.15);
    c_ZdcGainCorr->cd(iSlat+1)->SetBottomMargin(0.15);
    c_ZdcGainCorr->cd(iSlat+1)->SetTicks(1,1);
    c_ZdcGainCorr->cd(iSlat+1)->SetGrid(0,0);
    c_ZdcGainCorr->cd(iSlat+1)->SetLogy(1);
  }
  std::string figName = Form("../../figures/EventPlaneMaker/%s/ZdcGainQA_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_ZdcGainCorr->Print(figName.c_str());
  figName = Form("../../figures/EventPlaneMaker/%s/ZdcGainQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());

  string inputList = Form("../../data/EventPlaneMaker/%s/%s.list",globCons::str_mBeamType[beamType].c_str(),jobId.c_str());
  // string inputList = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/EventPlaneMaker/ZrZr200GeV_2018/%s.list",jobId.c_str());
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
	  // string inputFile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/EventPlaneMaker/%s/%s",globCons::str_mBeamType[beamType].c_str(),addfile.c_str());
	  TFile *file_InPut = TFile::Open(inputFile.c_str());
	  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
	  cout << "inputFile sets to: " << inputFile.c_str() << endl;
	  std::string str_mEastWest[2] = {"East","West"};
	  std::string str_mVertHori[2] = {"Vert","Hori"};

	  TH2F *h_mGainInPut[2][2][8];
	  double meanGain[2][2][8], rmsGain[2][2][8];
	  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
	  {
	    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
	    {
	      for(int iSlat = 0; iSlat < 8; ++iSlat)
	      {
		string histName = Form("h_mZdcGainCorr%s%sSlat%d",str_mEastWest[iEastWest].c_str(),str_mVertHori[iVertHori].c_str(),iSlat);
		h_mGainInPut[iEastWest][iVertHori][iSlat] = (TH2F*)file_InPut->Get(histName.c_str())->Clone();
		h_mGainInPut[iEastWest][iVertHori][iSlat]->GetYaxis()->SetTitle(addfile.c_str());
	      }
	    }
	  }

	  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
	  {
	    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
	    {
	      for(int iSlat = 0; iSlat < 8; ++iSlat)
	      {
		c_ZdcGainCorr->cd(iSlat+1);
		h_mGainInPut[iEastWest][iVertHori][iSlat]->ProjectionY()->DrawCopy("colz");
	      }
	      c_ZdcGainCorr->Update();
	      c_ZdcGainCorr->Print(figName.c_str());
	    }
	  }

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

  figName = Form("../../figures/EventPlaneMaker/%s/ZdcGainQA_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_ZdcGainCorr->Print(figName.c_str());
}
