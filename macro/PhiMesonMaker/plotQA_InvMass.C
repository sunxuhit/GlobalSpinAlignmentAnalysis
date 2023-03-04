#include <string>
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "../Utility/StSpinAlignmentCons.h"

void plotQA_InvMass(int energy = 6, int flag_ME = 1)
{
  string inputdir = Form("/global/homes/x/xusun/AuAu%s/SpinAlignment/Phi/Forest/",vmsa::mBeamEnergy[energy].c_str());
  string inputlist = Form("/global/homes/x/xusun/AuAu%s/SpinAlignment/Phi/List/Phi_%s_tree.list",vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[flag_ME].Data());

  TCanvas *c_InvMass = new TCanvas("c_InvMass","c_InvMass",10,10,1600,1600);
  c_InvMass->Divide(4,4);
  for(int i_pad = 0; i_pad < 16; ++i_pad)
  {
    c_InvMass->cd(i_pad+1);
    c_InvMass->cd(i_pad+1)->SetLeftMargin(0.1);
    c_InvMass->cd(i_pad+1)->SetBottomMargin(0.1);
    c_InvMass->cd(i_pad+1)->SetGrid(0,0);
    c_InvMass->cd(i_pad+1)->SetTicks(1,1);
  }
  string output_start = Form("../figures/InvMass_AuAu%s_%s_QA.pdf[",vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[flag_ME].Data());
  c_InvMass->Print(output_start.c_str());
  int i_pad = 0;
  int i_page = 1;

  string outputname = Form("../figures/InvMass_AuAu%s_%s_QA.pdf",vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[flag_ME].Data());
  if (!inputlist.empty())   // if input file is ok
  {
    cout << "Open input file list: " << inputlist.c_str() << endl;
    ifstream in(inputlist.c_str());  // input stream
    if(in)
    {
      cout << "input database file list is ok" << endl;
      char str[255];       // char array for each file name
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  string inputfile;
	  inputfile = str;
	  inputfile = inputdir+inputfile;
	  cout << "open file: " << inputfile.c_str() << endl;
	  TFile *File_InPut = TFile::Open(inputfile.c_str());
	  TH2F *h_mMass2_pt = (TH2F*)File_InPut->Get("Mass2_pt")->Clone();
	  TH1F *h_mInvMass = (TH1F*)h_mMass2_pt->ProjectionY()->Clone("h_mInvMass");

	  c_InvMass->cd(i_pad+1);
	  h_mInvMass->DrawCopy("hE");

	  i_pad++;
	  int NumOfTracks = h_mInvMass->GetEntries();
	  cout << "In page " << i_page << " pad " << i_pad << " with " << NumOfTracks << " tracks!" << endl;
	  cout << endl;
	  h_mMass2_pt->Reset();
	  h_mInvMass->Reset();
	  File_InPut->Close();
	  if(i_pad == 16)
	  {
	    i_pad = 0;
	    i_page++;
	    c_InvMass->Update();
	    c_InvMass->Print(outputname.c_str());
	  }
	}
      }
      cout << "Print last page" << endl;
      c_InvMass->Update();
      c_InvMass->Print(outputname.c_str());
    }
    else
    {
      cout << "WARNING: input database file input is problemtic" << endl;
    }
  }

  string output_stop = Form("../figures/InvMass_AuAu%s_%s_QA.pdf]",vmsa::mBeamEnergy[energy].c_str(),vmsa::MixEvent[flag_ME].Data());
  c_InvMass->Print(output_stop.c_str());
}
