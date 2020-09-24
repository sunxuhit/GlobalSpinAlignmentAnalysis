#include <string>
#include <vector>
#include "TStyle.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "../StRoot/StRunQAMaker/StRunQACons.h"

using namespace std;

static const string CutsQA[2] = {"Before","After"};

void plotQA_TriggerId(int energy = 0)
{
  gStyle->SetOptStat(111111);
  gStyle->SetStatX(0.95); gStyle->SetStatY(0.90);
  gStyle->SetStatW(0.35); gStyle->SetStatH(0.20);

  string JobId = "low";
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/test/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());

  vector<string> vecTriggerID;
  vecTriggerID.clear();
  if(energy == 0) 
  {
    vecTriggerID.push_back("450005");
    vecTriggerID.push_back("450015");
    vecTriggerID.push_back("450025");
    vecTriggerID.push_back("450050");
    vecTriggerID.push_back("450060");
  }
  if(energy == 1) 
  {
    vecTriggerID.push_back("580001");
    vecTriggerID.push_back("580021");
  }
  if(energy == 2) 
  {
    vecTriggerID.push_back("610001");
    vecTriggerID.push_back("610011");
    vecTriggerID.push_back("610021");
    vecTriggerID.push_back("610031");
    vecTriggerID.push_back("610041");
    vecTriggerID.push_back("610051");
  }

  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mTriggerID[2]; // 0: before cuts | 1: after cuts

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    string HistName = Form("h_mTriggerID_%s",CutsQA[i_cut].c_str());
    h_mTriggerID[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mTriggerID[i_cut]->SetLineColor(i_cut+1);
    h_mTriggerID[i_cut]->GetXaxis()->SetTitle("triggerID");
    for(int i_trigger = 0; i_trigger < vecTriggerID.size(); ++i_trigger)
    {
      h_mTriggerID[i_cut]->GetXaxis()->SetBinLabel(i_trigger+1,vecTriggerID[i_trigger].c_str());
    }
  }

  TCanvas *c_TriggerId = new TCanvas("c_TriggerId","c_TriggerId",10,10,1600,800);
  c_TriggerId->Divide(2,1);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_TriggerId->cd(i_pad+1);
    c_TriggerId->cd(i_pad+1)->SetLeftMargin(0.1);
    c_TriggerId->cd(i_pad+1)->SetRightMargin(0.1);
    c_TriggerId->cd(i_pad+1)->SetBottomMargin(0.1);
    c_TriggerId->cd(i_pad+1)->SetGrid(0,0);
    c_TriggerId->cd(i_pad+1)->SetTicks(1,1);
    c_TriggerId->cd(i_pad+1)->SetLogy(1);
  }
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    c_TriggerId->cd(i_cut+1);
    h_mTriggerID[i_cut]->Draw();
  }

    string FigName = Form("c_TriggerId_%s_%s.png",runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
    c_TriggerId->SaveAs(FigName.c_str());
}
