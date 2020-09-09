#include <string>
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
  string JobId = "B95F74FE855427A3F4DF102A1FDB0FC5";
  string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mTriggerID[2]; // 0: before cuts | 1: after cuts

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    string HistName = Form("h_mTriggerID_%s",CutsQA[i_cut].c_str());
    h_mTriggerID[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mTriggerID[i_cut]->SetLineColor(i_cut+1);
    h_mTriggerID[i_cut]->GetXaxis()->SetTitle("triggerID");
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
  }
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    c_TriggerId->cd(i_cut+1);
    h_mTriggerID[i_cut]->Draw();
  }
}
