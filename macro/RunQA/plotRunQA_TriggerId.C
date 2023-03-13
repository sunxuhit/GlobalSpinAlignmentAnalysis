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
#include "../../Utility/include/StSpinAlignmentCons.h"

using namespace std;

static const string CutStatus[2] = {"Bf","Af"};

void plotRunQA_TriggerId(int beamType = 0)
{
  gStyle->SetOptStat(111111);
  gStyle->SetStatX(0.95); gStyle->SetStatY(0.90);
  gStyle->SetStatW(0.35); gStyle->SetStatH(0.20);

  string inputfile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/%s/RunQA/file_%s_RunQA.root",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());

  vector<string> vecTriggerID;
  vecTriggerID.clear();
  if(beamType == 0 || beamType == 1) 
  {
    vecTriggerID.push_back("600001");
    vecTriggerID.push_back("600011");
    vecTriggerID.push_back("600021");
    vecTriggerID.push_back("600031");
  }

  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mTriggerId[2]; // 0: before cuts | 1: after cuts

  for(int iCut = 0; iCut < 2; ++iCut)
  {
    string HistName = Form("h_mTriggerId%s",CutStatus[iCut].c_str());
    h_mTriggerId[iCut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mTriggerId[iCut]->SetLineColor(iCut+1);
    h_mTriggerId[iCut]->GetXaxis()->SetTitle("triggerId");
    for(int iTrigger = 0; iTrigger < vecTriggerID.size(); ++iTrigger)
    {
      h_mTriggerId[iCut]->GetXaxis()->SetBinLabel(iTrigger+1,vecTriggerID[iTrigger].c_str());
    }
  }

  TCanvas *c_TriggerId = new TCanvas("c_TriggerId","c_TriggerId",10,10,1600,800);
  c_TriggerId->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_TriggerId->cd(i_pad+1);
    c_TriggerId->cd(i_pad+1)->SetLeftMargin(0.1);
    c_TriggerId->cd(i_pad+1)->SetRightMargin(0.1);
    c_TriggerId->cd(i_pad+1)->SetBottomMargin(0.1);
    c_TriggerId->cd(i_pad+1)->SetGrid(0,0);
    c_TriggerId->cd(i_pad+1)->SetTicks(1,1);
    c_TriggerId->cd(i_pad+1)->SetLogy(1);
  }
  for(int iCut = 0; iCut < 2; ++iCut)
  {
    c_TriggerId->cd(iCut+1);
    h_mTriggerId[iCut]->Draw();
  }

  string FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_TriggerId_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_TriggerId->SaveAs(FigName.c_str());
}
