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

void rejectPileUpEvent(int energy = 0)
{
  string JobId = "low";
  string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/test/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH2F *h_mRefMultGRefMult[2];
  TH2F *h_mRefMultTofMatch[2];
  TH2F *h_mRefMultTofHits[2];
  TH2F *h_mGRefMultTofMatch[2];
  TH2F *h_mGRefMultTofHits[2];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    // string HistName = Form("h_mRefMultGRefMult%s_trigger9",mCutsQA[i_cut].c_str());
    string HistName = Form("h_mRefMultGRefMult_%s",CutsQA[i_cut].c_str());
    h_mRefMultGRefMult[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mRefMultGRefMult[i_cut]->GetXaxis()->SetTitle("refMult");
    h_mRefMultGRefMult[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mRefMultGRefMult[i_cut]->GetYaxis()->SetTitle("grefMult");
    h_mRefMultGRefMult[i_cut]->GetYaxis()->SetRangeUser(0.0,800.0);

    HistName = Form("h_mRefMultTofMatch_%s",CutsQA[i_cut].c_str());
    h_mRefMultTofMatch[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mRefMultTofMatch[i_cut]->GetXaxis()->SetTitle("refMult");
    h_mRefMultTofMatch[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mRefMultTofMatch[i_cut]->GetYaxis()->SetTitle("ToFMatch");
    h_mRefMultTofMatch[i_cut]->GetYaxis()->SetRangeUser(0.0,1000.0);

    HistName = Form("h_mRefMultTofHits_%s",CutsQA[i_cut].c_str());
    h_mRefMultTofHits[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mRefMultTofHits[i_cut]->GetXaxis()->SetTitle("refMult");
    h_mRefMultTofHits[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mRefMultTofHits[i_cut]->GetYaxis()->SetTitle("ToFHits");
    h_mRefMultTofHits[i_cut]->GetYaxis()->SetRangeUser(0.0,3500.0);

    HistName = Form("h_mGRefMultTofMatch_%s",CutsQA[i_cut].c_str());
    h_mGRefMultTofMatch[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mGRefMultTofMatch[i_cut]->GetXaxis()->SetTitle("gRefMult");
    h_mGRefMultTofMatch[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mGRefMultTofMatch[i_cut]->GetYaxis()->SetTitle("ToFMatch");
    h_mGRefMultTofMatch[i_cut]->GetYaxis()->SetRangeUser(0.0,1000.0);

    HistName = Form("h_mGRefMultTofHits_%s",CutsQA[i_cut].c_str());
    h_mGRefMultTofHits[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mGRefMultTofHits[i_cut]->GetXaxis()->SetTitle("gRefMult");
    h_mGRefMultTofHits[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mGRefMultTofHits[i_cut]->GetYaxis()->SetTitle("ToFHits");
    h_mGRefMultTofHits[i_cut]->GetYaxis()->SetRangeUser(0.0,3500.0);
  }

  TCanvas *c_EventQA = new TCanvas("c_EventQA","c_EventQA",10,10,1500,600);
  c_EventQA->Divide(5,2);
  for(int i_pad = 0; i_pad < 10; ++i_pad)
  {
    c_EventQA->cd(i_pad+1);
    c_EventQA->cd(i_pad+1)->SetLeftMargin(0.1);
    c_EventQA->cd(i_pad+1)->SetRightMargin(0.1);
    c_EventQA->cd(i_pad+1)->SetBottomMargin(0.1);
    c_EventQA->cd(i_pad+1)->SetGrid(0,0);
    c_EventQA->cd(i_pad+1)->SetTicks(1,1);
    c_EventQA->cd(i_pad+1)->SetLogz(1);
  }

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    c_EventQA->cd(i_cut*5+1);
    h_mRefMultGRefMult[i_cut]->Draw("colz");

    c_EventQA->cd(i_cut*5+2);
    h_mRefMultTofMatch[i_cut]->Draw("colz");

    c_EventQA->cd(i_cut*5+3);
    h_mRefMultTofHits[i_cut]->Draw("colz");

    c_EventQA->cd(i_cut*5+4);
    h_mGRefMultTofMatch[i_cut]->Draw("colz");

    c_EventQA->cd(i_cut*5+5);
    h_mGRefMultTofHits[i_cut]->Draw("colz");
  }

  TCanvas *c_GRefMult_nBToFMatch = new TCanvas("c_GRefMult_nBToFMatch","c_GRefMult_nBToFMatch",10,10,900,900);
  c_GRefMult_nBToFMatch->Divide(3,3);
  TH1D *h_mGRefMultTofMatchProj_0 = (TH1D*)h_mGRefMultTofMatch[1]->ProjectionX("h_mGRefMultTofMatchProj_0",3,5);
  c_GRefMult_nBToFMatch->cd(1);
  c_GRefMult_nBToFMatch->cd(1)->SetLogy(1);
  h_mGRefMultTofMatchProj_0->Draw("hE");
}
