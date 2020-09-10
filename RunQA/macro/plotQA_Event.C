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

void plotQA_Event(int energy = 0)
{
  string JobId = "45056E11E5EFF5990EADD3B54A0E4FB6";
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/test/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mRefMult[2]; // 0: before cuts | 1: after cuts
  TH1F *h_mGRefMult[2]; // 0: before cuts | 1: after cuts
  TH1F *h_mCentrality9[2];
  TH2F *h_mRefMultTofMatch[2];
  TH2F *h_mRefMultTofHits[2];
  TH2F *h_mGRefMultTofMatch[2];
  TH2F *h_mGRefMultTofHits[2];
  TH2F *h_mVzVzVpd[2];
  TH1F *h_mDiffVzVzVpd[2];
  TH1F *h_mVertexZ[2];
  TH2F *h_mVertexXY[2];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    string HistName = Form("h_mRefMult_%s",CutsQA[i_cut].c_str());
    h_mRefMult[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mRefMult[i_cut]->SetLineColor(i_cut+1);
    h_mRefMult[i_cut]->GetXaxis()->SetTitle("refMult");

    string HistName = Form("h_mGRefMult_%s",CutsQA[i_cut].c_str());
    h_mGRefMult[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mGRefMult[i_cut]->SetLineColor(i_cut+1);
    h_mGRefMult[i_cut]->GetXaxis()->SetTitle("gRefMult");

    HistName = Form("h_mCentrality9_%s",CutsQA[i_cut].c_str());
    h_mCentrality9[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mCentrality9[i_cut]->SetLineColor(i_cut+1);
    h_mCentrality9[i_cut]->GetXaxis()->SetTitle("centrality");

    HistName = Form("h_mRefMultTofMatch_%s",CutsQA[i_cut].c_str());
    h_mRefMultTofMatch[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mRefMultTofMatch[i_cut]->GetXaxis()->SetTitle("refMult");
    h_mRefMultTofMatch[i_cut]->GetYaxis()->SetTitle("ToFMatch");

    HistName = Form("h_mRefMultTofHits_%s",CutsQA[i_cut].c_str());
    h_mRefMultTofHits[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mRefMultTofHits[i_cut]->GetXaxis()->SetTitle("refMult");
    h_mRefMultTofHits[i_cut]->GetYaxis()->SetTitle("ToFHits");

    HistName = Form("h_mGRefMultTofMatch_%s",CutsQA[i_cut].c_str());
    h_mGRefMultTofMatch[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mGRefMultTofMatch[i_cut]->GetXaxis()->SetTitle("gRefMult");
    h_mGRefMultTofMatch[i_cut]->GetYaxis()->SetTitle("ToFMatch");

    HistName = Form("h_mGRefMultTofHits_%s",CutsQA[i_cut].c_str());
    h_mGRefMultTofHits[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mGRefMultTofHits[i_cut]->GetXaxis()->SetTitle("gRefMult");
    h_mGRefMultTofHits[i_cut]->GetYaxis()->SetTitle("ToFHits");

    HistName = Form("h_mVertexXY_%s",CutsQA[i_cut].c_str());
    h_mVertexXY[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mVertexXY[i_cut]->GetXaxis()->SetTitle("vx");
    h_mVertexXY[i_cut]->GetYaxis()->SetTitle("vy");

    HistName = Form("h_mVertexZ_%s",CutsQA[i_cut].c_str());
    h_mVertexZ[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mVertexZ[i_cut]->GetXaxis()->SetTitle("vz");
    h_mVertexZ[i_cut]->SetLineColor(i_cut+1);

    HistName = Form("h_mVzVzVpd_%s",CutsQA[i_cut].c_str());
    h_mVzVzVpd[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mVzVzVpd[i_cut]->GetXaxis()->SetTitle("vz");
    h_mVzVzVpd[i_cut]->GetYaxis()->SetTitle("vzVpd");

    HistName = Form("h_mDiffVzVzVpd_%s",CutsQA[i_cut].c_str());
    h_mDiffVzVzVpd[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mDiffVzVzVpd[i_cut]->SetLineColor(i_cut+1);
    h_mDiffVzVzVpd[i_cut]->GetXaxis()->SetTitle("vz-vzVpd");
  }

  TF1 *f_tofHitsCut_low = new TF1("f_tofHitsCut_low","2.88*x-155",0,800);
  f_tofHitsCut_low->SetLineColor(2);
  f_tofHitsCut_low->SetLineWidth(2);
  f_tofHitsCut_low->SetLineStyle(2);

  TCanvas *c_EventQA[2];
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    string CanName = Form("c_EventQA_%s",CutsQA[i_cut].c_str());
    c_EventQA[i_cut] = new TCanvas(CanName.c_str(),CanName.c_str(),10,10,1600,800);
    c_EventQA[i_cut]->Divide(4,2);
    for(int i_pad = 0; i_pad < 8; ++i_pad)
    {
      c_EventQA[i_cut]->cd(i_pad+1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetLeftMargin(0.1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetRightMargin(0.1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetBottomMargin(0.1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetGrid(0,0);
      c_EventQA[i_cut]->cd(i_pad+1)->SetTicks(1,1);
    }

    c_EventQA[i_cut]->cd(1);
    c_EventQA[i_cut]->cd(1)->SetLogy();
    if(energy == 0)h_mGRefMult[i_cut]->Draw("hE");
    if(energy != 0)h_mRefMult[i_cut]->Draw("hE");

    c_EventQA[i_cut]->cd(2);
    c_EventQA[i_cut]->cd(2)->SetLogz();
    if(energy == 0)
    {
      h_mGRefMultTofMatch[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
      h_mGRefMultTofMatch[i_cut]->GetYaxis()->SetRangeUser(0.0,800.0);
      h_mGRefMultTofMatch[i_cut]->Draw("colz");
    }
    if(energy != 0)
    {
      h_mRefMultTofMatch[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
      h_mRefMultTofMatch[i_cut]->GetYaxis()->SetRangeUser(0.0,800.0);
      h_mRefMultTofMatch[i_cut]->Draw("colz");
    }

    c_EventQA[i_cut]->cd(3);
    c_EventQA[i_cut]->cd(3)->SetLogz();
    if(energy == 0)
    {
      h_mGRefMultTofHits[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
      h_mGRefMultTofHits[i_cut]->Draw("colz");
    }
    if(energy != 0)
    {
      h_mRefMultTofHits[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
      h_mRefMultTofHits[i_cut]->Draw("colz");
    }

    c_EventQA[i_cut]->cd(4);
    c_EventQA[i_cut]->cd(4)->SetLogy();
    h_mCentrality9[i_cut]->Draw("hE");

    c_EventQA[i_cut]->cd(5);
    c_EventQA[i_cut]->cd(5)->SetLogz();
    h_mVertexXY[i_cut]->Draw("colz");

    c_EventQA[i_cut]->cd(6);
    // c_EventQA[i_cut]->cd(6)->SetLogy();
    h_mVertexZ[i_cut]->Draw();

    c_EventQA[i_cut]->cd(7);
    c_EventQA[i_cut]->cd(7)->SetLogz();
    h_mVzVzVpd[i_cut]->Draw("colz");

    c_EventQA[i_cut]->cd(8);
    // c_EventQA[i_cut]->cd(8)->SetLogy();
    h_mDiffVzVzVpd[i_cut]->Draw();

    string FigName = Form("c_EventQA_%s_%s_%s.png",CutsQA[i_cut].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
    c_EventQA[i_cut]->SaveAs(FigName.c_str());
  }
}
