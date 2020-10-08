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

static const string mCutsQA[2] = {"Before","After"};

void plotQA_Event(int energy = 0)
{
  string JobId = "EEE0479FEE171BB7ACDE3FBF146413E7";
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/test/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mRefMult[2][10]; // 0: before cuts | 1: after cuts
  TH1F *h_mGRefMult[2][10]; // 0-8 for different triggerID | 9 for all triggers
  TH2F *h_mRefMultGRefMult[2][10];
  TH1F *h_mCentrality9[2][10];
  TH2F *h_mTofMatchRefMult[2][10];
  TH2F *h_mTofHitsRefMult[2][10];
  TH2F *h_mTofMatchGRefMult[2][10];
  TH2F *h_mTofHitsGRefMult[2][10];
  TH2F *h_mVzVzVpd[2][10];
  TH1F *h_mDiffVzVzVpd[2][10];
  TH1F *h_mVertexZ[2][10];
  TH2F *h_mVertexXY[2][10];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      std::string HistName = Form("h_mRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mRefMult[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mRefMult[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("refMult");

      HistName = Form("h_mGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mGRefMult[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());;
      h_mGRefMult[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mGRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("gRefMult");

      HistName = Form("h_mRefMultGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mRefMultGRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mRefMultGRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("refMult");
      h_mRefMultGRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("gRefMult");

      HistName = Form("h_mCentrality9%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mCentrality9[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mCentrality9[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mCentrality9[i_cut][i_trig]->GetXaxis()->SetTitle("centrality");

      HistName = Form("h_mTofMatchRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofMatchRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mTofMatchRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("tofMatch");
      h_mTofMatchRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("refMult");

      HistName = Form("h_mTofHitsRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofHitsRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mTofHitsRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("tofHits");
      h_mTofHitsRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("refMult");

      HistName = Form("h_mTofMatchGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofMatchGRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());;
      h_mTofMatchGRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("tofMatch");
      h_mTofMatchGRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("gRefMult");

      HistName = Form("h_mTofHitsGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofHitsGRefMult[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mTofHitsGRefMult[i_cut][i_trig]->GetXaxis()->SetTitle("tofHits");
      h_mTofHitsGRefMult[i_cut][i_trig]->GetYaxis()->SetTitle("gRefMult");

      HistName = Form("h_mVertexXY%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mVertexXY[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mVertexXY[i_cut][i_trig]->GetXaxis()->SetTitle("Vx");
      h_mVertexXY[i_cut][i_trig]->GetYaxis()->SetTitle("Vy");

      HistName = Form("h_mVertexZ%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mVertexZ[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mVertexZ[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mVertexZ[i_cut][i_trig]->GetXaxis()->SetTitle("Vz");

      HistName = Form("h_mVzVzVpd%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mVzVzVpd[i_cut][i_trig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mVzVzVpd[i_cut][i_trig]->GetXaxis()->SetTitle("Vz");
      h_mVzVzVpd[i_cut][i_trig]->GetYaxis()->SetTitle("VzVpd");

      HistName = Form("h_mDiffVzVzVpd%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mDiffVzVzVpd[i_cut][i_trig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mDiffVzVzVpd[i_cut][i_trig]->SetLineColor(i_cut+1);
      h_mDiffVzVzVpd[i_cut][i_trig]->GetXaxis()->SetTitle("Vz-VzVpd");
    }
  }

  TCanvas *c_EventQA[2];
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    string CanName = Form("c_EventQA_%s",mCutsQA[i_cut].c_str());
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
    if(energy == 0)h_mRefMult[i_cut][9]->Draw("hE");
    if(energy != 0)h_mRefMult[i_cut][9]->Draw("hE");

    c_EventQA[i_cut]->cd(2);
    c_EventQA[i_cut]->cd(2)->SetLogz();
    if(energy == 0)
    {
      h_mTofMatchRefMult[i_cut][9]->GetXaxis()->SetRangeUser(0.0,800.0);
      h_mTofMatchRefMult[i_cut][9]->GetYaxis()->SetRangeUser(0.0,800.0);
      h_mTofMatchRefMult[i_cut][9]->Draw("colz");
    }
    if(energy != 0)
    {
      h_mTofMatchRefMult[i_cut][9]->GetXaxis()->SetRangeUser(0.0,800.0);
      h_mTofMatchRefMult[i_cut][9]->GetYaxis()->SetRangeUser(0.0,800.0);
      h_mTofMatchRefMult[i_cut][9]->Draw("colz");
    }

    c_EventQA[i_cut]->cd(3);
    c_EventQA[i_cut]->cd(3)->SetLogz();
    if(energy == 0)
    {
      h_mTofHitsRefMult[i_cut][9]->GetYaxis()->SetRangeUser(0.0,800.0);
      h_mTofHitsRefMult[i_cut][9]->Draw("colz");
    }
    if(energy != 0)
    {
      h_mTofHitsRefMult[i_cut][9]->GetYaxis()->SetRangeUser(0.0,800.0);
      h_mTofHitsRefMult[i_cut][9]->Draw("colz");
    }

    c_EventQA[i_cut]->cd(4);
    c_EventQA[i_cut]->cd(4)->SetLogy();
    h_mCentrality9[i_cut][9]->Draw("hE");

    c_EventQA[i_cut]->cd(5);
    c_EventQA[i_cut]->cd(5)->SetLogz();
    h_mVertexXY[i_cut][9]->Draw("colz");

    c_EventQA[i_cut]->cd(6);
    // c_EventQA[i_cut]->cd(6)->SetLogy();
    h_mVertexZ[i_cut][9]->Draw();

    c_EventQA[i_cut]->cd(7);
    c_EventQA[i_cut]->cd(7)->SetLogz();
    h_mVzVzVpd[i_cut][9]->Draw("colz");

    c_EventQA[i_cut]->cd(8);
    // c_EventQA[i_cut]->cd(8)->SetLogy();
    h_mDiffVzVzVpd[i_cut][9]->Draw();

    string FigName = Form("./figures/c_EventQA_%s_%s_%s.png",mCutsQA[i_cut].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
    c_EventQA[i_cut]->SaveAs(FigName.c_str());
  }
}
