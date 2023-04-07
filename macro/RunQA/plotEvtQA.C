#include <string>
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
static const int numCuts = 2; // 0: before cuts | 1: after cuts
static const int numTriggerBins = 10; // 0-8 for different triggerID | 9 for all triggers

void plotEvtQA(int beamType = 2)
{
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/RunQA/%s/file_%s_RunQA.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mRefMult[numCuts][numTriggerBins];
  TH1F *h_mGRefMult[numCuts][numTriggerBins];
  TH2F *h_mRefMultGRefMult[numCuts][numTriggerBins];
  TH1F *h_mCentrality9[numCuts][numTriggerBins];
  TH2F *h_mTofMatchRefMult[numCuts][numTriggerBins];
  TH2F *h_mTofHitsRefMult[numCuts][numTriggerBins];
  TH2F *h_mTofMatchGRefMult[numCuts][numTriggerBins];
  TH2F *h_mTofHitsGRefMult[numCuts][numTriggerBins];
  TH2F *h_mVzVzVpd[numCuts][numTriggerBins];
  TH1F *h_mDiffVzVzVpd[numCuts][numTriggerBins];
  TH1F *h_mVertexZ[numCuts][numTriggerBins];
  TH2F *h_mVertexXY[numCuts][numTriggerBins];

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < numTriggerBins; ++iTrig)
    {
      string HistName = Form("h_mRefMult%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mRefMult[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mRefMult[iCut][iTrig]->SetLineColor(iCut+1);
      h_mRefMult[iCut][iTrig]->GetXaxis()->SetTitle("refMult");

      HistName = Form("h_mGRefMult%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mGRefMult[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());;
      h_mGRefMult[iCut][iTrig]->SetLineColor(iCut+1);
      h_mGRefMult[iCut][iTrig]->GetXaxis()->SetTitle("gRefMult");

      HistName = Form("h_mRefMultGRefMult%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mRefMultGRefMult[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mRefMultGRefMult[iCut][iTrig]->GetXaxis()->SetTitle("refMult");
      h_mRefMultGRefMult[iCut][iTrig]->GetYaxis()->SetTitle("gRefMult");

      HistName = Form("h_mCentrality9%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mCentrality9[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mCentrality9[iCut][iTrig]->SetLineColor(iCut+1);
      h_mCentrality9[iCut][iTrig]->GetXaxis()->SetTitle("centrality");

      HistName = Form("h_mTofMatchRefMult%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mTofMatchRefMult[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mTofMatchRefMult[iCut][iTrig]->GetXaxis()->SetTitle("tofMatch");
      h_mTofMatchRefMult[iCut][iTrig]->GetYaxis()->SetTitle("refMult");

      HistName = Form("h_mTofHitsRefMult%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mTofHitsRefMult[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mTofHitsRefMult[iCut][iTrig]->GetXaxis()->SetTitle("tofHits");
      h_mTofHitsRefMult[iCut][iTrig]->GetYaxis()->SetTitle("refMult");

      HistName = Form("h_mTofMatchGRefMult%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mTofMatchGRefMult[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());;
      h_mTofMatchGRefMult[iCut][iTrig]->GetXaxis()->SetTitle("tofMatch");
      h_mTofMatchGRefMult[iCut][iTrig]->GetYaxis()->SetTitle("gRefMult");

      HistName = Form("h_mTofHitsGRefMult%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mTofHitsGRefMult[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mTofHitsGRefMult[iCut][iTrig]->GetXaxis()->SetTitle("tofHits");
      h_mTofHitsGRefMult[iCut][iTrig]->GetYaxis()->SetTitle("gRefMult");

      HistName = Form("h_mVertexXY%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mVertexXY[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mVertexXY[iCut][iTrig]->GetXaxis()->SetTitle("Vx");
      h_mVertexXY[iCut][iTrig]->GetYaxis()->SetTitle("Vy");

      HistName = Form("h_mVertexZ%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mVertexZ[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mVertexZ[iCut][iTrig]->SetLineColor(iCut+1);
      h_mVertexZ[iCut][iTrig]->GetXaxis()->SetTitle("Vz");

      HistName = Form("h_mVzVzVpd%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mVzVzVpd[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());
      h_mVzVzVpd[iCut][iTrig]->GetXaxis()->SetTitle("Vz");
      h_mVzVzVpd[iCut][iTrig]->GetYaxis()->SetTitle("VzVpd");

      HistName = Form("h_mDiffVzVzVpd%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mDiffVzVzVpd[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());
      h_mDiffVzVzVpd[iCut][iTrig]->SetLineColor(iCut+1);
      h_mDiffVzVzVpd[iCut][iTrig]->GetXaxis()->SetTitle("Vz-VzVpd");
    }
  }

  TCanvas *c_EventQA = new TCanvas("c_EventQA","c_EventQA",10,10,1600,800);
  c_EventQA->Divide(4,2);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_EventQA->cd(i_pad+1);
    c_EventQA->cd(i_pad+1)->SetLeftMargin(0.1);
    c_EventQA->cd(i_pad+1)->SetRightMargin(0.1);
    c_EventQA->cd(i_pad+1)->SetBottomMargin(0.1);
    c_EventQA->cd(i_pad+1)->SetGrid(0,0);
    c_EventQA->cd(i_pad+1)->SetTicks(1,1);
  }
  std::string figName = Form("../../figures/RunQA/%s/EvtQA_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_EventQA->Print(figName.c_str());

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_EventQA->cd(1);
    c_EventQA->cd(1)->Clear();
    c_EventQA->cd(1)->SetLogy();
    h_mRefMult[iCut][9]->GetXaxis()->SetRangeUser(0.0,500.0);
    if(beamType == 2) h_mRefMult[iCut][9]->GetXaxis()->SetRangeUser(0.0,300.0);
    h_mRefMult[iCut][9]->Draw("hE");

    c_EventQA->cd(2);
    c_EventQA->cd(2)->Clear();
    c_EventQA->cd(2)->SetLogz();
    h_mTofMatchRefMult[iCut][9]->GetXaxis()->SetRangeUser(0.0,500.0);
    h_mTofMatchRefMult[iCut][9]->GetYaxis()->SetRangeUser(0.0,500.0);
    if(beamType == 2) h_mTofMatchRefMult[iCut][9]->GetXaxis()->SetRangeUser(0.0,300.0);
    if(beamType == 2) h_mTofMatchRefMult[iCut][9]->GetYaxis()->SetRangeUser(0.0,300.0);
    h_mTofMatchRefMult[iCut][9]->Draw("colz");

    c_EventQA->cd(3);
    c_EventQA->cd(3)->Clear();
    c_EventQA->cd(3)->SetLogz();
    h_mTofHitsRefMult[iCut][9]->GetXaxis()->SetRangeUser(0.0,500.0);
    h_mTofHitsRefMult[iCut][9]->GetYaxis()->SetRangeUser(0.0,500.0);
    if(beamType == 2) h_mTofHitsRefMult[iCut][9]->GetXaxis()->SetRangeUser(0.0,300.0);
    if(beamType == 2) h_mTofHitsRefMult[iCut][9]->GetYaxis()->SetRangeUser(0.0,300.0);
    h_mTofHitsRefMult[iCut][9]->Draw("colz");

    c_EventQA->cd(4);
    c_EventQA->cd(4)->Clear();
    c_EventQA->cd(4)->SetLogy();
    h_mCentrality9[iCut][9]->Draw("hE");

    c_EventQA->cd(5);
    c_EventQA->cd(5)->Clear();
    c_EventQA->cd(5)->SetLogz();
    h_mVertexXY[iCut][9]->Draw("colz");

    c_EventQA->cd(6);
    c_EventQA->cd(6)->Clear();
    c_EventQA->cd(6)->SetLogy();
    h_mVertexZ[iCut][9]->Draw();

    c_EventQA->cd(7);
    c_EventQA->cd(7)->Clear();
    c_EventQA->cd(7)->SetLogz();
    h_mVzVzVpd[iCut][9]->Draw("colz");

    c_EventQA->cd(8);
    c_EventQA->cd(8)->Clear();
    // c_EventQA->cd(8)->SetLogy();
    h_mDiffVzVzVpd[iCut][9]->Draw();

    figName = Form("../../figures/RunQA/%s/EvtQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EventQA->Update();
    c_EventQA->Print(figName.c_str());
  }

  figName = Form("../../figures/RunQA/%s/EvtQA_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_EventQA->Print(figName.c_str());
}
