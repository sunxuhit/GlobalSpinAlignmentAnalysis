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

void plotQA_RunbyRun(int energy = 0)
{
  string JobId = "EEE0479FEE171BB7ACDE3FBF146413E7";
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/test/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TProfile *p_mRefMult[2][10]; // 0: before cuts | 1: after cuts
  TProfile *p_mGRefMult[2][10]; // 0-8 for different triggerID | 9 for all triggers
  TProfile *p_mZdcX[2][10];
  TProfile *p_mVz[2][10];
  TProfile *p_mVr[2][10];

  TProfile *p_mGDca[2][10];
  TProfile *p_mNHitsFit[2][10];
  TProfile *p_mPrimPt[2][10];
  TProfile *p_mPrimEta[2][10];
  TProfile *p_mPrimPhi[2][10];
  TProfile *p_mGlobPt[2][10];
  TProfile *p_mGlobEta[2][10];
  TProfile *p_mGlobPhi[2][10];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      std::string ProName = Form("p_mRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mRefMult[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGRefMult[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mZdcX%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mZdcX[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mVz%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mVz[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mVr%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mVr[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGDca%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGDca[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mNHitsFit%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsFit[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimPt[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimEta[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimPhi[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGlobPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobPt[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGlobEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobEta[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGlobPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobPhi[i_cut][i_trig] = (TProfile*)File_InPut->Get(ProName.c_str());
    }
  }

  const int numOfTriggers = 5;
  const int triggerID[numOfTriggers] = {450005,450015,450025,450050,450060};
  const int MarkerColor[numOfTriggers] = {7,6,4,2,16};

  TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);

  //---------------------
  TCanvas *c_RunQA_RefMult = new TCanvas("c_RunQA_RefMult","c_RunQA_RefMult",10,10,800,400);
  c_RunQA_RefMult->cd()->SetLeftMargin(0.1);
  c_RunQA_RefMult->cd()->SetRightMargin(0.1);
  c_RunQA_RefMult->cd()->SetBottomMargin(0.1);
  c_RunQA_RefMult->cd()->SetGrid(0,0);
  c_RunQA_RefMult->cd()->SetTicks(1,1);

  p_mRefMult[1][9]->SetTitle("refMult vs. runIndex");
  p_mRefMult[1][9]->SetStats(0);
  p_mRefMult[1][9]->SetMarkerColor(1);
  p_mRefMult[1][9]->SetMarkerStyle(20);
  p_mRefMult[1][9]->SetMarkerSize(1.0);
  p_mRefMult[1][9]->SetLineColor(1);
  p_mRefMult[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mRefMult[1][9]->GetYaxis()->SetTitle("<refMult>");
  p_mRefMult[1][9]->GetYaxis()->SetRangeUser(0,400);
  p_mRefMult[1][9]->Draw("pE");
  leg->AddEntry(p_mRefMult[1][9],"All Triggers","P");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mRefMult[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mRefMult[1][i_trig]->SetMarkerStyle(24);
    p_mRefMult[1][i_trig]->SetMarkerSize(0.8);
    p_mRefMult[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mRefMult[1][i_trig]->GetEntries() > 0) p_mRefMult[1][i_trig]->Draw("pE same");
    leg->AddEntry(p_mRefMult[1][i_trig],to_string(triggerID[i_trig]).c_str(),"P");
  }

  leg->Draw("same");
  c_RunQA_RefMult->SaveAs("./figures/c_RunQA_RefMult.eps");
  //---------------------
  TCanvas *c_RunQA_gRefMult = new TCanvas("c_RunQA_gRefMult","c_RunQA_gRefMult",10,10,800,400);
  c_RunQA_gRefMult->cd()->SetLeftMargin(0.1);
  c_RunQA_gRefMult->cd()->SetRightMargin(0.1);
  c_RunQA_gRefMult->cd()->SetBottomMargin(0.1);
  c_RunQA_gRefMult->cd()->SetGrid(0,0);
  c_RunQA_gRefMult->cd()->SetTicks(1,1);

  p_mGRefMult[1][9]->SetTitle("grefMult vs. runIndex");
  p_mGRefMult[1][9]->SetStats(0);
  p_mGRefMult[1][9]->SetMarkerColor(1);
  p_mGRefMult[1][9]->SetMarkerStyle(20);
  p_mGRefMult[1][9]->SetMarkerSize(1.0);
  p_mGRefMult[1][9]->SetLineColor(1);
  p_mGRefMult[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGRefMult[1][9]->GetYaxis()->SetTitle("<grefMult>");
  p_mGRefMult[1][9]->GetYaxis()->SetRangeUser(0,400);
  p_mGRefMult[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGRefMult[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGRefMult[1][i_trig]->SetMarkerStyle(24);
    p_mGRefMult[1][i_trig]->SetMarkerSize(0.8);
    p_mGRefMult[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGRefMult[1][i_trig]->GetEntries() > 0) p_mGRefMult[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_gRefMult->SaveAs("./figures/c_RunQA_gRefMult.eps");
  //---------------------
  TCanvas *c_RunQA_ZdcX = new TCanvas("c_RunQA_ZdcX","c_RunQA_ZdcX",10,10,800,400);
  c_RunQA_ZdcX->cd()->SetLeftMargin(0.1);
  c_RunQA_ZdcX->cd()->SetRightMargin(0.1);
  c_RunQA_ZdcX->cd()->SetBottomMargin(0.1);
  c_RunQA_ZdcX->cd()->SetGrid(0,0);
  c_RunQA_ZdcX->cd()->SetTicks(1,1);

  p_mZdcX[1][9]->SetTitle("ZdcX vs. runIndex");
  p_mZdcX[1][9]->SetStats(0);
  p_mZdcX[1][9]->SetMarkerColor(1);
  p_mZdcX[1][9]->SetMarkerStyle(20);
  p_mZdcX[1][9]->SetMarkerSize(1.0);
  p_mZdcX[1][9]->SetLineColor(1);
  p_mZdcX[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mZdcX[1][9]->GetYaxis()->SetTitle("<ZdcX>");
  // p_mZdcX[1][9]->GetYaxis()->SetRangeUser(0,400);
  p_mZdcX[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mZdcX[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mZdcX[1][i_trig]->SetMarkerStyle(24);
    p_mZdcX[1][i_trig]->SetMarkerSize(0.8);
    p_mZdcX[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mZdcX[1][i_trig]->GetEntries() > 0) p_mZdcX[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_ZdcX->SaveAs("./figures/c_RunQA_ZdcX.eps");
  //---------------------
  TCanvas *c_RunQA_Vz = new TCanvas("c_RunQA_Vz","c_RunQA_Vz",10,10,800,400);
  c_RunQA_Vz->cd()->SetLeftMargin(0.1);
  c_RunQA_Vz->cd()->SetRightMargin(0.1);
  c_RunQA_Vz->cd()->SetBottomMargin(0.1);
  c_RunQA_Vz->cd()->SetGrid(0,0);
  c_RunQA_Vz->cd()->SetTicks(1,1);

  p_mVz[1][9]->SetTitle("Vz vs. runIndex");
  p_mVz[1][9]->SetStats(0);
  p_mVz[1][9]->SetMarkerColor(1);
  p_mVz[1][9]->SetMarkerStyle(20);
  p_mVz[1][9]->SetMarkerSize(1.0);
  p_mVz[1][9]->SetLineColor(1);
  p_mVz[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mVz[1][9]->GetYaxis()->SetTitle("<Vz>");
  p_mVz[1][9]->GetYaxis()->SetRangeUser(-3.0,3.0);
  p_mVz[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mVz[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mVz[1][i_trig]->SetMarkerStyle(24);
    p_mVz[1][i_trig]->SetMarkerSize(0.8);
    p_mVz[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mVz[1][i_trig]->GetEntries() > 0) p_mVz[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_Vz->SaveAs("./figures/c_RunQA_Vz.eps");
  //---------------------
  TCanvas *c_RunQA_Vr = new TCanvas("c_RunQA_Vr","c_RunQA_Vr",10,10,800,400);
  c_RunQA_Vr->cd()->SetLeftMargin(0.1);
  c_RunQA_Vr->cd()->SetRightMargin(0.1);
  c_RunQA_Vr->cd()->SetBottomMargin(0.1);
  c_RunQA_Vr->cd()->SetGrid(0,0);
  c_RunQA_Vr->cd()->SetTicks(1,1);

  p_mVr[1][9]->SetTitle("Vr vs. runIndex");
  p_mVr[1][9]->SetStats(0);
  p_mVr[1][9]->SetMarkerColor(1);
  p_mVr[1][9]->SetMarkerStyle(20);
  p_mVr[1][9]->SetMarkerSize(1.0);
  p_mVr[1][9]->SetLineColor(1);
  p_mVr[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mVr[1][9]->GetYaxis()->SetTitle("<Vr>");
  p_mVr[1][9]->GetYaxis()->SetRangeUser(0.1,0.5);
  p_mVr[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mVr[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mVr[1][i_trig]->SetMarkerStyle(24);
    p_mVr[1][i_trig]->SetMarkerSize(0.8);
    p_mVr[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mVr[1][i_trig]->GetEntries() > 0) p_mVr[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_Vr->SaveAs("./figures/c_RunQA_Vr.eps");
  //---------------------
  TCanvas *c_RunQA_gDca = new TCanvas("c_RunQA_gDca","c_RunQA_gDca",10,10,800,400);
  c_RunQA_gDca->cd()->SetLeftMargin(0.1);
  c_RunQA_gDca->cd()->SetRightMargin(0.1);
  c_RunQA_gDca->cd()->SetBottomMargin(0.1);
  c_RunQA_gDca->cd()->SetGrid(0,0);
  c_RunQA_gDca->cd()->SetTicks(1,1);

  p_mGDca[1][9]->SetTitle("gDca vs. runIndex");
  p_mGDca[1][9]->SetStats(0);
  p_mGDca[1][9]->SetMarkerColor(1);
  p_mGDca[1][9]->SetMarkerStyle(20);
  p_mGDca[1][9]->SetMarkerSize(1.0);
  p_mGDca[1][9]->SetLineColor(1);
  p_mGDca[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGDca[1][9]->GetYaxis()->SetTitle("<gDca>");
  p_mGDca[1][9]->GetYaxis()->SetRangeUser(0.2,0.6);
  p_mGDca[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGDca[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGDca[1][i_trig]->SetMarkerStyle(24);
    p_mGDca[1][i_trig]->SetMarkerSize(0.8);
    p_mGDca[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGDca[1][i_trig]->GetEntries() > 0) p_mGDca[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_gDca->SaveAs("./figures/c_RunQA_gDca.eps");
  //---------------------
  TCanvas *c_RunQA_nHitsFit = new TCanvas("c_RunQA_nHitsFit","c_RunQA_nHitsFit",10,10,800,400);
  c_RunQA_nHitsFit->cd()->SetLeftMargin(0.1);
  c_RunQA_nHitsFit->cd()->SetRightMargin(0.1);
  c_RunQA_nHitsFit->cd()->SetBottomMargin(0.1);
  c_RunQA_nHitsFit->cd()->SetGrid(0,0);
  c_RunQA_nHitsFit->cd()->SetTicks(1,1);

  p_mNHitsFit[1][9]->SetTitle("nHitsFit vs. runIndex");
  p_mNHitsFit[1][9]->SetStats(0);
  p_mNHitsFit[1][9]->SetMarkerColor(1);
  p_mNHitsFit[1][9]->SetMarkerStyle(20);
  p_mNHitsFit[1][9]->SetMarkerSize(1.0);
  p_mNHitsFit[1][9]->SetLineColor(1);
  p_mNHitsFit[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mNHitsFit[1][9]->GetYaxis()->SetTitle("<nHitsFit>");
  p_mNHitsFit[1][9]->GetYaxis()->SetRangeUser(25,40);
  p_mNHitsFit[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mNHitsFit[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mNHitsFit[1][i_trig]->SetMarkerStyle(24);
    p_mNHitsFit[1][i_trig]->SetMarkerSize(0.8);
    p_mNHitsFit[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mNHitsFit[1][i_trig]->GetEntries() > 0) p_mNHitsFit[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_nHitsFit->SaveAs("./figures/c_RunQA_nHitsFit.eps");
  //---------------------
  TCanvas *c_RunQA_primPt = new TCanvas("c_RunQA_primPt","c_RunQA_primPt",10,10,800,400);
  c_RunQA_primPt->cd()->SetLeftMargin(0.1);
  c_RunQA_primPt->cd()->SetRightMargin(0.1);
  c_RunQA_primPt->cd()->SetBottomMargin(0.1);
  c_RunQA_primPt->cd()->SetGrid(0,0);
  c_RunQA_primPt->cd()->SetTicks(1,1);

  p_mPrimPt[1][9]->SetTitle("primPt vs. runIndex");
  p_mPrimPt[1][9]->SetStats(0);
  p_mPrimPt[1][9]->SetMarkerColor(1);
  p_mPrimPt[1][9]->SetMarkerStyle(20);
  p_mPrimPt[1][9]->SetMarkerSize(1.0);
  p_mPrimPt[1][9]->SetLineColor(1);
  p_mPrimPt[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mPrimPt[1][9]->GetYaxis()->SetTitle("<p_{T}^{prim}>");
  p_mPrimPt[1][9]->GetYaxis()->SetRangeUser(0.4,0.8);
  p_mPrimPt[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mPrimPt[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mPrimPt[1][i_trig]->SetMarkerStyle(24);
    p_mPrimPt[1][i_trig]->SetMarkerSize(0.8);
    p_mPrimPt[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mPrimPt[1][i_trig]->GetEntries() > 0) p_mPrimPt[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_primPt->SaveAs("./figures/c_RunQA_primPt.eps");
  //---------------------
  TCanvas *c_RunQA_primEta = new TCanvas("c_RunQA_primEta","c_RunQA_primEta",10,10,800,400);
  c_RunQA_primEta->cd()->SetLeftMargin(0.1);
  c_RunQA_primEta->cd()->SetRightMargin(0.1);
  c_RunQA_primEta->cd()->SetBottomMargin(0.1);
  c_RunQA_primEta->cd()->SetGrid(0,0);
  c_RunQA_primEta->cd()->SetTicks(1,1);

  p_mPrimEta[1][9]->SetTitle("primEta vs. runIndex");
  p_mPrimEta[1][9]->SetStats(0);
  p_mPrimEta[1][9]->SetMarkerColor(1);
  p_mPrimEta[1][9]->SetMarkerStyle(20);
  p_mPrimEta[1][9]->SetMarkerSize(1.0);
  p_mPrimEta[1][9]->SetLineColor(1);
  p_mPrimEta[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mPrimEta[1][9]->GetYaxis()->SetTitle("<#eta^{prim}>");
  p_mPrimEta[1][9]->GetYaxis()->SetRangeUser(-0.05,0.10);
  p_mPrimEta[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mPrimEta[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mPrimEta[1][i_trig]->SetMarkerStyle(24);
    p_mPrimEta[1][i_trig]->SetMarkerSize(0.8);
    p_mPrimEta[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mPrimEta[1][i_trig]->GetEntries() > 0) p_mPrimEta[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_primEta->SaveAs("./figures/c_RunQA_primEta.eps");
  //---------------------
  TCanvas *c_RunQA_primPhi = new TCanvas("c_RunQA_primPhi","c_RunQA_primPhi",10,10,800,400);
  c_RunQA_primPhi->cd()->SetLeftMargin(0.1);
  c_RunQA_primPhi->cd()->SetRightMargin(0.1);
  c_RunQA_primPhi->cd()->SetBottomMargin(0.1);
  c_RunQA_primPhi->cd()->SetGrid(0,0);
  c_RunQA_primPhi->cd()->SetTicks(1,1);

  p_mPrimPhi[1][9]->SetTitle("primPhi vs. runIndex");
  p_mPrimPhi[1][9]->SetStats(0);
  p_mPrimPhi[1][9]->SetMarkerColor(1);
  p_mPrimPhi[1][9]->SetMarkerStyle(20);
  p_mPrimPhi[1][9]->SetMarkerSize(1.0);
  p_mPrimPhi[1][9]->SetLineColor(1);
  p_mPrimPhi[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mPrimPhi[1][9]->GetYaxis()->SetTitle("<#phi^{prim}>");
  p_mPrimPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mPrimPhi[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mPrimPhi[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mPrimPhi[1][i_trig]->SetMarkerStyle(24);
    p_mPrimPhi[1][i_trig]->SetMarkerSize(0.8);
    p_mPrimPhi[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mPrimPhi[1][i_trig]->GetEntries() > 0) p_mPrimPhi[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_primPhi->SaveAs("./figures/c_RunQA_primPhi.eps");
  //---------------------
  TCanvas *c_RunQA_globPt = new TCanvas("c_RunQA_globPt","c_RunQA_globPt",10,10,800,400);
  c_RunQA_globPt->cd()->SetLeftMargin(0.1);
  c_RunQA_globPt->cd()->SetRightMargin(0.1);
  c_RunQA_globPt->cd()->SetBottomMargin(0.1);
  c_RunQA_globPt->cd()->SetGrid(0,0);
  c_RunQA_globPt->cd()->SetTicks(1,1);

  p_mGlobPt[1][9]->SetTitle("globPt vs. runIndex");
  p_mGlobPt[1][9]->SetStats(0);
  p_mGlobPt[1][9]->SetMarkerColor(1);
  p_mGlobPt[1][9]->SetMarkerStyle(20);
  p_mGlobPt[1][9]->SetMarkerSize(1.0);
  p_mGlobPt[1][9]->SetLineColor(1);
  p_mGlobPt[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGlobPt[1][9]->GetYaxis()->SetTitle("<p_{T}^{glob}>");
  p_mGlobPt[1][9]->GetYaxis()->SetRangeUser(0.4,0.8);
  p_mGlobPt[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGlobPt[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGlobPt[1][i_trig]->SetMarkerStyle(24);
    p_mGlobPt[1][i_trig]->SetMarkerSize(0.8);
    p_mGlobPt[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGlobPt[1][i_trig]->GetEntries() > 0) p_mGlobPt[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_globPt->SaveAs("./figures/c_RunQA_globPt.eps");
  //---------------------
  TCanvas *c_RunQA_globEta = new TCanvas("c_RunQA_globEta","c_RunQA_globEta",10,10,800,400);
  c_RunQA_globEta->cd()->SetLeftMargin(0.1);
  c_RunQA_globEta->cd()->SetRightMargin(0.1);
  c_RunQA_globEta->cd()->SetBottomMargin(0.1);
  c_RunQA_globEta->cd()->SetGrid(0,0);
  c_RunQA_globEta->cd()->SetTicks(1,1);

  p_mGlobEta[1][9]->SetTitle("globEta vs. runIndex");
  p_mGlobEta[1][9]->SetStats(0);
  p_mGlobEta[1][9]->SetMarkerColor(1);
  p_mGlobEta[1][9]->SetMarkerStyle(20);
  p_mGlobEta[1][9]->SetMarkerSize(1.0);
  p_mGlobEta[1][9]->SetLineColor(1);
  p_mGlobEta[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGlobEta[1][9]->GetYaxis()->SetTitle("<#eta^{glob}>");
  p_mGlobEta[1][9]->GetYaxis()->SetRangeUser(-0.05,0.10);
  p_mGlobEta[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGlobEta[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGlobEta[1][i_trig]->SetMarkerStyle(24);
    p_mGlobEta[1][i_trig]->SetMarkerSize(0.8);
    p_mGlobEta[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGlobEta[1][i_trig]->GetEntries() > 0) p_mGlobEta[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_globEta->SaveAs("./figures/c_RunQA_globEta.eps");
  //---------------------
  TCanvas *c_RunQA_globPhi = new TCanvas("c_RunQA_globPhi","c_RunQA_globPhi",10,10,800,400);
  c_RunQA_globPhi->cd()->SetLeftMargin(0.1);
  c_RunQA_globPhi->cd()->SetRightMargin(0.1);
  c_RunQA_globPhi->cd()->SetBottomMargin(0.1);
  c_RunQA_globPhi->cd()->SetGrid(0,0);
  c_RunQA_globPhi->cd()->SetTicks(1,1);

  p_mGlobPhi[1][9]->SetTitle("globPhi vs. runIndex");
  p_mGlobPhi[1][9]->SetStats(0);
  p_mGlobPhi[1][9]->SetMarkerColor(1);
  p_mGlobPhi[1][9]->SetMarkerStyle(20);
  p_mGlobPhi[1][9]->SetMarkerSize(1.0);
  p_mGlobPhi[1][9]->SetLineColor(1);
  p_mGlobPhi[1][9]->GetXaxis()->SetTitle("runIndex");
  p_mGlobPhi[1][9]->GetYaxis()->SetTitle("<#phi^{glob}>");
  p_mGlobPhi[1][9]->GetYaxis()->SetRangeUser(-0.1,0.50);
  p_mGlobPhi[1][9]->Draw("pE");

  for(int i_trig = 0; i_trig < numOfTriggers; ++i_trig)
  {
    p_mGlobPhi[1][i_trig]->SetMarkerColor(MarkerColor[i_trig]);
    p_mGlobPhi[1][i_trig]->SetMarkerStyle(24);
    p_mGlobPhi[1][i_trig]->SetMarkerSize(0.8);
    p_mGlobPhi[1][i_trig]->SetLineColor(MarkerColor[i_trig]);
    if(p_mGlobPhi[1][i_trig]->GetEntries() > 0) p_mGlobPhi[1][i_trig]->Draw("pE same");
  }

  leg->Draw("same");
  c_RunQA_globPhi->SaveAs("./figures/c_RunQA_globPhi.eps");
}
