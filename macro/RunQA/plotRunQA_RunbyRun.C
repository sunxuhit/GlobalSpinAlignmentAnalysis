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

void plotRunQA_RunbyRun(int beamType = 0)
{
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/%s/RunQA/file_%s_RunQA.root",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TProfile *p_mRefMult[numCuts][numTriggerBins];
  TProfile *p_mGRefMult[numCuts][numTriggerBins];
  TProfile *p_mZdcX[numCuts][numTriggerBins];
  TProfile *p_mVz[numCuts][numTriggerBins];
  TProfile *p_mVr[numCuts][numTriggerBins];

  TProfile *p_mGDca[numCuts][numTriggerBins];
  TProfile *p_mNHitsFit[numCuts][numTriggerBins];
  TProfile *p_mPrimPt[numCuts][numTriggerBins];
  TProfile *p_mPrimEta[numCuts][numTriggerBins];
  TProfile *p_mPrimPhi[numCuts][numTriggerBins];
  TProfile *p_mGlobPt[numCuts][numTriggerBins];
  TProfile *p_mGlobEta[numCuts][numTriggerBins];
  TProfile *p_mGlobPhi[numCuts][numTriggerBins];

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < numTriggerBins; ++iTrig)
    {
      std::string ProName = Form("p_mRefMult%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mRefMult[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGRefMult%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mGRefMult[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mZdcX%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mZdcX[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mVz%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mVz[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mVr%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mVr[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGDca%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mGDca[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mNHitsFit%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mNHitsFit[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimPt%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mPrimPt[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimEta%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mPrimEta[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mPrimPhi%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mPrimPhi[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGlobPt%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mGlobPt[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGlobEta%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mGlobEta[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());

      ProName = Form("p_mGlobPhi%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      p_mGlobPhi[iCut][iTrig] = (TProfile*)File_InPut->Get(ProName.c_str());
    }
  }

  const int numOfTriggers = 4;
  const int triggerID[numOfTriggers] = {600001,600011,600021,600031};
  const int MarkerColor[numOfTriggers] = {7,6,4,2};

  TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);
  string FigName;

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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mRefMult[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mRefMult[1][iTrig]->SetMarkerStyle(24);
    p_mRefMult[1][iTrig]->SetMarkerSize(0.8);
    p_mRefMult[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mRefMult[1][iTrig]->GetEntries() > 0) p_mRefMult[1][iTrig]->Draw("pE same");
    leg->AddEntry(p_mRefMult[1][iTrig],to_string(triggerID[iTrig]).c_str(),"P");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_RefMult_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_RefMult->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mGRefMult[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mGRefMult[1][iTrig]->SetMarkerStyle(24);
    p_mGRefMult[1][iTrig]->SetMarkerSize(0.8);
    p_mGRefMult[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mGRefMult[1][iTrig]->GetEntries() > 0) p_mGRefMult[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_gRefMult_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_gRefMult->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mZdcX[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mZdcX[1][iTrig]->SetMarkerStyle(24);
    p_mZdcX[1][iTrig]->SetMarkerSize(0.8);
    p_mZdcX[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mZdcX[1][iTrig]->GetEntries() > 0) p_mZdcX[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_ZdcX_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_ZdcX->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mVz[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mVz[1][iTrig]->SetMarkerStyle(24);
    p_mVz[1][iTrig]->SetMarkerSize(0.8);
    p_mVz[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mVz[1][iTrig]->GetEntries() > 0) p_mVz[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_Vz_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_Vz->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mVr[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mVr[1][iTrig]->SetMarkerStyle(24);
    p_mVr[1][iTrig]->SetMarkerSize(0.8);
    p_mVr[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mVr[1][iTrig]->GetEntries() > 0) p_mVr[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_Vr_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_Vr->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mGDca[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mGDca[1][iTrig]->SetMarkerStyle(24);
    p_mGDca[1][iTrig]->SetMarkerSize(0.8);
    p_mGDca[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mGDca[1][iTrig]->GetEntries() > 0) p_mGDca[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_gDca_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_gDca->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mNHitsFit[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mNHitsFit[1][iTrig]->SetMarkerStyle(24);
    p_mNHitsFit[1][iTrig]->SetMarkerSize(0.8);
    p_mNHitsFit[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mNHitsFit[1][iTrig]->GetEntries() > 0) p_mNHitsFit[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_nHitsFit_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_nHitsFit->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mPrimPt[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mPrimPt[1][iTrig]->SetMarkerStyle(24);
    p_mPrimPt[1][iTrig]->SetMarkerSize(0.8);
    p_mPrimPt[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mPrimPt[1][iTrig]->GetEntries() > 0) p_mPrimPt[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_primPt_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_primPt->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mPrimEta[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mPrimEta[1][iTrig]->SetMarkerStyle(24);
    p_mPrimEta[1][iTrig]->SetMarkerSize(0.8);
    p_mPrimEta[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mPrimEta[1][iTrig]->GetEntries() > 0) p_mPrimEta[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_primEta_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_primEta->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mPrimPhi[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mPrimPhi[1][iTrig]->SetMarkerStyle(24);
    p_mPrimPhi[1][iTrig]->SetMarkerSize(0.8);
    p_mPrimPhi[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mPrimPhi[1][iTrig]->GetEntries() > 0) p_mPrimPhi[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_primPhi_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_primPhi->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mGlobPt[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mGlobPt[1][iTrig]->SetMarkerStyle(24);
    p_mGlobPt[1][iTrig]->SetMarkerSize(0.8);
    p_mGlobPt[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mGlobPt[1][iTrig]->GetEntries() > 0) p_mGlobPt[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_globPt_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_globPt->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mGlobEta[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mGlobEta[1][iTrig]->SetMarkerStyle(24);
    p_mGlobEta[1][iTrig]->SetMarkerSize(0.8);
    p_mGlobEta[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mGlobEta[1][iTrig]->GetEntries() > 0) p_mGlobEta[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_globEta_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_globEta->SaveAs(FigName.c_str());
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

  for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
  {
    p_mGlobPhi[1][iTrig]->SetMarkerColor(MarkerColor[iTrig]);
    p_mGlobPhi[1][iTrig]->SetMarkerStyle(24);
    p_mGlobPhi[1][iTrig]->SetMarkerSize(0.8);
    p_mGlobPhi[1][iTrig]->SetLineColor(MarkerColor[iTrig]);
    if(p_mGlobPhi[1][iTrig]->GetEntries() > 0) p_mGlobPhi[1][iTrig]->Draw("pE same");
  }

  leg->Draw("same");
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_RunQA_globPhi_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_RunQA_globPhi->SaveAs(FigName.c_str());
}
