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

void plotRunQA(int beamType = 2)
{
  vector<string> vecTriggerID;
  vector<int> vecMarkerColor;
  vecTriggerID.clear();
  vecMarkerColor.clear();
  if(beamType == 0 || beamType == 1) 
  {
    vecTriggerID.push_back("600001");
    vecTriggerID.push_back("600011");
    vecTriggerID.push_back("600021");
    vecTriggerID.push_back("600031");
    vecMarkerColor.push_back(7);
    vecMarkerColor.push_back(6);
    vecMarkerColor.push_back(4);
    vecMarkerColor.push_back(2);
  }
  if(beamType == 2) 
  {
    vecTriggerID.push_back("620052");
    vecMarkerColor.push_back(4);
  }
  const int numOfTriggers = (int)vecTriggerID.size();

  string inputfile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/RunQA/%s/file_%s_RunQA.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TH1F *h_mTriggerId[2]; // 0: before cuts | 1: after cuts

  for(int iCut = 0; iCut < 2; ++iCut)
  {
    string histName = Form("h_mTriggerId%s",CutStatus[iCut].c_str());
    h_mTriggerId[iCut] = (TH1F*)File_InPut->Get(histName.c_str());
    h_mTriggerId[iCut]->SetLineColor(iCut+1);
    h_mTriggerId[iCut]->GetXaxis()->SetTitle("triggerId");
    for(int iTrigger = 0; iTrigger < vecTriggerID.size(); ++iTrigger)
    {
      h_mTriggerId[iCut]->GetXaxis()->SetBinLabel(iTrigger+1,vecTriggerID[iTrigger].c_str());
    }
  }

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



  TCanvas *c_RunQA = new TCanvas("c_RunQA","c_RunQA",10,10,800,400);
  c_RunQA->cd()->SetLeftMargin(0.15);
  c_RunQA->cd()->SetRightMargin(0.15);
  c_RunQA->cd()->SetBottomMargin(0.15);
  c_RunQA->cd()->SetTicks(1,1);
  c_RunQA->cd()->SetGrid(0,0);

  std::string figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_RunQA->Print(figName.c_str());

  TLegend *leg; 
  if(beamType == 0 || beamType == 1) leg = new TLegend(0.7,0.6,0.9,0.9);
  if(beamType == 2) leg = new TLegend(0.7,0.75,0.9,0.9);

  { // refMult
    c_RunQA->Clear();
    p_mRefMult[1][9]->SetTitle("refMult vs. runIndex");
    p_mRefMult[1][9]->SetStats(0);
    p_mRefMult[1][9]->SetMarkerColor(1);
    p_mRefMult[1][9]->SetMarkerStyle(20);
    p_mRefMult[1][9]->SetMarkerSize(1.0);
    p_mRefMult[1][9]->SetLineColor(1);
    p_mRefMult[1][9]->GetXaxis()->SetTitle("runIndex");
    p_mRefMult[1][9]->GetYaxis()->SetTitle("<refMult>");
    p_mRefMult[1][9]->GetYaxis()->SetRangeUser(0,400);
    if(beamType == 2) p_mRefMult[1][9]->GetYaxis()->SetRangeUser(0,200);
    p_mRefMult[1][9]->Draw("pE");
    leg->AddEntry(p_mRefMult[1][9],"All Triggers","P");

    for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
    {
      p_mRefMult[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mRefMult[1][iTrig]->SetMarkerStyle(24);
      p_mRefMult[1][iTrig]->SetMarkerSize(0.8);
      p_mRefMult[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mRefMult[1][iTrig]->GetEntries() > 0) p_mRefMult[1][iTrig]->Draw("pE same");
      leg->AddEntry(p_mRefMult[1][iTrig],vecTriggerID[iTrig].c_str(),"P");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    // global refMult
    c_RunQA->Clear();
    p_mGRefMult[1][9]->SetTitle("grefMult vs. runIndex");
    p_mGRefMult[1][9]->SetStats(0);
    p_mGRefMult[1][9]->SetMarkerColor(1);
    p_mGRefMult[1][9]->SetMarkerStyle(20);
    p_mGRefMult[1][9]->SetMarkerSize(1.0);
    p_mGRefMult[1][9]->SetLineColor(1);
    p_mGRefMult[1][9]->GetXaxis()->SetTitle("runIndex");
    p_mGRefMult[1][9]->GetYaxis()->SetTitle("<grefMult>");
    p_mGRefMult[1][9]->GetYaxis()->SetRangeUser(0,400);
    if(beamType == 2) p_mGRefMult[1][9]->GetYaxis()->SetRangeUser(0,200);
    p_mGRefMult[1][9]->Draw("pE");

    for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
    {
      p_mGRefMult[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mGRefMult[1][iTrig]->SetMarkerStyle(24);
      p_mGRefMult[1][iTrig]->SetMarkerSize(0.8);
      p_mGRefMult[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mGRefMult[1][iTrig]->GetEntries() > 0) p_mGRefMult[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  if(beamType == 0 || beamType == 1)
  {
    c_RunQA->Clear();
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
      p_mZdcX[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mZdcX[1][iTrig]->SetMarkerStyle(24);
      p_mZdcX[1][iTrig]->SetMarkerSize(0.8);
      p_mZdcX[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mZdcX[1][iTrig]->GetEntries() > 0) p_mZdcX[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    c_RunQA->Clear();
    p_mVz[1][9]->SetTitle("Vz vs. runIndex");
    p_mVz[1][9]->SetStats(0);
    p_mVz[1][9]->SetMarkerColor(1);
    p_mVz[1][9]->SetMarkerStyle(20);
    p_mVz[1][9]->SetMarkerSize(1.0);
    p_mVz[1][9]->SetLineColor(1);
    p_mVz[1][9]->GetXaxis()->SetTitle("runIndex");
    p_mVz[1][9]->GetYaxis()->SetTitle("<Vz>");
    p_mVz[1][9]->GetYaxis()->SetRangeUser(-3.0,3.0);
    if(beamType == 2) p_mVz[1][9]->GetYaxis()->SetRangeUser(195,205);
    p_mVz[1][9]->Draw("pE");

    for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
    {
      p_mVz[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mVz[1][iTrig]->SetMarkerStyle(24);
      p_mVz[1][iTrig]->SetMarkerSize(0.8);
      p_mVz[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mVz[1][iTrig]->GetEntries() > 0) p_mVz[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    c_RunQA->Clear();
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
      p_mVr[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mVr[1][iTrig]->SetMarkerStyle(24);
      p_mVr[1][iTrig]->SetMarkerSize(0.8);
      p_mVr[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mVr[1][iTrig]->GetEntries() > 0) p_mVr[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    c_RunQA->Clear();
    p_mGDca[1][9]->SetTitle("gDca vs. runIndex");
    p_mGDca[1][9]->SetStats(0);
    p_mGDca[1][9]->SetMarkerColor(1);
    p_mGDca[1][9]->SetMarkerStyle(20);
    p_mGDca[1][9]->SetMarkerSize(1.0);
    p_mGDca[1][9]->SetLineColor(1);
    p_mGDca[1][9]->GetXaxis()->SetTitle("runIndex");
    p_mGDca[1][9]->GetYaxis()->SetTitle("<gDca>");
    p_mGDca[1][9]->GetYaxis()->SetRangeUser(0.2,1.0);
    p_mGDca[1][9]->Draw("pE");

    for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
    {
      p_mGDca[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mGDca[1][iTrig]->SetMarkerStyle(24);
      p_mGDca[1][iTrig]->SetMarkerSize(0.8);
      p_mGDca[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mGDca[1][iTrig]->GetEntries() > 0) p_mGDca[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    c_RunQA->Clear();
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
      p_mNHitsFit[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mNHitsFit[1][iTrig]->SetMarkerStyle(24);
      p_mNHitsFit[1][iTrig]->SetMarkerSize(0.8);
      p_mNHitsFit[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mNHitsFit[1][iTrig]->GetEntries() > 0) p_mNHitsFit[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  //---------------------
  {
    c_RunQA->Clear();
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
      p_mPrimPt[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mPrimPt[1][iTrig]->SetMarkerStyle(24);
      p_mPrimPt[1][iTrig]->SetMarkerSize(0.8);
      p_mPrimPt[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mPrimPt[1][iTrig]->GetEntries() > 0) p_mPrimPt[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    c_RunQA->Clear();
    p_mPrimEta[1][9]->SetTitle("primEta vs. runIndex");
    p_mPrimEta[1][9]->SetStats(0);
    p_mPrimEta[1][9]->SetMarkerColor(1);
    p_mPrimEta[1][9]->SetMarkerStyle(20);
    p_mPrimEta[1][9]->SetMarkerSize(1.0);
    p_mPrimEta[1][9]->SetLineColor(1);
    p_mPrimEta[1][9]->GetXaxis()->SetTitle("runIndex");
    p_mPrimEta[1][9]->GetYaxis()->SetTitle("<#eta^{prim}>");
    p_mPrimEta[1][9]->GetYaxis()->SetRangeUser(-0.05,0.10);
    if(beamType == 2) p_mPrimEta[1][9]->GetYaxis()->SetRangeUser(-2.0,0.0);
    p_mPrimEta[1][9]->Draw("pE");

    for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
    {
      p_mPrimEta[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mPrimEta[1][iTrig]->SetMarkerStyle(24);
      p_mPrimEta[1][iTrig]->SetMarkerSize(0.8);
      p_mPrimEta[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mPrimEta[1][iTrig]->GetEntries() > 0) p_mPrimEta[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    c_RunQA->Clear();
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
      p_mPrimPhi[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mPrimPhi[1][iTrig]->SetMarkerStyle(24);
      p_mPrimPhi[1][iTrig]->SetMarkerSize(0.8);
      p_mPrimPhi[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mPrimPhi[1][iTrig]->GetEntries() > 0) p_mPrimPhi[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    c_RunQA->Clear();
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
      p_mGlobPt[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mGlobPt[1][iTrig]->SetMarkerStyle(24);
      p_mGlobPt[1][iTrig]->SetMarkerSize(0.8);
      p_mGlobPt[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mGlobPt[1][iTrig]->GetEntries() > 0) p_mGlobPt[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    c_RunQA->Clear();
    p_mGlobEta[1][9]->SetTitle("globEta vs. runIndex");
    p_mGlobEta[1][9]->SetStats(0);
    p_mGlobEta[1][9]->SetMarkerColor(1);
    p_mGlobEta[1][9]->SetMarkerStyle(20);
    p_mGlobEta[1][9]->SetMarkerSize(1.0);
    p_mGlobEta[1][9]->SetLineColor(1);
    p_mGlobEta[1][9]->GetXaxis()->SetTitle("runIndex");
    p_mGlobEta[1][9]->GetYaxis()->SetTitle("<#eta^{glob}>");
    p_mGlobEta[1][9]->GetYaxis()->SetRangeUser(-0.05,0.10);
    if(beamType == 2) p_mGlobEta[1][9]->GetYaxis()->SetRangeUser(-2.0,0.0);
    p_mGlobEta[1][9]->Draw("pE");

    for(int iTrig = 0; iTrig < numOfTriggers; ++iTrig)
    {
      p_mGlobEta[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mGlobEta[1][iTrig]->SetMarkerStyle(24);
      p_mGlobEta[1][iTrig]->SetMarkerSize(0.8);
      p_mGlobEta[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mGlobEta[1][iTrig]->GetEntries() > 0) p_mGlobEta[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
  }

  {
    c_RunQA->Clear();
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
      p_mGlobPhi[1][iTrig]->SetMarkerColor(vecMarkerColor[iTrig]);
      p_mGlobPhi[1][iTrig]->SetMarkerStyle(24);
      p_mGlobPhi[1][iTrig]->SetMarkerSize(0.8);
      p_mGlobPhi[1][iTrig]->SetLineColor(vecMarkerColor[iTrig]);
      if(p_mGlobPhi[1][iTrig]->GetEntries() > 0) p_mGlobPhi[1][iTrig]->Draw("pE same");
    }
    leg->Draw("same");

    figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_RunQA->Update();
    c_RunQA->Print(figName.c_str());
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
  figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TriggerId->Update();
  c_TriggerId->Print(figName.c_str());

  figName = Form("../../figures/RunQA/%s/RunQA_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_RunQA->Print(figName.c_str());
}
