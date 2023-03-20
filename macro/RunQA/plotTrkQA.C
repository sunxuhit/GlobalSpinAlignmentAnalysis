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

void plotTrkQA(int beamType = 2)
{
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/RunQA/%s/file_%s_RunQA.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mPrimPt[numCuts][numTriggerBins]; // Track Kinematics
  TH1F *h_mPrimEta[numCuts][numTriggerBins];
  TH1F *h_mPrimPhi[numCuts][numTriggerBins];
  TH1F *h_mGlobPt[numCuts][numTriggerBins];
  TH1F *h_mGlobEta[numCuts][numTriggerBins];
  TH1F *h_mGlobPhi[numCuts][numTriggerBins];

  TH1F *h_mDca[numCuts][numTriggerBins]; // Track Quality
  TH1F *h_mNHitsFit[numCuts][numTriggerBins];
  TH1F *h_mNHitsRatio[numCuts][numTriggerBins];
  TH1F *h_mNHitsDEdx[numCuts][numTriggerBins];

  TH2F *h_mDEdxMom[numCuts][numTriggerBins]; // Track PID
  TH2F *h_mMass2Mom[numCuts][numTriggerBins];
  TH2F *h_mBetaMom[numCuts][numTriggerBins];

  TH1F *h_mPrimEtaEpFull[numCuts][numTriggerBins]; // test StAnalysisCuts
  TH1F *h_mPrimEtaEpEast[numCuts][numTriggerBins];
  TH1F *h_mPrimEtaEpWest[numCuts][numTriggerBins];
  TH1F *h_mPrimEtaFlowFull[numCuts][numTriggerBins];
  TH1F *h_mPrimEtaFlowEast[numCuts][numTriggerBins];
  TH1F *h_mPrimEtaFlowWest[numCuts][numTriggerBins];
  TH2F *h_mPrimEtaNSigKaonFull[numCuts][numTriggerBins];
  TH2F *h_mPrimEtaNSigKaonEast[numCuts][numTriggerBins];
  TH2F *h_mPrimEtaNSigKaonWest[numCuts][numTriggerBins];

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < numTriggerBins; ++iTrig)
    {
      std::string histName = Form("h_mPrimPt%sTrigger%d",CutStatus[iCut].c_str(),iTrig); // Track Kinematics
      h_mPrimPt[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEta%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEta[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimPhi%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimPhi[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mGlobPt%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mGlobPt[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mGlobEta%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mGlobEta[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mGlobPhi%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mGlobPhi[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mDca%sTrigger%d",CutStatus[iCut].c_str(),iTrig); // Track Quality
      h_mDca[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mNHitsFit%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mNHitsFit[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mNHitsRatio%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mNHitsRatio[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mNHitsDEdx%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mNHitsDEdx[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mDEdxMom%sTrigger%d",CutStatus[iCut].c_str(),iTrig); // Track PID
      h_mDEdxMom[iCut][iTrig] = (TH2F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mMass2Mom%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mMass2Mom[iCut][iTrig] = (TH2F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mBetaMom%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mBetaMom[iCut][iTrig] = (TH2F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEtaEpFull%sTrigger%d",CutStatus[iCut].c_str(),iTrig); // test StAnalysisCuts
      h_mPrimEtaEpFull[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEtaEpEast%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEtaEpEast[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEtaEpWest%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEtaEpWest[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEtaFlowFull%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEtaFlowFull[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEtaFlowEast%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEtaFlowEast[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEtaFlowWest%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEtaFlowWest[iCut][iTrig] = (TH1F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEtaNSigKaonFull%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEtaNSigKaonFull[iCut][iTrig] = (TH2F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEtaNSigKaonEast%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEtaNSigKaonEast[iCut][iTrig] = (TH2F*)File_InPut->Get(histName.c_str());

      histName = Form("h_mPrimEtaNSigKaonWest%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEtaNSigKaonWest[iCut][iTrig] = (TH2F*)File_InPut->Get(histName.c_str());
    }
  }

  TCanvas *c_TrackQA = new TCanvas("c_TrackQA","c_TrackQA",10,10,1200,800);
  c_TrackQA->Divide(3,2);
  for(int i_pad = 0; i_pad < 6; ++i_pad)
  {
    c_TrackQA->cd(i_pad+1);
    c_TrackQA->cd(i_pad+1)->SetLeftMargin(0.1);
    c_TrackQA->cd(i_pad+1)->SetRightMargin(0.1);
    c_TrackQA->cd(i_pad+1)->SetBottomMargin(0.1);
    c_TrackQA->cd(i_pad+1)->SetGrid(0,0);
    c_TrackQA->cd(i_pad+1)->SetTicks(1,1);
  }
  std::string figName = Form("../../figures/RunQA/%s/TrkQA_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TrackQA->Print(figName.c_str());

  // primMom
  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_TrackQA->cd(iCut*3+1);
    c_TrackQA->cd(iCut*3+1)->Clear();
    c_TrackQA->cd(iCut*3+1)->SetLogy();
    h_mPrimPt[iCut][9]->GetXaxis()->SetTitle("prim p_{T}");
    h_mPrimPt[iCut][9]->Draw("hE");

    c_TrackQA->cd(iCut*3+2);
    c_TrackQA->cd(iCut*3+2)->Clear();
    h_mPrimEta[iCut][9]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEta[iCut][9]->Draw("hE");

    c_TrackQA->cd(iCut*3+3);
    c_TrackQA->cd(iCut*3+3)->Clear();;
    h_mPrimPhi[iCut][9]->GetXaxis()->SetTitle("prim #phi");
    h_mPrimPhi[iCut][9]->Draw("hE");
  }
  figName = Form("../../figures/RunQA/%s/TrkQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TrackQA->Update();
  c_TrackQA->Print(figName.c_str());

  // globMom
  for(int iCut = 0; iCut < 2; ++iCut)
  {
    c_TrackQA->cd(iCut*3+1);
    c_TrackQA->cd(iCut*3+1)->Clear();
    c_TrackQA->cd(iCut*3+1)->SetLogy();
    h_mGlobPt[iCut][9]->GetXaxis()->SetTitle("glob p_{T}");
    h_mGlobPt[iCut][9]->Draw("hE");

    c_TrackQA->cd(iCut*3+2);
    c_TrackQA->cd(iCut*3+2)->Clear();
    h_mGlobEta[iCut][9]->GetXaxis()->SetTitle("glob #eta");
    h_mGlobEta[iCut][9]->Draw("hE");

    c_TrackQA->cd(iCut*3+3);
    c_TrackQA->cd(iCut*3+3)->Clear();
    h_mGlobPhi[iCut][9]->GetXaxis()->SetTitle("glob #phi");
    h_mGlobPhi[iCut][9]->Draw("hE");
  }
  figName = Form("../../figures/RunQA/%s/TrkQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TrackQA->Update();
  c_TrackQA->Print(figName.c_str());

  // Track Quality
  TCanvas *c_TrackQA_Quality = new TCanvas("c_TrackQA_Quality","c_TrackQA_Quality",10,10,1600,800);
  c_TrackQA_Quality->Divide(4,2);
  for(int i_pad = 0; i_pad < 8; ++i_pad)
  {
    c_TrackQA_Quality->cd(i_pad+1);
    c_TrackQA_Quality->cd(i_pad+1)->SetLeftMargin(0.1);
    c_TrackQA_Quality->cd(i_pad+1)->SetRightMargin(0.1);
    c_TrackQA_Quality->cd(i_pad+1)->SetBottomMargin(0.1);
    c_TrackQA_Quality->cd(i_pad+1)->SetGrid(0,0);
    c_TrackQA_Quality->cd(i_pad+1)->SetTicks(1,1);
  }

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_TrackQA_Quality->cd(iCut*4+1)->SetLogy();
    h_mDca[iCut][9]->GetXaxis()->SetTitle("dca");
    h_mDca[iCut][9]->Draw("hE");

    c_TrackQA_Quality->cd(iCut*4+2);
    h_mNHitsFit[iCut][9]->Draw("hE");
    h_mNHitsFit[iCut][9]->GetXaxis()->SetTitle("nHitsFit");

    c_TrackQA_Quality->cd(iCut*4+3);
    h_mNHitsRatio[iCut][9]->Draw("hE");
    h_mNHitsRatio[iCut][9]->GetXaxis()->SetTitle("nHitsRatio");

    c_TrackQA_Quality->cd(iCut*4+4);
    h_mNHitsDEdx[iCut][9]->Draw("hE");
    h_mNHitsDEdx[iCut][9]->GetXaxis()->SetTitle("nHitsDedx");
  }
  figName = Form("../../figures/RunQA/%s/TrkQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TrackQA_Quality->Update();
  c_TrackQA_Quality->Print(figName.c_str());

  // PID
  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_TrackQA->cd(iCut*3+1);
    c_TrackQA->cd(iCut*3+1)->Clear();
    c_TrackQA->cd(iCut*3+1)->SetLogy(0);
    c_TrackQA->cd(iCut*3+1)->SetLogz();
    h_mDEdxMom[iCut][9]->GetXaxis()->SetTitle("p*q");
    h_mDEdxMom[iCut][9]->GetYaxis()->SetTitle("dE/dx");
    h_mDEdxMom[iCut][9]->Draw("colz");

    c_TrackQA->cd(iCut*3+2);
    c_TrackQA->cd(iCut*3+2)->Clear();
    c_TrackQA->cd(iCut*3+2)->SetLogy(0);
    c_TrackQA->cd(iCut*3+2)->SetLogz();
    h_mMass2Mom[iCut][9]->GetXaxis()->SetTitle("p*q");
    h_mMass2Mom[iCut][9]->GetYaxis()->SetTitle("m^{2}");
    h_mMass2Mom[iCut][9]->Draw("colz");

    c_TrackQA->cd(iCut*3+3);
    c_TrackQA->cd(iCut*3+3)->Clear();
    c_TrackQA->cd(iCut*3+3)->SetLogy(0);
    c_TrackQA->cd(iCut*3+3)->SetLogz();
    h_mBetaMom[iCut][9]->GetXaxis()->SetTitle("p*q");
    h_mBetaMom[iCut][9]->GetYaxis()->SetTitle("1/#beta");
    h_mBetaMom[iCut][9]->Draw("colz");
  }
  figName = Form("../../figures/RunQA/%s/TrkQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TrackQA->Update();
  c_TrackQA->Print(figName.c_str());

  // Eta EP Cut
  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_TrackQA->cd(iCut*3+1);
    c_TrackQA->cd(iCut*3+1)->Clear();
    c_TrackQA->cd(iCut*3+1)->SetLogy();
    h_mPrimEtaEpFull[iCut][0]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEtaEpFull[iCut][0]->Draw("hE");

    c_TrackQA->cd(iCut*3+2);
    c_TrackQA->cd(iCut*3+2)->Clear();
    h_mPrimEtaEpEast[iCut][0]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEtaEpEast[iCut][0]->Draw("hE");

    c_TrackQA->cd(iCut*3+3);
    c_TrackQA->cd(iCut*3+3)->Clear();;
    h_mPrimEtaEpWest[iCut][0]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEtaEpWest[iCut][0]->Draw("hE");
  }
  figName = Form("../../figures/RunQA/%s/TrkQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TrackQA->Update();
  c_TrackQA->Print(figName.c_str());

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_TrackQA->cd(iCut*3+1);
    c_TrackQA->cd(iCut*3+1)->Clear();
    c_TrackQA->cd(iCut*3+1)->SetLogy();
    h_mPrimEtaFlowFull[iCut][0]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEtaFlowFull[iCut][0]->Draw("hE");

    c_TrackQA->cd(iCut*3+2);
    c_TrackQA->cd(iCut*3+2)->Clear();
    h_mPrimEtaFlowEast[iCut][0]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEtaFlowEast[iCut][0]->Draw("hE");

    c_TrackQA->cd(iCut*3+3);
    c_TrackQA->cd(iCut*3+3)->Clear();;
    h_mPrimEtaFlowWest[iCut][0]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEtaFlowWest[iCut][0]->Draw("hE");
  }
  figName = Form("../../figures/RunQA/%s/TrkQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TrackQA->Update();
  c_TrackQA->Print(figName.c_str());

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_TrackQA->cd(iCut*3+1);
    c_TrackQA->cd(iCut*3+1)->Clear();
    c_TrackQA->cd(iCut*3+1)->SetLogy(0);;
    h_mPrimEtaNSigKaonFull[iCut][0]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEtaNSigKaonFull[iCut][0]->GetYaxis()->SetTitle("n#sigma_{K}");
    h_mPrimEtaNSigKaonFull[iCut][0]->Draw("colz");

    c_TrackQA->cd(iCut*3+2);
    c_TrackQA->cd(iCut*3+2)->Clear();
    c_TrackQA->cd(iCut*3+2)->SetLogy(0);;
    h_mPrimEtaNSigKaonEast[iCut][0]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEtaNSigKaonEast[iCut][0]->GetYaxis()->SetTitle("n#sigma_{K}");
    h_mPrimEtaNSigKaonEast[iCut][0]->Draw("colz");

    c_TrackQA->cd(iCut*3+3);
    c_TrackQA->cd(iCut*3+3)->Clear();;
    c_TrackQA->cd(iCut*3+3)->SetLogy(0);;
    h_mPrimEtaNSigKaonWest[iCut][0]->GetXaxis()->SetTitle("prim #eta");
    h_mPrimEtaNSigKaonWest[iCut][0]->GetYaxis()->SetTitle("n#sigma_{K}");
    h_mPrimEtaNSigKaonWest[iCut][0]->Draw("colz");
  }
  figName = Form("../../figures/RunQA/%s/TrkQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TrackQA->Update();
  c_TrackQA->Print(figName.c_str());

  figName = Form("../../figures/RunQA/%s/TrkQA_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_TrackQA->Print(figName.c_str());
}
