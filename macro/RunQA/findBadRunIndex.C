#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"

#include "../../Utility/include/StSpinAlignmentCons.h"
#include "../../Utility/include/draw.h"

using namespace std;

static const string CutStatus[2] = {"Bf","Af"};
static const int numCuts = 2; // 0: before cuts | 1: after cuts
static const int numTriggerBins = 5; // 0-3 for different triggerID | 4 for all triggers

void findMean(TProfile *p_runQA, double &mean, double &sigma);
bool isBadRun(double val, double mean, double sigma);
void plotGoodRunRange(double runIndexStart, double runIndexStop, double mean, double sigma);

int findBadRunIndex(int beamType = 2)
{
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/RunQA/%s/file_RunQA_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
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


  string outputfile = Form("../../Utility/RunIndex/%s/badRunIndexUnSorted_%s.txt",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  ofstream file_badRunIndex;
  file_badRunIndex.open(outputfile.c_str());
  if (!file_badRunIndex.is_open()) 
  {
    std::cout << "failed to open " << outputfile << '\n';
    return -1;
  } 

  const double runIndexStart = -0.5;
  const double runIndexStop  = 299.5;

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
  const int numOfUsedTrigs = (int)vecTriggerID.size();

  TCanvas *c_BadRun = new TCanvas("c_BadRun","c_BadRun",10,10,800,400);
  c_BadRun->cd()->SetLeftMargin(0.15);
  c_BadRun->cd()->SetRightMargin(0.15);
  c_BadRun->cd()->SetBottomMargin(0.15);
  c_BadRun->cd()->SetTicks(1,1);
  c_BadRun->cd()->SetGrid(0,0);

  std::string figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Print(figName.c_str());


  //-------------refMult----------------
  double meanRefMult = 0.0;
  double sigmaRefMult = 0.0;
  findMean(p_mRefMult[1][numTriggerBins-1],meanRefMult,sigmaRefMult);

  // TCanvas *c_RunQA_RefMult = new TCanvas("c_RunQA_RefMult","c_RunQA_RefMult",10,10,800,400);
  // c_RunQA_RefMult->cd()->SetLeftMargin(0.1);
  // c_RunQA_RefMult->cd()->SetRightMargin(0.1);
  // c_RunQA_RefMult->cd()->SetBottomMargin(0.1);
  // c_RunQA_RefMult->cd()->SetGrid(0,0);
  // c_RunQA_RefMult->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mRefMult[1][numTriggerBins-1]->SetTitle("refMult vs. runIndex");
  p_mRefMult[1][numTriggerBins-1]->SetStats(0);
  p_mRefMult[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mRefMult[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mRefMult[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mRefMult[1][numTriggerBins-1]->SetLineColor(1);
  p_mRefMult[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mRefMult[1][numTriggerBins-1]->GetYaxis()->SetTitle("<refMult>");
  p_mRefMult[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(0,400);
  if(beamType == 2) p_mRefMult[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(40,80);
  p_mRefMult[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanRefMult,sigmaRefMult);

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_RefMult_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_RefMult->SaveAs(FigName.c_str());

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  for(int i_run = 0; i_run < p_mRefMult[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mRefMult[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mRefMult[1][numTriggerBins-1]->GetBinContent(i_run+1),meanRefMult,sigmaRefMult))
    {
      cout << "bad runIndex from refMult: " << p_mRefMult[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mRefMult[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------refMult----------------

  //-------------grefMult----------------
  double meanGRefMult = 0.0;
  double sigmaGRefMult = 0.0;
  findMean(p_mGRefMult[1][numTriggerBins-1],meanGRefMult,sigmaGRefMult);

  // TCanvas *c_RunQA_gRefMult = new TCanvas("c_RunQA_gRefMult","c_RunQA_gRefMult",10,10,800,400);
  // c_RunQA_gRefMult->cd()->SetLeftMargin(0.1);
  // c_RunQA_gRefMult->cd()->SetRightMargin(0.1);
  // c_RunQA_gRefMult->cd()->SetBottomMargin(0.1);
  // c_RunQA_gRefMult->cd()->SetGrid(0,0);
  // c_RunQA_gRefMult->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mGRefMult[1][numTriggerBins-1]->SetTitle("grefMult vs. runIndex");
  p_mGRefMult[1][numTriggerBins-1]->SetStats(0);
  p_mGRefMult[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mGRefMult[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mGRefMult[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mGRefMult[1][numTriggerBins-1]->SetLineColor(1);
  p_mGRefMult[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mGRefMult[1][numTriggerBins-1]->GetYaxis()->SetTitle("<grefMult>");
  p_mGRefMult[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(0,400);
  if(beamType == 2) p_mGRefMult[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(0,200);
  p_mGRefMult[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanGRefMult,sigmaGRefMult);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_gRefMult_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_gRefMult->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mGRefMult[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mGRefMult[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mGRefMult[1][numTriggerBins-1]->GetBinContent(i_run+1),meanGRefMult,sigmaGRefMult))
    {
      cout << "bad runIndex from grefMult: " << p_mGRefMult[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGRefMult[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------grefMult----------------

  //-------------ZdcX----------------
  if(beamType !=2)
  {
    double meanZdcX = 0.0;
    double sigmaZdcX = 0.0;
    findMean(p_mZdcX[1][numTriggerBins-1],meanZdcX,sigmaZdcX);

    // TCanvas *c_RunQA_ZdcX = new TCanvas("c_RunQA_ZdcX","c_RunQA_ZdcX",10,10,800,400);
    // c_RunQA_ZdcX->cd()->SetLeftMargin(0.1);
    // c_RunQA_ZdcX->cd()->SetRightMargin(0.1);
    // c_RunQA_ZdcX->cd()->SetBottomMargin(0.1);
    // c_RunQA_ZdcX->cd()->SetGrid(0,0);
    // c_RunQA_ZdcX->cd()->SetTicks(1,1);

    c_BadRun->Clear();
    p_mZdcX[1][numTriggerBins-1]->SetTitle("ZdcX vs. runIndex");
    p_mZdcX[1][numTriggerBins-1]->SetStats(0);
    p_mZdcX[1][numTriggerBins-1]->SetMarkerColor(1);
    p_mZdcX[1][numTriggerBins-1]->SetMarkerStyle(20);
    p_mZdcX[1][numTriggerBins-1]->SetMarkerSize(1.0);
    p_mZdcX[1][numTriggerBins-1]->SetLineColor(1);
    p_mZdcX[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
    p_mZdcX[1][numTriggerBins-1]->GetYaxis()->SetTitle("<ZdcX>");
    // p_mZdcX[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(0,60000);
    p_mZdcX[1][numTriggerBins-1]->Draw("pE");

    plotGoodRunRange(runIndexStart,runIndexStop,meanZdcX,sigmaZdcX);

    figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_BadRun->Update();
    c_BadRun->Print(figName.c_str());

    // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_ZdcX_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    // c_RunQA_ZdcX->SaveAs(FigName.c_str());

    for(int i_run = 0; i_run < p_mZdcX[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
    {
      if(p_mZdcX[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mZdcX[1][numTriggerBins-1]->GetBinContent(i_run+1),meanZdcX,sigmaZdcX))
      {
	cout << "bad runIndex from ZdcX: " << p_mZdcX[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
	file_badRunIndex << p_mZdcX[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
      }
    }
  }
  //-------------ZdcX----------------

  //-------------Vz----------------
  double meanVz = 0.0;
  double sigmaVz = 0.0;
  findMean(p_mVz[1][numTriggerBins-1],meanVz,sigmaVz);

  // TCanvas *c_RunQA_Vz = new TCanvas("c_RunQA_Vz","c_RunQA_Vz",10,10,800,400);
  // c_RunQA_Vz->cd()->SetLeftMargin(0.1);
  // c_RunQA_Vz->cd()->SetRightMargin(0.1);
  // c_RunQA_Vz->cd()->SetBottomMargin(0.1);
  // c_RunQA_Vz->cd()->SetGrid(0,0);
  // c_RunQA_Vz->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mVz[1][numTriggerBins-1]->SetTitle("Vz vs. runIndex");
  p_mVz[1][numTriggerBins-1]->SetStats(0);
  p_mVz[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mVz[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mVz[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mVz[1][numTriggerBins-1]->SetLineColor(1);
  p_mVz[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mVz[1][numTriggerBins-1]->GetYaxis()->SetTitle("<Vz>");
  p_mVz[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(-3.0,3.0);
  if(beamType == 2) p_mVz[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(200.5,201);
  p_mVz[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanVz,sigmaVz);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_Vz_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_Vz->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mVz[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mVz[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mVz[1][numTriggerBins-1]->GetBinContent(i_run+1),meanVz,sigmaVz))
    {
      cout << "bad runIndex from Vz: " << p_mVz[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mVz[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------Vz----------------

  //-------------Vr----------------
  double meanVr = 0.0;
  double sigmaVr = 0.0;
  findMean(p_mVr[1][numTriggerBins-1],meanVr,sigmaVr);

  // TCanvas *c_RunQA_Vr = new TCanvas("c_RunQA_Vr","c_RunQA_Vr",10,10,800,400);
  // c_RunQA_Vr->cd()->SetLeftMargin(0.1);
  // c_RunQA_Vr->cd()->SetRightMargin(0.1);
  // c_RunQA_Vr->cd()->SetBottomMargin(0.1);
  // c_RunQA_Vr->cd()->SetGrid(0,0);
  // c_RunQA_Vr->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mVr[1][numTriggerBins-1]->SetTitle("Vr vs. runIndex");
  p_mVr[1][numTriggerBins-1]->SetStats(0);
  p_mVr[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mVr[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mVr[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mVr[1][numTriggerBins-1]->SetLineColor(1);
  p_mVr[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mVr[1][numTriggerBins-1]->GetYaxis()->SetTitle("<Vr>");
  p_mVr[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(0.1,0.5);
  p_mVr[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanVr,sigmaVr);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_Vr_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_Vr->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mVr[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mVr[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mVr[1][numTriggerBins-1]->GetBinContent(i_run+1),meanVr,sigmaVr))
    {
      cout << "bad runIndex from Vr: " << p_mVr[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mVr[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------Vr----------------

  //-------------gDca----------------
  double meanGDca = 0.0;
  double sigmaGDca = 0.0;
  findMean(p_mGDca[1][numTriggerBins-1],meanGDca,sigmaGDca);

  // TCanvas *c_RunQA_gDca = new TCanvas("c_RunQA_gDca","c_RunQA_gDca",10,10,800,400);
  // c_RunQA_gDca->cd()->SetLeftMargin(0.1);
  // c_RunQA_gDca->cd()->SetRightMargin(0.1);
  // c_RunQA_gDca->cd()->SetBottomMargin(0.1);
  // c_RunQA_gDca->cd()->SetGrid(0,0);
  // c_RunQA_gDca->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mGDca[1][numTriggerBins-1]->SetTitle("gDca vs. runIndex");
  p_mGDca[1][numTriggerBins-1]->SetStats(0);
  p_mGDca[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mGDca[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mGDca[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mGDca[1][numTriggerBins-1]->SetLineColor(1);
  p_mGDca[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mGDca[1][numTriggerBins-1]->GetYaxis()->SetTitle("<gDca>");
  p_mGDca[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(0.2,1.0);
  p_mGDca[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanGDca,sigmaGDca);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_gDca_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_gDca->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mGDca[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mGDca[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mGDca[1][numTriggerBins-1]->GetBinContent(i_run+1),meanGDca,sigmaGDca))
    {
      cout << "bad runIndex from GDca: " << p_mGDca[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGDca[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------gDca----------------

  //-------------nHitsFit----------------
  double meanNHitsFit = 0.0;
  double sigmaNHitsFit = 0.0;
  findMean(p_mNHitsFit[1][numTriggerBins-1],meanNHitsFit,sigmaNHitsFit);

  // TCanvas *c_RunQA_nHitsFit = new TCanvas("c_RunQA_nHitsFit","c_RunQA_nHitsFit",10,10,800,400);
  // c_RunQA_nHitsFit->cd()->SetLeftMargin(0.1);
  // c_RunQA_nHitsFit->cd()->SetRightMargin(0.1);
  // c_RunQA_nHitsFit->cd()->SetBottomMargin(0.1);
  // c_RunQA_nHitsFit->cd()->SetGrid(0,0);
  // c_RunQA_nHitsFit->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mNHitsFit[1][numTriggerBins-1]->SetTitle("nHitsFit vs. runIndex");
  p_mNHitsFit[1][numTriggerBins-1]->SetStats(0);
  p_mNHitsFit[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mNHitsFit[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mNHitsFit[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mNHitsFit[1][numTriggerBins-1]->SetLineColor(1);
  p_mNHitsFit[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mNHitsFit[1][numTriggerBins-1]->GetYaxis()->SetTitle("<nHitsFit>");
  p_mNHitsFit[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(35,37);
  p_mNHitsFit[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanNHitsFit,sigmaNHitsFit);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_nHitsFit_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_nHitsFit->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mNHitsFit[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mNHitsFit[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mNHitsFit[1][numTriggerBins-1]->GetBinContent(i_run+1),meanNHitsFit,sigmaNHitsFit))
    {
      cout << "bad runIndex from NHitsFit: " << p_mNHitsFit[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mNHitsFit[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------nHitsFit----------------

  //-------------primPt----------------
  double meanPrimPt = 0.0;
  double sigmaPrimPt = 0.0;
  findMean(p_mPrimPt[1][numTriggerBins-1],meanPrimPt,sigmaPrimPt);

  // TCanvas *c_RunQA_primPt = new TCanvas("c_RunQA_primPt","c_RunQA_primPt",10,10,800,400);
  // c_RunQA_primPt->cd()->SetLeftMargin(0.1);
  // c_RunQA_primPt->cd()->SetRightMargin(0.1);
  // c_RunQA_primPt->cd()->SetBottomMargin(0.1);
  // c_RunQA_primPt->cd()->SetGrid(0,0);
  // c_RunQA_primPt->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mPrimPt[1][numTriggerBins-1]->SetTitle("primPt vs. runIndex");
  p_mPrimPt[1][numTriggerBins-1]->SetStats(0);
  p_mPrimPt[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mPrimPt[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mPrimPt[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mPrimPt[1][numTriggerBins-1]->SetLineColor(1);
  p_mPrimPt[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mPrimPt[1][numTriggerBins-1]->GetYaxis()->SetTitle("<p_{T}^{prim}>");
  p_mPrimPt[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(0.63,0.65);
  p_mPrimPt[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanPrimPt,sigmaPrimPt);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_primPt_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_primPt->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mPrimPt[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mPrimPt[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mPrimPt[1][numTriggerBins-1]->GetBinContent(i_run+1),meanPrimPt,sigmaPrimPt))
    {
      cout << "bad runIndex from PrimPt: " << p_mPrimPt[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mPrimPt[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------primPt----------------

  //-------------primEta----------------
  double meanPrimEta = 0.0;
  double sigmaPrimEta = 0.0;
  findMean(p_mPrimEta[1][numTriggerBins-1],meanPrimEta,sigmaPrimEta);

  // TCanvas *c_RunQA_primEta = new TCanvas("c_RunQA_primEta","c_RunQA_primEta",10,10,800,400);
  // c_RunQA_primEta->cd()->SetLeftMargin(0.1);
  // c_RunQA_primEta->cd()->SetRightMargin(0.1);
  // c_RunQA_primEta->cd()->SetBottomMargin(0.1);
  // c_RunQA_primEta->cd()->SetGrid(0,0);
  // c_RunQA_primEta->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mPrimEta[1][numTriggerBins-1]->SetTitle("primEta vs. runIndex");
  p_mPrimEta[1][numTriggerBins-1]->SetStats(0);
  p_mPrimEta[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mPrimEta[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mPrimEta[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mPrimEta[1][numTriggerBins-1]->SetLineColor(1);
  p_mPrimEta[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mPrimEta[1][numTriggerBins-1]->GetYaxis()->SetTitle("<#eta^{prim}>");
  p_mPrimEta[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(-0.05,0.10);
  if(beamType == 2) p_mPrimEta[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(-1.05,-1.0);
  p_mPrimEta[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanPrimEta,sigmaPrimEta);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_primEta_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_primEta->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mPrimEta[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mPrimEta[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mPrimEta[1][numTriggerBins-1]->GetBinContent(i_run+1),meanPrimEta,sigmaPrimEta))
    {
      cout << "bad runIndex from PrimEta: " << p_mPrimEta[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mPrimEta[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------primEta----------------

  //-------------primPhi----------------
  double meanPrimPhi = 0.0;
  double sigmaPrimPhi = 0.0;
  findMean(p_mPrimPhi[1][numTriggerBins-1],meanPrimPhi,sigmaPrimPhi);

  // TCanvas *c_RunQA_primPhi = new TCanvas("c_RunQA_primPhi","c_RunQA_primPhi",10,10,800,400);
  // c_RunQA_primPhi->cd()->SetLeftMargin(0.1);
  // c_RunQA_primPhi->cd()->SetRightMargin(0.1);
  // c_RunQA_primPhi->cd()->SetBottomMargin(0.1);
  // c_RunQA_primPhi->cd()->SetGrid(0,0);
  // c_RunQA_primPhi->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mPrimPhi[1][numTriggerBins-1]->SetTitle("primPhi vs. runIndex");
  p_mPrimPhi[1][numTriggerBins-1]->SetStats(0);
  p_mPrimPhi[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mPrimPhi[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mPrimPhi[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mPrimPhi[1][numTriggerBins-1]->SetLineColor(1);
  p_mPrimPhi[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mPrimPhi[1][numTriggerBins-1]->GetYaxis()->SetTitle("<#phi^{prim}>");
  p_mPrimPhi[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(-0.05,0.05);
  p_mPrimPhi[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanPrimPhi,sigmaPrimPhi);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_primPhi_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_primPhi->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mPrimPhi[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mPrimPhi[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mPrimPhi[1][numTriggerBins-1]->GetBinContent(i_run+1),meanPrimPhi,sigmaPrimPhi))
    {
      cout << "bad runIndex from PrimPhi: " << p_mPrimPhi[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mPrimPhi[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------primPhi----------------

  //-------------globPt----------------
  double meanGlobPt = 0.0;
  double sigmaGlobPt = 0.0;
  findMean(p_mGlobPt[1][numTriggerBins-1],meanGlobPt,sigmaGlobPt);

  // TCanvas *c_RunQA_globPt = new TCanvas("c_RunQA_globPt","c_RunQA_globPt",10,10,800,400);
  // c_RunQA_globPt->cd()->SetLeftMargin(0.1);
  // c_RunQA_globPt->cd()->SetRightMargin(0.1);
  // c_RunQA_globPt->cd()->SetBottomMargin(0.1);
  // c_RunQA_globPt->cd()->SetGrid(0,0);
  // c_RunQA_globPt->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mGlobPt[1][numTriggerBins-1]->SetTitle("globPt vs. runIndex");
  p_mGlobPt[1][numTriggerBins-1]->SetStats(0);
  p_mGlobPt[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mGlobPt[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mGlobPt[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mGlobPt[1][numTriggerBins-1]->SetLineColor(1);
  p_mGlobPt[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mGlobPt[1][numTriggerBins-1]->GetYaxis()->SetTitle("<p_{T}^{glob}>");
  p_mGlobPt[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(0.6,0.7);
  p_mGlobPt[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanGlobPt,sigmaGlobPt);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_globPt_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_globPt->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mGlobPt[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mGlobPt[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mGlobPt[1][numTriggerBins-1]->GetBinContent(i_run+1),meanGlobPt,sigmaGlobPt))
    {
      cout << "bad runIndex from GlobPt: " << p_mGlobPt[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGlobPt[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------globPt----------------

  //-------------globEta----------------
  double meanGlobEta = 0.0;
  double sigmaGlobEta = 0.0;
  findMean(p_mGlobEta[1][numTriggerBins-1],meanGlobEta,sigmaGlobEta);

  // TCanvas *c_RunQA_globEta = new TCanvas("c_RunQA_globEta","c_RunQA_globEta",10,10,800,400);
  // c_RunQA_globEta->cd()->SetLeftMargin(0.1);
  // c_RunQA_globEta->cd()->SetRightMargin(0.1);
  // c_RunQA_globEta->cd()->SetBottomMargin(0.1);
  // c_RunQA_globEta->cd()->SetGrid(0,0);
  // c_RunQA_globEta->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mGlobEta[1][numTriggerBins-1]->SetTitle("globEta vs. runIndex");
  p_mGlobEta[1][numTriggerBins-1]->SetStats(0);
  p_mGlobEta[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mGlobEta[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mGlobEta[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mGlobEta[1][numTriggerBins-1]->SetLineColor(1);
  p_mGlobEta[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mGlobEta[1][numTriggerBins-1]->GetYaxis()->SetTitle("<#eta^{glob}>");
  p_mGlobEta[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(-0.05,0.10);
  if(beamType == 2) p_mGlobEta[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(-1.05,-1.0);
  p_mGlobEta[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanGlobEta,sigmaGlobEta);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_globEta_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_globEta->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mGlobEta[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mGlobEta[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mGlobEta[1][numTriggerBins-1]->GetBinContent(i_run+1),meanGlobEta,sigmaGlobEta))
    {
      cout << "bad runIndex from GlobEta: " << p_mGlobEta[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGlobEta[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------globEta----------------

  //-------------globPhi----------------
  double meanGlobPhi = 0.0;
  double sigmaGlobPhi = 0.0;
  findMean(p_mGlobPhi[1][numTriggerBins-1],meanGlobPhi,sigmaGlobPhi);

  // TCanvas *c_RunQA_globPhi = new TCanvas("c_RunQA_globPhi","c_RunQA_globPhi",10,10,800,400);
  // c_RunQA_globPhi->cd()->SetLeftMargin(0.1);
  // c_RunQA_globPhi->cd()->SetRightMargin(0.1);
  // c_RunQA_globPhi->cd()->SetBottomMargin(0.1);
  // c_RunQA_globPhi->cd()->SetGrid(0,0);
  // c_RunQA_globPhi->cd()->SetTicks(1,1);

  c_BadRun->Clear();
  p_mGlobPhi[1][numTriggerBins-1]->SetTitle("globPhi vs. runIndex");
  p_mGlobPhi[1][numTriggerBins-1]->SetStats(0);
  p_mGlobPhi[1][numTriggerBins-1]->SetMarkerColor(1);
  p_mGlobPhi[1][numTriggerBins-1]->SetMarkerStyle(20);
  p_mGlobPhi[1][numTriggerBins-1]->SetMarkerSize(1.0);
  p_mGlobPhi[1][numTriggerBins-1]->SetLineColor(1);
  p_mGlobPhi[1][numTriggerBins-1]->GetXaxis()->SetTitle("runIndex");
  p_mGlobPhi[1][numTriggerBins-1]->GetYaxis()->SetTitle("<#phi^{glob}>");
  p_mGlobPhi[1][numTriggerBins-1]->GetYaxis()->SetRangeUser(-0.05,0.05);
  p_mGlobPhi[1][numTriggerBins-1]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanGlobPhi,sigmaGlobPhi);

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Update();
  c_BadRun->Print(figName.c_str());

  // FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/RunQA/%s/c_RunQA_globPhi_badRunIndex_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  // c_RunQA_globPhi->SaveAs(FigName.c_str());

  for(int i_run = 0; i_run < p_mGlobPhi[1][numTriggerBins-1]->GetNbinsX(); ++i_run)
  {
    if(p_mGlobPhi[1][numTriggerBins-1]->GetBinError(i_run+1) > 0 && isBadRun(p_mGlobPhi[1][numTriggerBins-1]->GetBinContent(i_run+1),meanGlobPhi,sigmaGlobPhi))
    {
      cout << "bad runIndex from GlobPhi: " << p_mGlobPhi[1][numTriggerBins-1]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGlobPhi[1][numTriggerBins-1]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------globPhi----------------

  figName = Form("../../figures/RunQA/%s/BadRun_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_BadRun->Print(figName.c_str());

  file_badRunIndex.close();
  return 1;
}

void findMean(TProfile *p_runQA, double &mean, double &sigma)
{
  double sum = 0.0;
  double sumSquared = 0.0;
  int counter = 0;

  for(int i_run = 0; i_run < p_runQA->GetNbinsX(); ++i_run)
  {
    if(p_runQA->GetBinError(i_run+1) > 0)
    {
      sum += p_runQA->GetBinContent(i_run+1);
      sumSquared += p_runQA->GetBinContent(i_run+1)*p_runQA->GetBinContent(i_run+1);
      counter++;
    }
  }

  mean = sum/counter;
  sigma = sqrt((sumSquared-(double)counter*mean*mean)/(double)(counter-1));
}

bool isBadRun(double val, double mean, double sigma)
{
  if( val >= mean-3.0*sigma && val <= mean+3.0*sigma) return false;

  return true;
}

void plotGoodRunRange(double runIndexStart, double runIndexStop, double mean, double sigma)
{
  PlotLine(runIndexStart, runIndexStop, mean, mean, 2, 2, 2);
  PlotLine(runIndexStart, runIndexStop, mean+3.0*sigma, mean+3.0*sigma, 4, 2, 2);
  PlotLine(runIndexStart, runIndexStop, mean-3.0*sigma, mean-3.0*sigma, 4, 2, 2);
}
