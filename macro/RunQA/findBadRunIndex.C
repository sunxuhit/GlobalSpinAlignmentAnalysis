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

#include "../StRoot/StRunQAMaker/StRunQACons.h"
#include "./draw.h"

using namespace std;

static const string mCutsQA[2] = {"Before","After"};

void findMean(TProfile *p_runQA, double &mean, double &sigma);
bool isBadRun(double val, double mean, double sigma);
void plotGoodRunRange(double runIndexStart, double runIndexStop, double mean, double sigma);

int findBadRunIndex(int energy = 0)
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

  string outputfile = Form("../StRoot/StRunQAUtility/RunIndex/badRunIndexUnSorted_%s.txt",runQA::mBeamEnergy[energy].c_str());
  ofstream file_badRunIndex;
  file_badRunIndex.open(outputfile.c_str());
  if (!file_badRunIndex.is_open()) 
  {
    std::cout << "failed to open " << outputfile << '\n';
    return -1;
  } 

  const double runIndexStart = -0.5;
  const double runIndexStop  = 3999.5;

  const int numOfTriggers = 5;
  const int triggerID[numOfTriggers] = {450005,450015,450025,450050,450060};
  const int MarkerColor[numOfTriggers] = {7,6,4,2,16};

  //-------------refMult----------------
  double meanRefMult = 0.0;
  double sigmaRefMult = 0.0;
  findMean(p_mRefMult[1][9],meanRefMult,sigmaRefMult);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanRefMult,sigmaRefMult);

  c_RunQA_RefMult->SaveAs("./figures/c_RunQA_RefMult_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mRefMult[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mRefMult[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mRefMult[1][9]->GetBinContent(i_run+1),meanRefMult,sigmaRefMult))
    {
      cout << "bad runIndex from refMult: " << p_mRefMult[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mRefMult[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------refMult----------------

  //-------------grefMult----------------
  double meanGRefMult = 0.0;
  double sigmaGRefMult = 0.0;
  findMean(p_mGRefMult[1][9],meanGRefMult,sigmaGRefMult);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanGRefMult,sigmaGRefMult);

  c_RunQA_gRefMult->SaveAs("./figures/c_RunQA_gRefMult_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mGRefMult[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mGRefMult[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mGRefMult[1][9]->GetBinContent(i_run+1),meanGRefMult,sigmaGRefMult))
    {
      cout << "bad runIndex from grefMult: " << p_mGRefMult[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGRefMult[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------grefMult----------------

  //-------------ZdcX----------------
  double meanZdcX = 0.0;
  double sigmaZdcX = 0.0;
  findMean(p_mZdcX[1][9],meanZdcX,sigmaZdcX);

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
  p_mZdcX[1][9]->GetYaxis()->SetRangeUser(0,60000);
  p_mZdcX[1][9]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanZdcX,sigmaZdcX);

  c_RunQA_ZdcX->SaveAs("./figures/c_RunQA_ZdcX_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mZdcX[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mZdcX[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mZdcX[1][9]->GetBinContent(i_run+1),meanZdcX,sigmaZdcX))
    {
      cout << "bad runIndex from ZdcX: " << p_mZdcX[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mZdcX[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------ZdcX----------------

  //-------------Vz----------------
  double meanVz = 0.0;
  double sigmaVz = 0.0;
  findMean(p_mVz[1][9],meanVz,sigmaVz);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanVz,sigmaVz);

  c_RunQA_Vz->SaveAs("./figures/c_RunQA_Vz_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mVz[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mVz[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mVz[1][9]->GetBinContent(i_run+1),meanVz,sigmaVz))
    {
      cout << "bad runIndex from Vz: " << p_mVz[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mVz[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------Vz----------------

  //-------------Vr----------------
  double meanVr = 0.0;
  double sigmaVr = 0.0;
  findMean(p_mVr[1][9],meanVr,sigmaVr);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanVr,sigmaVr);

  c_RunQA_Vr->SaveAs("./figures/c_RunQA_Vr_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mVr[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mVr[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mVr[1][9]->GetBinContent(i_run+1),meanVr,sigmaVr))
    {
      cout << "bad runIndex from Vr: " << p_mVr[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mVr[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------Vr----------------

  //-------------gDca----------------
  double meanGDca = 0.0;
  double sigmaGDca = 0.0;
  findMean(p_mGDca[1][9],meanGDca,sigmaGDca);

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
  p_mGDca[1][9]->GetYaxis()->SetRangeUser(0.1,0.7);
  p_mGDca[1][9]->Draw("pE");

  plotGoodRunRange(runIndexStart,runIndexStop,meanGDca,sigmaGDca);

  c_RunQA_gDca->SaveAs("./figures/c_RunQA_gDca_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mGDca[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mGDca[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mGDca[1][9]->GetBinContent(i_run+1),meanGDca,sigmaGDca))
    {
      cout << "bad runIndex from GDca: " << p_mGDca[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGDca[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------gDca----------------

  //-------------nHitsFit----------------
  double meanNHitsFit = 0.0;
  double sigmaNHitsFit = 0.0;
  findMean(p_mNHitsFit[1][9],meanNHitsFit,sigmaNHitsFit);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanNHitsFit,sigmaNHitsFit);

  c_RunQA_nHitsFit->SaveAs("./figures/c_RunQA_nHitsFit_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mNHitsFit[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mNHitsFit[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mNHitsFit[1][9]->GetBinContent(i_run+1),meanNHitsFit,sigmaNHitsFit))
    {
      cout << "bad runIndex from NHitsFit: " << p_mNHitsFit[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mNHitsFit[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------nHitsFit----------------

  //-------------primPt----------------
  double meanPrimPt = 0.0;
  double sigmaPrimPt = 0.0;
  findMean(p_mPrimPt[1][9],meanPrimPt,sigmaPrimPt);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanPrimPt,sigmaPrimPt);

  c_RunQA_primPt->SaveAs("./figures/c_RunQA_primPt_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mPrimPt[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mPrimPt[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mPrimPt[1][9]->GetBinContent(i_run+1),meanPrimPt,sigmaPrimPt))
    {
      cout << "bad runIndex from PrimPt: " << p_mPrimPt[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mPrimPt[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------primPt----------------

  //-------------primEta----------------
  double meanPrimEta = 0.0;
  double sigmaPrimEta = 0.0;
  findMean(p_mPrimEta[1][9],meanPrimEta,sigmaPrimEta);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanPrimEta,sigmaPrimEta);

  c_RunQA_primEta->SaveAs("./figures/c_RunQA_primEta_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mPrimEta[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mPrimEta[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mPrimEta[1][9]->GetBinContent(i_run+1),meanPrimEta,sigmaPrimEta))
    {
      cout << "bad runIndex from PrimEta: " << p_mPrimEta[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mPrimEta[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------primEta----------------

  //-------------primPhi----------------
  double meanPrimPhi = 0.0;
  double sigmaPrimPhi = 0.0;
  findMean(p_mPrimPhi[1][9],meanPrimPhi,sigmaPrimPhi);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanPrimPhi,sigmaPrimPhi);

  c_RunQA_primPhi->SaveAs("./figures/c_RunQA_primPhi_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mPrimPhi[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mPrimPhi[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mPrimPhi[1][9]->GetBinContent(i_run+1),meanPrimPhi,sigmaPrimPhi))
    {
      cout << "bad runIndex from PrimPhi: " << p_mPrimPhi[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mPrimPhi[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------primPhi----------------

  //-------------globPt----------------
  double meanGlobPt = 0.0;
  double sigmaGlobPt = 0.0;
  findMean(p_mGlobPt[1][9],meanGlobPt,sigmaGlobPt);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanGlobPt,sigmaGlobPt);

  c_RunQA_globPt->SaveAs("./figures/c_RunQA_globPt_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mGlobPt[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mGlobPt[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mGlobPt[1][9]->GetBinContent(i_run+1),meanGlobPt,sigmaGlobPt))
    {
      cout << "bad runIndex from GlobPt: " << p_mGlobPt[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGlobPt[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------globPt----------------

  //-------------globEta----------------
  double meanGlobEta = 0.0;
  double sigmaGlobEta = 0.0;
  findMean(p_mGlobEta[1][9],meanGlobEta,sigmaGlobEta);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanGlobEta,sigmaGlobEta);

  c_RunQA_globEta->SaveAs("./figures/c_RunQA_globEta_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mGlobEta[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mGlobEta[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mGlobEta[1][9]->GetBinContent(i_run+1),meanGlobEta,sigmaGlobEta))
    {
      cout << "bad runIndex from GlobEta: " << p_mGlobEta[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGlobEta[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------globEta----------------

  //-------------globPhi----------------
  double meanGlobPhi = 0.0;
  double sigmaGlobPhi = 0.0;
  findMean(p_mGlobPhi[1][9],meanGlobPhi,sigmaGlobPhi);

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

  plotGoodRunRange(runIndexStart,runIndexStop,meanGlobPhi,sigmaGlobPhi);

  c_RunQA_globPhi->SaveAs("./figures/c_RunQA_globPhi_badRunIndex.eps");

  for(int i_run = 0; i_run < p_mGlobPhi[1][9]->GetNbinsX(); ++i_run)
  {
    if(p_mGlobPhi[1][9]->GetBinError(i_run+1) > 0 && isBadRun(p_mGlobPhi[1][9]->GetBinContent(i_run+1),meanGlobPhi,sigmaGlobPhi))
    {
      cout << "bad runIndex from GlobPhi: " << p_mGlobPhi[1][9]->GetBinCenter(i_run+1) << endl;
      file_badRunIndex << p_mGlobPhi[1][9]->GetBinCenter(i_run+1) << std::endl;
    }
  }
  //-------------globPhi----------------

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
  PlotLine(runIndexStart, runIndexStop, mean+3.0*sigma, mean+3.0*sigma, 4, 2, 1);
  PlotLine(runIndexStart, runIndexStop, mean-3.0*sigma, mean-3.0*sigma, 4, 2, 1);
}
