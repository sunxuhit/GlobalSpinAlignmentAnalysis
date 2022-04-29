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

void rejectPileUpEvent_nTofHit(int energy = 0)
{
  string JobId = "low";
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/test/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/RunQA/merged_file/nTofHitsCut/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
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

  TGraphAsymmErrors *g_UpperBand = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_LowerBand = new TGraphAsymmErrors();
  TH1D *h_projTofHit[100];
  int numTofHit = 0;
  int count = 0;

  string outputname = Form("./figures/projTofHit.pdf");
  TCanvas *c_projTofHit= new TCanvas("c_projTofHit","c_projTofHit",10,10,900,900);
  c_projTofHit->Divide(3,3);
  for(int i_pad = 0; i_pad < 9; ++i_pad)
  {
    c_projTofHit->cd(i_pad+1);
    c_projTofHit->cd(i_pad+1)->SetLeftMargin(0.1);
    c_projTofHit->cd(i_pad+1)->SetRightMargin(0.1);
    c_projTofHit->cd(i_pad+1)->SetBottomMargin(0.1);
    c_projTofHit->cd(i_pad+1)->SetGrid(0,0);
    c_projTofHit->cd(i_pad+1)->SetTicks(1,1);
    c_projTofHit->cd(i_pad+1)->SetLogy(1);
  }

  string output_start = Form("./figures/projTofHit.pdf[");
  c_projTofHit->Print(output_start.c_str());

  while(numTofHit < 3000)
  // while(numTofHit < 100)
  {
    int deltaN = 2;
    if(numTofHit > 100) deltaN = 10;
    if(numTofHit > 300) deltaN = 20;
    if(numTofHit > 500) deltaN = 50;
    if(numTofHit > 1000) deltaN = 100;
    int proj_start = numTofHit;
    int proj_stop  = numTofHit + deltaN;
    int bin_start = h_mTofHitsGRefMult[1][9]->GetXaxis()->FindBin(proj_start);
    int bin_stop  = h_mTofHitsGRefMult[1][9]->GetXaxis()->FindBin(proj_stop);
    // cout << "proj_start = " << proj_start << ", proj_stop = " << proj_stop << endl;
    // cout << "bin_start = " << bin_start << ", bin_stop = " << bin_stop << endl;

    string HistName = Form("h_projTofHit_%d",count);
    h_projTofHit[count] = (TH1D*)h_mTofHitsGRefMult[1][9]->ProjectionY(HistName.c_str(),bin_start,bin_stop);
    string title = Form("nTofHit = [%d,%d]",proj_start,proj_stop);
    h_projTofHit[count]->SetTitle(title.c_str());
    h_projTofHit[count]->SetMarkerStyle(24);
    h_projTofHit[count]->SetMarkerColor(kGray+2);
    h_projTofHit[count]->SetMarkerSize(1.2);
    if(numTofHit < 100) h_projTofHit[count]->GetXaxis()->SetRangeUser(0,100.0);
    else if(numTofHit < 200) h_projTofHit[count]->GetXaxis()->SetRangeUser(0,200.0);
    else if(numTofHit < 500) h_projTofHit[count]->GetXaxis()->SetRangeUser(0,400.0);
    else h_projTofHit[count]->GetXaxis()->SetRangeUser(0,800.0);

    int nPad = count%9;
    c_projTofHit->cd(nPad+1);
    h_projTofHit[count]->Draw("pE");
    h_projTofHit[count]->Fit("gaus");
    TF1 *f_gaus = h_projTofHit[count]->GetFunction("gaus");
    double mean = f_gaus->GetParameter(1);
    double sigma = f_gaus->GetParameter(2);
    // cout << "Norm = " << f_gaus->GetParameter(0) << ", Mean = " << f_gaus->GetParameter(1) << ", Sigma = " << f_gaus->GetParameter(2) << endl;
    g_UpperBand->SetPoint(count,0.5*(proj_start+proj_stop),mean+4.0*sigma);
    g_LowerBand->SetPoint(count,0.5*(proj_start+proj_stop),mean-4.0*sigma);

    if(nPad == 8)
    {
      c_projTofHit->Update();
      c_projTofHit->Print(outputname.c_str());
      for(int i_pad = 0; i_pad < 9; ++i_pad)
      {
	c_projTofHit->cd(i_pad+1)->Clear();
      }
    }

    numTofHit = proj_stop + 1;
    count++;
    // cout << "proj_start = " << proj_start << ", proj_stop = " << proj_stop << ", numTofHit = " << numTofHit << endl;
  }

  c_projTofHit->Update();
  c_projTofHit->Print(outputname.c_str());

  string output_stop = Form("./figures/projTofHit.pdf]");
  c_projTofHit->Print(output_stop.c_str());

  cout << "count = " << count << endl;
  g_UpperBand->RemovePoint(g_UpperBand->GetN()-1);
  g_UpperBand->RemovePoint(g_UpperBand->GetN()-1);
  g_UpperBand->RemovePoint(0);
  g_LowerBand->RemovePoint(g_LowerBand->GetN()-1);
  g_LowerBand->RemovePoint(g_LowerBand->GetN()-1);
  g_LowerBand->RemovePoint(0);

  TCanvas *c_PileUp = new TCanvas("c_PileUp","c_PileUp",10,10,800,800);
  c_PileUp->cd();
  c_PileUp->cd()->SetLeftMargin(0.15);
  c_PileUp->cd()->SetBottomMargin(0.15);
  c_PileUp->cd()->SetGrid(0,0);
  c_PileUp->cd()->SetTicks(1,1);
  c_PileUp->cd()->SetLogz(1);
  TH1D *h_play = new TH1D("h_play","h_play",4000,-0.5,3999.5);
  for(int i_bin = 0; i_bin < 4000; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-100.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("gRefMult vs. nTofHits");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("nTofHits");
  h_play->GetYaxis()->SetTitle("gRefMult");
  h_play->GetYaxis()->SetRangeUser(-30.0,1000.0);
  h_play->Draw("hE");

  h_mTofHitsGRefMult[1][9]->Draw("colz same");

  g_UpperBand->SetMarkerStyle(24);
  g_UpperBand->SetMarkerColor(2);
  g_UpperBand->SetMarkerSize(0.8);
  g_UpperBand->Draw("pE Same");
  TF1 *f_UpperBand = new TF1("f_UpperBand","pol5",0,4000);
  f_UpperBand->SetRange(0.0,3000.0);
  g_UpperBand->Fit(f_UpperBand,"RN");
  TF1 *f_UpperExt = new TF1("f_UpperExt","pol1",2500,4000);
  f_UpperExt->SetRange(2500.0,3000.0);
  g_UpperBand->Fit(f_UpperExt,"RN");

  f_UpperBand->SetLineColor(kRed+2);
  f_UpperBand->SetLineStyle(1);
  f_UpperBand->SetRange(0.0,2800.0);
  f_UpperBand->Draw("l Same");

  f_UpperExt->SetLineColor(kRed+2);
  f_UpperExt->SetLineStyle(2);
  f_UpperExt->SetRange(2500.0,4000.0);
  f_UpperExt->Draw("l same");

  cout << "UpperBand: " << endl;
  cout << "Pol5 fit" << endl;
  cout << "p0Upper = " << f_UpperBand->GetParameter(0) << endl;
  cout << "p1Upper = " << f_UpperBand->GetParameter(1) << endl;
  cout << "p2Upper = " << f_UpperBand->GetParameter(2) << endl;
  cout << "p3Upper = " << f_UpperBand->GetParameter(3) << endl;
  cout << "p4Upper = " << f_UpperBand->GetParameter(4) << endl;
  cout << "p5Upper = " << f_UpperBand->GetParameter(5) << endl;
  cout << "Pol1 fit" << endl;
  cout << "p0UpperExt = " << f_UpperExt->GetParameter(0) << endl;
  cout << "p1UpperExt = " << f_UpperExt->GetParameter(1) << endl;

  g_LowerBand->SetMarkerStyle(24);
  g_LowerBand->SetMarkerColor(kGray+2);
  g_LowerBand->SetMarkerSize(0.8);
  g_LowerBand->Draw("PE Same");
  TF1 *f_LowerBand = new TF1("f_LowerBand","pol5",0,4000);
  f_LowerBand->SetRange(50.0,3000.0);
  g_LowerBand->Fit(f_LowerBand,"RN");
  TF1 *f_LowerExt = new TF1("f_LowerExt","pol1",2500,4000);
  f_LowerExt->SetRange(2500.0,3000.0);
  g_LowerBand->Fit(f_LowerExt,"RN");

  f_LowerBand->SetLineColor(kGray+2);
  f_LowerBand->SetLineStyle(1);
  f_LowerBand->SetRange(0.0,2800.0);
  f_LowerBand->Draw("l Same");

  f_LowerExt->SetLineColor(kGray+2);
  f_LowerExt->SetLineStyle(2);
  f_LowerExt->SetRange(2500.0,4000.0);
  f_LowerExt->Draw("l same");

  cout << "LowerBand: " << endl;
  cout << "Pol5 fit" << endl;
  cout << "p0Lower = " << f_LowerBand->GetParameter(0) << endl;
  cout << "p1Lower = " << f_LowerBand->GetParameter(1) << endl;
  cout << "p2Lower = " << f_LowerBand->GetParameter(2) << endl;
  cout << "p3Lower = " << f_LowerBand->GetParameter(3) << endl;
  cout << "p4Lower = " << f_LowerBand->GetParameter(4) << endl;
  cout << "p5Lower = " << f_LowerBand->GetParameter(5) << endl;
  cout << "Pol1 fit" << endl;
  cout << "p0LowerExt = " << f_LowerExt->GetParameter(0) << endl;
  cout << "p1LowerExt = " << f_LowerExt->GetParameter(1) << endl;

  // c_PileUp->SaveAs("c_PileUp_nTofHit.eps");
  c_PileUp->SaveAs("./figures/c_PileUp_nTofHit.png");
}
