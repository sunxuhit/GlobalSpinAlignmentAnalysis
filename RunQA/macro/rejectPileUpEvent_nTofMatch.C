#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "../StRoot/StRunQAMaker/StRunQACons.h"

using namespace std;

static const string mCutsQA[2] = {"Before","After"};

double doubleGaus(double *var, double *par)
{
  double x = var[0];

  double norm1  = par[0];
  double mean1  = par[1];
  double sigma1 = par[2];
  double power1 = -0.5*(x-mean1)*(x-mean1)/(sigma1*sigma1);
  double y1 = norm1*exp(power1)/(sigma1*sqrt(TMath::TwoPi()));
  // double y1 = norm1*exp(power1);

  double norm2  = par[3];
  double mean2  = par[4];
  double sigma2 = par[5];
  double power2 = -0.5*(x-mean2)*(x-mean2)/(sigma2*sigma2);
  double y2 = norm2*exp(power2)/(sigma2*sqrt(TMath::TwoPi()));
  // double y2 = norm2*exp(power2);

  double y = y1 + y2;

  return y;
}

void rejectPileUpEvent_nTofMatch(int energy = 0)
{
  string JobId = "low";
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/test/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/RunQA/merged_file/nTofMatchCut/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
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
  TH1D *h_projTofMatch[100];
  int numTofMatch = 0;
  int count = 0;

  string outputname = Form("./figures/projTofMatch.pdf");
  TCanvas *c_projTofMatch= new TCanvas("c_projTofMatch","c_projTofMatch",10,10,900,900);
  c_projTofMatch->Divide(3,3);
  for(int i_pad = 0; i_pad < 9; ++i_pad)
  {
    c_projTofMatch->cd(i_pad+1);
    c_projTofMatch->cd(i_pad+1)->SetLeftMargin(0.1);
    c_projTofMatch->cd(i_pad+1)->SetRightMargin(0.1);
    c_projTofMatch->cd(i_pad+1)->SetBottomMargin(0.1);
    c_projTofMatch->cd(i_pad+1)->SetGrid(0,0);
    c_projTofMatch->cd(i_pad+1)->SetTicks(1,1);
    c_projTofMatch->cd(i_pad+1)->SetLogy(1);
  }

  string output_start = Form("./figures/projTofMatch.pdf[");
  c_projTofMatch->Print(output_start.c_str());

  while(numTofMatch < 800)
  {
    int deltaN = 2;
    if(numTofMatch > 100) deltaN = 10;
    if(numTofMatch > 400) deltaN = 20;
    if(numTofMatch > 600) deltaN = 50;
    int proj_start = numTofMatch;
    int proj_stop  = numTofMatch + deltaN;
    int bin_start = h_mTofMatchGRefMult[1][9]->GetXaxis()->FindBin(proj_start);
    int bin_stop  = h_mTofMatchGRefMult[1][9]->GetXaxis()->FindBin(proj_stop);
    // cout << "proj_start = " << proj_start << ", proj_stop = " << proj_stop << endl;
    // cout << "bin_start = " << bin_start << ", bin_stop = " << bin_stop << endl;

    string HistName = Form("h_projTofMatch_%d",count);
    h_projTofMatch[count] = (TH1D*)h_mTofMatchGRefMult[1][9]->ProjectionY(HistName.c_str(),bin_start,bin_stop);
    string title = Form("nTofMatch = [%d,%d]",proj_start,proj_stop);
    h_projTofMatch[count]->SetTitle(title.c_str());
    h_projTofMatch[count]->SetMarkerStyle(24);
    h_projTofMatch[count]->SetMarkerColor(kGray+2);
    h_projTofMatch[count]->SetMarkerSize(1.2);
    if(numTofMatch < 10) h_projTofMatch[count]->GetXaxis()->SetRangeUser(0,200.0);
    else if(numTofMatch < 100) h_projTofMatch[count]->GetXaxis()->SetRangeUser(0,600.0);
    else h_projTofMatch[count]->GetXaxis()->SetRangeUser(0,800.0);

    int nPad = count%9;
    c_projTofMatch->cd(nPad+1);
    h_projTofMatch[count]->Draw("pE");
    if(h_projTofMatch[count]->GetEntries() > 0)
    {
      if(numTofMatch < 400) 
      {
	h_projTofMatch[count]->Fit("gaus");
	TF1 *f_gaus = h_projTofMatch[count]->GetFunction("gaus");
	double mean = f_gaus->GetParameter(1);
	double sigma = f_gaus->GetParameter(2);
	// cout << "Norm = " << f_gaus->GetParameter(0) << ", Mean = " << f_gaus->GetParameter(1) << ", Sigma = " << f_gaus->GetParameter(2) << endl;
	g_UpperBand->SetPoint(count,0.5*(proj_start+proj_stop),mean+4.0*sigma);
	g_LowerBand->SetPoint(count,0.5*(proj_start+proj_stop),mean-4.0*sigma);
      }
      /*
      else if(numTofMatch < 450)
      {
	TF1 *f_gaus = new TF1("f_gaus","gaus",0,1000.0);
	f_gaus->SetRange(numTofMatch-50,800.0);
	h_projTofMatch[count]->Fit(f_gaus,"R");
	double mean = f_gaus->GetParameter(1);
	double sigma = f_gaus->GetParameter(2);
	// cout << "Norm = " << f_gaus->GetParameter(0) << ", Mean = " << f_gaus->GetParameter(1) << ", Sigma = " << f_gaus->GetParameter(2) << endl;
	g_UpperBand->SetPoint(count,0.5*(proj_start+proj_stop),mean+4.0*sigma);
	g_LowerBand->SetPoint(count,0.5*(proj_start+proj_stop),mean-3.5*sigma);
      }
      */
      else if(numTofMatch < 500)
      {
	TF1 *f_gaus = new TF1("f_gaus",doubleGaus,0,1000.0,6);
	for(int i_par = 0; i_par < 6; ++i_par)
	{
	  f_gaus->ReleaseParameter(i_par);
	}
	f_gaus->SetParameter(0,100.0);
	f_gaus->SetParameter(1,numTofMatch-50);
	f_gaus->SetParLimits(1,300.0,numTofMatch);
	f_gaus->SetParameter(2,10.0);
	f_gaus->SetParameter(3,100.0);
	f_gaus->SetParameter(4,numTofMatch+50);
	f_gaus->SetParLimits(4,numTofMatch,700.0);
	f_gaus->SetParameter(5,20.0);
	h_projTofMatch[count]->Fit(f_gaus);
	double mean = f_gaus->GetParameter(4);
	double sigma = f_gaus->GetParameter(5);
	// cout << "Norm = " << f_gaus->GetParameter(0) << ", Mean = " << f_gaus->GetParameter(1) << ", Sigma = " << f_gaus->GetParameter(2) << endl;
	g_UpperBand->SetPoint(count,0.5*(proj_start+proj_stop),mean+4.0*sigma);
	g_LowerBand->SetPoint(count,0.5*(proj_start+proj_stop),mean-4.0*sigma);
      }
    }

    if(nPad == 8)
    {
      c_projTofMatch->Update();
      c_projTofMatch->Print(outputname.c_str());
      for(int i_pad = 0; i_pad < 9; ++i_pad)
      {
	c_projTofMatch->cd(i_pad+1)->Clear();
      }
    }

    numTofMatch = proj_stop + 1;
    count++;
    // cout << "proj_start = " << proj_start << ", proj_stop = " << proj_stop << ", numTofMatch = " << numTofMatch << endl;
  }

  c_projTofMatch->Update();
  c_projTofMatch->Print(outputname.c_str());

  string output_stop = Form("./figures/projTofMatch.pdf]");
  c_projTofMatch->Print(output_stop.c_str());

  cout << "count = " << count << endl;

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
  h_play->SetTitle("gRefMult vs. nTofMatch");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("nTofMatch");
  h_play->GetXaxis()->SetRangeUser(-30.0,800.0);
  h_play->GetYaxis()->SetTitle("gRefMult");
  h_play->GetYaxis()->SetRangeUser(-30.0,800.0);
  h_play->Draw("hE");

  h_mTofMatchGRefMult[1][9]->Draw("colz same");

  g_UpperBand->SetMarkerStyle(24);
  g_UpperBand->SetMarkerColor(2);
  g_UpperBand->SetMarkerSize(0.8);
  g_UpperBand->Draw("pE Same");

  TF1 *f_UpperBand = new TF1("f_UpperBand","pol5",0,4000);
  f_UpperBand->SetRange(0.0,500.0);
  g_UpperBand->Fit(f_UpperBand,"RN");
  f_UpperBand->SetLineColor(kRed+2);
  f_UpperBand->SetLineStyle(1);
  f_UpperBand->SetRange(0.0,500.0);
  f_UpperBand->Draw("l Same");

  TF1 *f_UpperExt = new TF1("f_UpperExt","pol1",0,800);
  f_UpperExt->SetRange(400.0,550.0);
  g_UpperBand->Fit(f_UpperExt,"RN");
  f_UpperExt->SetLineColor(kRed+2);
  f_UpperExt->SetLineStyle(2);
  f_UpperExt->SetRange(400.0,800.0);
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
  f_LowerBand->SetRange(0.0,500.0);
  g_LowerBand->Fit(f_LowerBand,"RN");
  f_LowerBand->SetLineColor(kGray+2);
  f_LowerBand->SetLineStyle(1);
  f_LowerBand->SetRange(0.0,500.0);
  f_LowerBand->Draw("l Same");

  TF1 *f_LowerExt = new TF1("f_LowerExt","pol1",0,800);
  f_LowerExt->SetRange(400.0,550.0);
  g_LowerBand->Fit(f_LowerExt,"RN");
  f_LowerExt->SetLineColor(kGray+2);
  f_LowerExt->SetLineStyle(2);
  f_LowerExt->SetRange(400.0,800.0);
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

  // c_PileUp->SaveAs("./figures/c_PileUp_nTofMatch.eps");
  c_PileUp->SaveAs("./figures/c_PileUp_nTofMatch.png");
}
