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

static const string CutsQA[2] = {"Before","After"};

void rejectPileUpEvent_nTofHit(int energy = 0)
{
  string JobId = "low";
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/Data/SpinAlignment/AuAu%s/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  // string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/test/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH2F *h_mRefMultGRefMult[2];
  TH2F *h_mRefMultTofMatch[2];
  TH2F *h_mRefMultTofHits[2];
  TH2F *h_mGRefMultTofMatch[2];
  TH2F *h_mGRefMultTofHits[2];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    // string HistName = Form("h_mRefMultGRefMult%s_trigger9",mCutsQA[i_cut].c_str());
    string HistName = Form("h_mRefMultGRefMult_%s",CutsQA[i_cut].c_str());
    h_mRefMultGRefMult[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mRefMultGRefMult[i_cut]->GetXaxis()->SetTitle("refMult");
    h_mRefMultGRefMult[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mRefMultGRefMult[i_cut]->GetYaxis()->SetTitle("grefMult");
    h_mRefMultGRefMult[i_cut]->GetYaxis()->SetRangeUser(0.0,800.0);

    HistName = Form("h_mRefMultTofMatch_%s",CutsQA[i_cut].c_str());
    h_mRefMultTofMatch[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mRefMultTofMatch[i_cut]->GetXaxis()->SetTitle("refMult");
    h_mRefMultTofMatch[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mRefMultTofMatch[i_cut]->GetYaxis()->SetTitle("ToFMatch");
    h_mRefMultTofMatch[i_cut]->GetYaxis()->SetRangeUser(0.0,1000.0);

    HistName = Form("h_mRefMultTofHits_%s",CutsQA[i_cut].c_str());
    h_mRefMultTofHits[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mRefMultTofHits[i_cut]->GetXaxis()->SetTitle("refMult");
    h_mRefMultTofHits[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mRefMultTofHits[i_cut]->GetYaxis()->SetTitle("ToFHits");
    h_mRefMultTofHits[i_cut]->GetYaxis()->SetRangeUser(0.0,3500.0);

    HistName = Form("h_mGRefMultTofMatch_%s",CutsQA[i_cut].c_str());
    h_mGRefMultTofMatch[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mGRefMultTofMatch[i_cut]->GetXaxis()->SetTitle("gRefMult");
    h_mGRefMultTofMatch[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mGRefMultTofMatch[i_cut]->GetYaxis()->SetTitle("ToFMatch");
    h_mGRefMultTofMatch[i_cut]->GetYaxis()->SetRangeUser(0.0,1000.0);

    HistName = Form("h_mGRefMultTofHits_%s",CutsQA[i_cut].c_str());
    h_mGRefMultTofHits[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    h_mGRefMultTofHits[i_cut]->GetXaxis()->SetTitle("gRefMult");
    h_mGRefMultTofHits[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mGRefMultTofHits[i_cut]->GetYaxis()->SetTitle("ToFHits");
    h_mGRefMultTofHits[i_cut]->GetYaxis()->SetRangeUser(0.0,3500.0);
  }

  /*
  TCanvas *c_EventQA = new TCanvas("c_EventQA","c_EventQA",10,10,1500,600);
  c_EventQA->Divide(5,2);
  for(int i_pad = 0; i_pad < 10; ++i_pad)
  {
    c_EventQA->cd(i_pad+1);
    c_EventQA->cd(i_pad+1)->SetLeftMargin(0.1);
    c_EventQA->cd(i_pad+1)->SetRightMargin(0.1);
    c_EventQA->cd(i_pad+1)->SetBottomMargin(0.1);
    c_EventQA->cd(i_pad+1)->SetGrid(0,0);
    c_EventQA->cd(i_pad+1)->SetTicks(1,1);
    c_EventQA->cd(i_pad+1)->SetLogz(1);
  }

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    c_EventQA->cd(i_cut*5+1);
    h_mRefMultGRefMult[i_cut]->Draw("colz");

    c_EventQA->cd(i_cut*5+2);
    h_mRefMultTofMatch[i_cut]->Draw("colz");

    c_EventQA->cd(i_cut*5+3);
    h_mRefMultTofHits[i_cut]->Draw("colz");

    c_EventQA->cd(i_cut*5+4);
    h_mGRefMultTofMatch[i_cut]->Draw("colz");

    c_EventQA->cd(i_cut*5+5);
    h_mGRefMultTofHits[i_cut]->Draw("colz");
  }
  */

  TGraphAsymmErrors *g_UpperBand = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_LowerBand = new TGraphAsymmErrors();
  TH1D *h_projTofHit[100];
  int numTofHit = 0;
  int count = 0;

  string outputname = Form("./projTofHit.pdf");
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

  string output_start = Form("./projTofHit.pdf[");
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
    int bin_start = h_mGRefMultTofHits[1]->GetYaxis()->FindBin(proj_start);
    int bin_stop  = h_mGRefMultTofHits[1]->GetYaxis()->FindBin(proj_stop);
    // cout << "proj_start = " << proj_start << ", proj_stop = " << proj_stop << endl;
    // cout << "bin_start = " << bin_start << ", bin_stop = " << bin_stop << endl;

    string HistName = Form("h_projTofHit_%d",count);
    h_projTofHit[count] = (TH1D*)h_mGRefMultTofHits[1]->ProjectionX(HistName.c_str(),bin_start,bin_stop);
    string title = Form("nTofHit = [%d,%d]",proj_start,proj_stop);
    h_projTofHit[count]->SetTitle(title.c_str());
    h_projTofHit[count]->SetMarkerStyle(24);
    h_projTofHit[count]->SetMarkerColor(kGray+2);
    h_projTofHit[count]->SetMarkerSize(1.2);
    if(numTofHit < 100) h_projTofHit[count]->GetXaxis()->SetRangeUser(0,100.0);
    else if(numTofHit < 200) h_projTofHit[count]->GetXaxis()->SetRangeUser(0,200.0);
    else if(numTofHit < 500) h_projTofHit[count]->GetXaxis()->SetRangeUser(0,400.0);

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

  string output_stop = Form("./projTofHit.pdf]");
  c_projTofHit->Print(output_stop.c_str());

  cout << "count = " << count << endl;

  TCanvas *c_PileUp = new TCanvas("c_PileUp","c_PileUp",10,10,800,800);
  c_PileUp->cd();
  c_PileUp->cd()->SetLeftMargin(0.1);
  c_PileUp->cd()->SetRightMargin(0.1);
  c_PileUp->cd()->SetBottomMargin(0.1);
  c_PileUp->cd()->SetGrid(0,0);
  c_PileUp->cd()->SetTicks(1,1);
  TH1D *h_play = new TH1D("h_play","h_play",4000,-0.5,3999.5);
  for(int i_bin = 0; i_bin < 4000; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-100.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->GetXaxis()->SetTitle("nTofHits");
  h_play->GetYaxis()->SetTitle("gRefMult");
  h_play->GetYaxis()->SetRangeUser(-30.0,1000.0);
  h_play->Draw("hE");

  g_UpperBand->SetMarkerStyle(20);
  g_UpperBand->SetMarkerColor(kRed+2);
  g_UpperBand->SetMarkerSize(0.8);
  g_UpperBand->Draw("pE Same");
  g_LowerBand->SetMarkerStyle(20);
  g_LowerBand->SetMarkerColor(kGray+2);
  g_LowerBand->SetMarkerSize(0.8);
  g_LowerBand->Draw("PE Same");
}
