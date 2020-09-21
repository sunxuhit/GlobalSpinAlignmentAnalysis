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

void rejectPileUpEvent_nTofMatch(int energy = 0)
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

  TH1D *h_projTofMatch[300];
  int numTofMatch = 0;
  int count = 0;

  string outputname = Form("./projTofMatch.pdf");
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

  string output_start = Form("./projTofMatch.pdf[");
  c_projTofMatch->Print(output_start.c_str());

  while(numTofMatch < 1000)
  {
    int deltaN = 2;
    if(numTofMatch > 100) deltaN = 10;
    if(numTofMatch > 300) deltaN = 20;
    if(numTofMatch > 500) deltaN = 50;
    int proj_start = numTofMatch;
    int proj_stop  = numTofMatch + deltaN;
    int bin_start = h_mGRefMultTofMatch[1]->GetYaxis()->FindBin(proj_start);
    int bin_stop  = h_mGRefMultTofMatch[1]->GetYaxis()->FindBin(proj_stop);
    // cout << "proj_start = " << proj_start << ", proj_stop = " << proj_stop << endl;
    // cout << "bin_start = " << bin_start << ", bin_stop = " << bin_stop << endl;

    string HistName = Form("h_projTofMatch_%d",count);
    h_projTofMatch[count] = (TH1D*)h_mGRefMultTofMatch[1]->ProjectionX(HistName.c_str(),bin_start,bin_stop);
    string title = Form("nTofMatch = [%d,%d]",proj_start,proj_stop);
    h_projTofMatch[count]->SetTitle(title.c_str());
    h_projTofMatch[count]->SetMarkerStyle(24);
    h_projTofMatch[count]->SetMarkerColor(kGray+2);
    h_projTofMatch[count]->SetMarkerSize(1.2);

    int nPad = count%9;
    c_projTofMatch->cd(nPad+1);
    h_projTofMatch[count]->Draw("pE");

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

  string output_stop = Form("./projTofMatch.pdf]");
  c_projTofMatch->Print(output_stop.c_str());

  cout << "count = " << count << endl;
}
