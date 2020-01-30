#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "../StRoot/StVecMesonMaker/StVecMesonCons.h"

using namespace std;

static const string CutsQA[2] = {"Before","After"};

void plotQA_Event(int energy = 2)
{
  string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/QA/file_%s_QA.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mRefMult[2]; // 0: before cuts | 1: after cuts
  TH1F *h_mCentrality9[2];
  TH2F *h_mRefMultTofMatch[2];
  TH2F *h_mRefMultTofHits[2];
  TH2F *h_mVzVzVpd[2];
  TH1F *h_mDiffVzVzVpd[2];
  TH1F *h_mVertexZ[2];
  TH2F *h_mVertexXY[2];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    string HistName = Form("h_mRefMult_%s",CutsQA[i_cut].c_str());
    h_mRefMult[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mRefMult[i_cut]->SetLineColor(i_cut+1);

    HistName = Form("h_mCentrality9_%s",CutsQA[i_cut].c_str());
    h_mCentrality9[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mCentrality9[i_cut]->SetLineColor(i_cut+1);

    HistName = Form("h_mRefMultTofMatch_%s",CutsQA[i_cut].c_str());
    h_mRefMultTofMatch[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    // h_mRefMultTofMatch[i_cut]->SetLineColor(i_cut+1);

    HistName = Form("h_mRefMultTofHits_%s",CutsQA[i_cut].c_str());
    h_mRefMultTofHits[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    // h_mRefMultTofHits[i_cut]->SetLineColor(i_cut+1);

    HistName = Form("h_mVertexXY_%s",CutsQA[i_cut].c_str());
    h_mVertexXY[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    // h_mVertexXY[i_cut]->SetLineColor(i_cut+1);

    HistName = Form("h_mVertexZ_%s",CutsQA[i_cut].c_str());
    h_mVertexZ[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mVertexZ[i_cut]->SetLineColor(i_cut+1);

    HistName = Form("h_mVzVzVpd_%s",CutsQA[i_cut].c_str());
    h_mVzVzVpd[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
    // h_mVzVzVpd[i_cut]->SetLineColor(i_cut+1);

    HistName = Form("h_mDiffVzVzVpd_%s",CutsQA[i_cut].c_str());
    h_mDiffVzVzVpd[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
    h_mDiffVzVzVpd[i_cut]->SetLineColor(i_cut+1);
  }

  TF1 *f_tofCut_low = new TF1("f_tofCut_low","0.75*x-20",0,800);
  f_tofCut_low->SetLineColor(2);
  f_tofCut_low->SetLineWidth(2);
  f_tofCut_low->SetLineStyle(2);

  TF1 *f_tofCut_up  = new TF1("f_tofCut_up","1.80*x+15",0,800);
  f_tofCut_up->SetLineColor(2);
  f_tofCut_up->SetLineWidth(2);
  f_tofCut_up->SetLineStyle(2);

  TCanvas *c_EventQA[2];
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    string CanName = Form("c_EventQA_%s",CutsQA[i_cut].c_str());
    c_EventQA[i_cut] = new TCanvas(CanName.c_str(),CanName.c_str(),10,10,1600,800);
    c_EventQA[i_cut]->Divide(4,2);
    for(int i_pad = 0; i_pad < 8; ++i_pad)
    {
      c_EventQA[i_cut]->cd(i_pad+1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetLeftMargin(0.1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetRightMargin(0.1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetBottomMargin(0.1);
      c_EventQA[i_cut]->cd(i_pad+1)->SetGrid(0,0);
      c_EventQA[i_cut]->cd(i_pad+1)->SetTicks(1,1);
    }

    c_EventQA[i_cut]->cd(1);
    c_EventQA[i_cut]->cd(1)->SetLogy();
    h_mRefMult[i_cut]->Draw("hE");

    c_EventQA[i_cut]->cd(2);
    c_EventQA[i_cut]->cd(2)->SetLogz();
    h_mRefMultTofMatch[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mRefMultTofMatch[i_cut]->GetYaxis()->SetRangeUser(0.0,800.0);
    h_mRefMultTofMatch[i_cut]->Draw("colz");
    if(i_cut == 0)
    {
      f_tofCut_low->Draw("l same");
      f_tofCut_up->Draw("l same");
    }

    c_EventQA[i_cut]->cd(3);
    c_EventQA[i_cut]->cd(3)->SetLogz();
    h_mRefMultTofHits[i_cut]->GetXaxis()->SetRangeUser(0.0,800.0);
    h_mRefMultTofHits[i_cut]->Draw("colz");

    c_EventQA[i_cut]->cd(4);
    c_EventQA[i_cut]->cd(4)->SetLogy();
    h_mCentrality9[i_cut]->Draw("hE");

    c_EventQA[i_cut]->cd(5);
    c_EventQA[i_cut]->cd(5)->SetLogz();
    h_mVertexXY[i_cut]->Draw("colz");

    c_EventQA[i_cut]->cd(6);
    // c_EventQA[i_cut]->cd(6)->SetLogy();
    h_mVertexZ[i_cut]->Draw();

    c_EventQA[i_cut]->cd(7);
    c_EventQA[i_cut]->cd(7)->SetLogz();
    h_mVzVzVpd[i_cut]->Draw("colz");

    c_EventQA[i_cut]->cd(8);
    // c_EventQA[i_cut]->cd(8)->SetLogy();
    h_mDiffVzVzVpd[i_cut]->Draw();
  }
}
