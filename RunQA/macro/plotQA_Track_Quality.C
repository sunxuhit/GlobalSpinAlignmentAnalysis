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

void plotQA_Track_Quality(int energy = 2)
{
  // string JobId = "FAEC35EE309DEE5B0A7F0F9B57AA6008";
  string JobId = "B1C42134A640995406F8513F28A6447D";
  string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/QA/file_%s_QA_%s.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str(),JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mDca[2];
  TH1F *h_mNHitsFit[2];
  TH1F *h_mNHitsRatio[2];
  TH1F *h_mNHitsDEdx[2];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    std::string HistName = Form("h_mDca_%s",CutsQA[i_cut].c_str());
    h_mDca[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mNHitsFit_%s",CutsQA[i_cut].c_str());
    h_mNHitsFit[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mNHitsRatio_%s",CutsQA[i_cut].c_str());
    h_mNHitsRatio[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mNHitsDEdx_%s",CutsQA[i_cut].c_str());
    h_mNHitsDEdx[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
  }

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

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    c_TrackQA_Quality->cd(i_cut*4+1);
    h_mDca[i_cut]->GetXaxis()->SetTitle("dca");
    h_mDca[i_cut]->Draw("hE");

    c_TrackQA_Quality->cd(i_cut*4+2);
    h_mNHitsFit[i_cut]->Draw("hE");
    h_mNHitsFit[i_cut]->GetXaxis()->SetTitle("nHitsFit");

    c_TrackQA_Quality->cd(i_cut*4+3);
    h_mNHitsRatio[i_cut]->Draw("hE");
    h_mNHitsRatio[i_cut]->GetXaxis()->SetTitle("nHitsRatio");

    c_TrackQA_Quality->cd(i_cut*4+4);
    h_mNHitsDEdx[i_cut]->Draw("hE");
    h_mNHitsDEdx[i_cut]->GetXaxis()->SetTitle("nHitsDedx");
  }

  string FigName = Form("c_TrackQA_Quality_%s.eps",vmsa::mBeamEnergy[energy].c_str());
  c_TrackQA_Quality->SaveAs(FigName.c_str());
}
