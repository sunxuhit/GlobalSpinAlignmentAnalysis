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

void plotQA_Track_Kinematics(int energy = 2)
{
  string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/QA/file_%s_QA.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mPrimPt[2]; // 0: before cuts | 1: after cuts
  TH1F *h_mPrimEta[2];
  TH1F *h_mPrimPhi[2];
  TH1F *h_mGlobPt[2];
  TH1F *h_mGlobEta[2];
  TH1F *h_mGlobPhi[2];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    std::string HistName = Form("h_mPrimPt_%s",CutsQA[i_cut].c_str());
    h_mPrimPt[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mPrimEta_%s",CutsQA[i_cut].c_str());
    h_mPrimEta[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mPrimPhi_%s",CutsQA[i_cut].c_str());
    h_mPrimPhi[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mGlobPt_%s",CutsQA[i_cut].c_str());
    h_mGlobPt[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mGlobEta_%s",CutsQA[i_cut].c_str());
    h_mGlobEta[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mGlobPhi_%s",CutsQA[i_cut].c_str());
    h_mGlobPhi[i_cut] = (TH1F*)File_InPut->Get(HistName.c_str());
  }

  TCanvas *c_TrackQA_prim = new TCanvas("c_TrackQA_prim","c_TrackQA_prim",10,10,1200,800);
  c_TrackQA_prim->Divide(3,2);
  for(int i_pad = 0; i_pad < 6; ++i_pad)
  {
    c_TrackQA_prim->cd(i_pad+1);
    c_TrackQA_prim->cd(i_pad+1)->SetLeftMargin(0.1);
    c_TrackQA_prim->cd(i_pad+1)->SetRightMargin(0.1);
    c_TrackQA_prim->cd(i_pad+1)->SetBottomMargin(0.1);
    c_TrackQA_prim->cd(i_pad+1)->SetGrid(0,0);
    c_TrackQA_prim->cd(i_pad+1)->SetTicks(1,1);
  }

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    c_TrackQA_prim->cd(i_cut*3+1);
    h_mPrimPt[i_cut]->Draw("hE");

    c_TrackQA_prim->cd(i_cut*3+2);
    h_mPrimEta[i_cut]->Draw("hE");

    c_TrackQA_prim->cd(i_cut*3+3);
    h_mPrimPhi[i_cut]->Draw("hE");
  }

  TCanvas *c_TrackQA_glob = new TCanvas("c_TrackQA_glob","c_TrackQA_glob",10,10,1200,800);
  c_TrackQA_glob->Divide(3,2);
  for(int i_pad = 0; i_pad < 6; ++i_pad)
  {
    c_TrackQA_glob->cd(i_pad+1);
    c_TrackQA_glob->cd(i_pad+1)->SetLeftMargin(0.1);
    c_TrackQA_glob->cd(i_pad+1)->SetRightMargin(0.1);
    c_TrackQA_glob->cd(i_pad+1)->SetBottomMargin(0.1);
    c_TrackQA_glob->cd(i_pad+1)->SetGrid(0,0);
    c_TrackQA_glob->cd(i_pad+1)->SetTicks(1,1);
  }

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    c_TrackQA_glob->cd(i_cut*3+1);
    h_mGlobPt[i_cut]->Draw("hE");

    c_TrackQA_glob->cd(i_cut*3+2);
    h_mGlobEta[i_cut]->Draw("hE");

    c_TrackQA_glob->cd(i_cut*3+3);
    h_mGlobPhi[i_cut]->Draw("hE");
  }
}
