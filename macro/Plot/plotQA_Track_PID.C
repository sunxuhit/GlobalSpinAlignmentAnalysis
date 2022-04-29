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

void plotQA_Track_PID(int energy = 2)
{
  string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/QA/file_%s_QA.root",vmsa::mBeamEnergy[energy].c_str(),vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH2F *h_mDEdxMom[2];
  TH2F *h_mBetaMom[2];
  TH2F *h_mMass2Mom[2];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    std::string HistName = Form("h_mDEdxMom_%s",CutsQA[i_cut].c_str());
    h_mDEdxMom[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mMass2Mom_%s",CutsQA[i_cut].c_str());
    h_mMass2Mom[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());

    HistName = Form("h_mBetaMom_%s",CutsQA[i_cut].c_str());
    h_mBetaMom[i_cut] = (TH2F*)File_InPut->Get(HistName.c_str());
  }

  TCanvas *c_TrackQA_PID = new TCanvas("c_TrackQA_PID","c_TrackQA_PID",10,10,1200,800);
  c_TrackQA_PID->Divide(3,2);
  for(int i_pad = 0; i_pad < 6; ++i_pad)
  {
    c_TrackQA_PID->cd(i_pad+1);
    c_TrackQA_PID->cd(i_pad+1)->SetLeftMargin(0.1);
    c_TrackQA_PID->cd(i_pad+1)->SetRightMargin(0.1);
    c_TrackQA_PID->cd(i_pad+1)->SetBottomMargin(0.1);
    c_TrackQA_PID->cd(i_pad+1)->SetGrid(0,0);
    c_TrackQA_PID->cd(i_pad+1)->SetTicks(1,1);
    c_TrackQA_PID->cd(i_pad+1)->SetLogz();
  }

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    c_TrackQA_PID->cd(i_cut*3+1);
    h_mDEdxMom[i_cut]->Draw("colz");

    c_TrackQA_PID->cd(i_cut*3+2);
    h_mMass2Mom[i_cut]->Draw("colz");

    c_TrackQA_PID->cd(i_cut*3+3);
    h_mBetaMom[i_cut]->Draw("colz");
  }
}
