#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "../../Utility/include/StSpinAlignmentCons.h"

using namespace std;

static const string CutStatus[2] = {"Bf","Af"};
static const int numCuts = 2; // 0: before cuts | 1: after cuts
static const int numTriggerBins = 10; // 0-8 for different triggerID | 9 for all triggers

void plotRunQA_TrackQuality(int beamType = 0)
{
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/%s/RunQA/file_%s_RunQA.root",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mDca[numCuts][numTriggerBins];
  TH1F *h_mNHitsFit[numCuts][numTriggerBins];
  TH1F *h_mNHitsRatio[numCuts][numTriggerBins];
  TH1F *h_mNHitsDEdx[numCuts][numTriggerBins];

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < numTriggerBins; ++iTrig)
    {
      string HistName = Form("h_mDca%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mDca[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mNHitsFit%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mNHitsFit[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mNHitsRatio%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mNHitsRatio[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mNHitsDEdx%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mNHitsDEdx[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());
    }
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

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_TrackQA_Quality->cd(iCut*4+1)->SetLogy();
    h_mDca[iCut][9]->GetXaxis()->SetTitle("dca");
    h_mDca[iCut][9]->Draw("hE");

    c_TrackQA_Quality->cd(iCut*4+2);
    h_mNHitsFit[iCut][9]->Draw("hE");
    h_mNHitsFit[iCut][9]->GetXaxis()->SetTitle("nHitsFit");

    c_TrackQA_Quality->cd(iCut*4+3);
    h_mNHitsRatio[iCut][9]->Draw("hE");
    h_mNHitsRatio[iCut][9]->GetXaxis()->SetTitle("nHitsRatio");

    c_TrackQA_Quality->cd(iCut*4+4);
    h_mNHitsDEdx[iCut][9]->Draw("hE");
    h_mNHitsDEdx[iCut][9]->GetXaxis()->SetTitle("nHitsDedx");
  }

  string FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_TrackQA_Quality_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_TrackQA_Quality->SaveAs(FigName.c_str());
}
