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

void plotRunQA_TrackKinematics(int beamType = 0)
{
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/%s/RunQA/file_%s_RunQA.root",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH1F *h_mPrimPt[numCuts][numTriggerBins];
  TH1F *h_mPrimEta[numCuts][numTriggerBins];
  TH1F *h_mPrimPhi[numCuts][numTriggerBins];
  TH1F *h_mGlobPt[numCuts][numTriggerBins];
  TH1F *h_mGlobEta[numCuts][numTriggerBins];
  TH1F *h_mGlobPhi[numCuts][numTriggerBins];

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < numTriggerBins; ++iTrig)
    {
      std::string HistName = Form("h_mPrimPt%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimPt[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mPrimEta%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimEta[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mPrimPhi%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mPrimPhi[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mGlobPt%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mGlobPt[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mGlobEta%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mGlobEta[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mGlobPhi%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mGlobPhi[iCut][iTrig] = (TH1F*)File_InPut->Get(HistName.c_str());
    }
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

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_TrackQA_prim->cd(iCut*3+1)->SetLogy();
    c_TrackQA_prim->cd(iCut*3+1);
    h_mPrimPt[iCut][9]->Draw("hE");

    c_TrackQA_prim->cd(iCut*3+2);
    h_mPrimEta[iCut][9]->Draw("hE");

    c_TrackQA_prim->cd(iCut*3+3);
    h_mPrimPhi[iCut][9]->Draw("hE");
  }
  string FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_TrackQA_primKinematics_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_TrackQA_prim->SaveAs(FigName.c_str());

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

  for(int iCut = 0; iCut < 2; ++iCut)
  {
    c_TrackQA_glob->cd(iCut*3+1)->SetLogy();
    h_mGlobPt[iCut][9]->Draw("hE");

    c_TrackQA_glob->cd(iCut*3+2);
    h_mGlobEta[iCut][9]->Draw("hE");

    c_TrackQA_glob->cd(iCut*3+3);
    h_mGlobPhi[iCut][9]->Draw("hE");
  }
  FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_TrackQA_globKinematics_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_TrackQA_glob->SaveAs(FigName.c_str());
}
