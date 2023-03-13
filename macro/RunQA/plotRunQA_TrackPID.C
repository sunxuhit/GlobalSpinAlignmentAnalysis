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

void plotRunQA_TrackPID(int beamType = 0)
{
  string inputfile = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/data/%s/RunQA/file_%s_RunQA.root",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());
  TH2F *h_mDEdxMom[numCuts][numTriggerBins];
  TH2F *h_mMass2Mom[numCuts][numTriggerBins];
  TH2F *h_mBetaMom[numCuts][numTriggerBins];

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < numTriggerBins; ++iTrig)
    {
      string HistName = Form("h_mDEdxMom%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mDEdxMom[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mMass2Mom%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mMass2Mom[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());

      HistName = Form("h_mBetaMom%sTrigger%d",CutStatus[iCut].c_str(),iTrig);
      h_mBetaMom[iCut][iTrig] = (TH2F*)File_InPut->Get(HistName.c_str());
    }
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

  for(int iCut = 0; iCut < numCuts; ++iCut)
  {
    c_TrackQA_PID->cd(iCut*3+1);
    h_mDEdxMom[iCut][9]->Draw("colz");

    c_TrackQA_PID->cd(iCut*3+2);
    h_mMass2Mom[iCut][9]->Draw("colz");

    c_TrackQA_PID->cd(iCut*3+3);
    h_mBetaMom[iCut][9]->Draw("colz");
  }

  string FigName = Form("/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/figures/%s/RunQA/c_TrackQA_PID_%s.pdf",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  c_TrackQA_PID->SaveAs(FigName.c_str());
}
