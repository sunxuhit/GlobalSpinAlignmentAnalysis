#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getEpdResolution(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumShiftCorr = 20;
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_EpResolution_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // Shift Correction for Full EP
  TProfile *p_mEpdSubEp1Res = (TProfile*)file_InPut->Get("p_mEpdSubEp1Res");;

  {
    TCanvas *c_EpdSubEp1Res = new TCanvas("c_EpdSubEp1Res","c_EpdSubEp1Res",10,10,800,800);
    c_EpdSubEp1Res->cd()->SetLeftMargin(0.15);
    c_EpdSubEp1Res->cd()->SetRightMargin(0.15);
    c_EpdSubEp1Res->cd()->SetBottomMargin(0.15);
    c_EpdSubEp1Res->cd()->SetTicks(1,1);
    c_EpdSubEp1Res->cd()->SetGrid(0,0);
    p_mEpdSubEp1Res->GetXaxis()->SetTitle("Centrality 9 Bins");
    p_mEpdSubEp1Res->GetYaxis()->SetTitle("cos(#Psi_{1}^{West}-#Psi_{1}^{East})");
    p_mEpdSubEp1Res->GetYaxis()->SetTitleOffset(1.25);
    p_mEpdSubEp1Res->GetYaxis()->SetNdivisions(505);
    p_mEpdSubEp1Res->GetYaxis()->SetLabelSize(0.03);
    p_mEpdSubEp1Res->Draw("PE");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/EpdSubEp1Resolution_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdSubEp1Res->SaveAs(figName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ShiftPar/file_EpdResolutoin_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  p_mEpdSubEp1Res->Write();
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mEpdEp1ShiftEast[mNumCentrality]; // shift EP
  TH2F *h_mEpdEp1ShiftWest[mNumCentrality];
  TH2F *h_mEpdEp1ShiftFull[mNumCentrality]; // Qwest-QEast
  TH2F *h_mEpdEp1ShiftFullCorr[mNumCentrality];
  TH2F *h_mEpdEp1ShiftCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1ShiftEastCent%d",iCent);
    h_mEpdEp1ShiftEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ShiftWestCent%d",iCent);
    h_mEpdEp1ShiftWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ShiftFullCent%d",iCent);
    h_mEpdEp1ShiftFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ShiftFullCorrCent%d",iCent);
    h_mEpdEp1ShiftFullCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ShiftCorrCent%d",iCent);
    h_mEpdEp1ShiftCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  {
    TCanvas *c_EpdEp1ShiftDist = new TCanvas("c_EpdEp1ShiftDist","c_EpdEp1ShiftDist",10,10,800,1200);
    c_EpdEp1ShiftDist->Divide(2,3);
    for(int iPad = 0; iPad < 6; ++iPad)
    {
      c_EpdEp1ShiftDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetTicks(1,1);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/%s/EventPlaneMaker/EpdShiftFullEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1ShiftDist->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      figName = Form("../../figures/%s/EventPlaneMaker/EpdShiftFullEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_EpdEp1ShiftDist->cd(1)->Clear(); c_EpdEp1ShiftDist->cd(1); h_mEpdEp1ShiftEast[iCent]->ProjectionY()->Draw();
      c_EpdEp1ShiftDist->cd(2)->Clear(); c_EpdEp1ShiftDist->cd(2); h_mEpdEp1ShiftWest[iCent]->ProjectionY()->Draw();
      c_EpdEp1ShiftDist->cd(3)->Clear(); c_EpdEp1ShiftDist->cd(3); h_mEpdEp1ShiftFull[iCent]->ProjectionY()->Draw();
      c_EpdEp1ShiftDist->cd(4)->Clear(); c_EpdEp1ShiftDist->cd(4); h_mEpdEp1ShiftFullCorr[iCent]->ProjectionY()->Draw();
      c_EpdEp1ShiftDist->cd(5)->Clear(); c_EpdEp1ShiftDist->cd(5); h_mEpdEp1ShiftCorr[iCent]->Draw("colz");
      c_EpdEp1ShiftDist->Update();
      c_EpdEp1ShiftDist->Print(figName.c_str());
    }
    figName = Form("../../figures/%s/EventPlaneMaker/EpdShiftFullEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1ShiftDist->Print(figName.c_str());
  }

  string outputFileShiftEp = Form("../../data/%s/EventPlaneMaker/file_EpdShiftFullEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Shift EP: " << outputFileShiftEp.c_str() << endl;
  TFile *file_OutPutShiftEp = new TFile(outputFileShiftEp.c_str(),"RECREATE");
  file_OutPutShiftEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1ShiftEast[iCent]->Write();
    h_mEpdEp1ShiftWest[iCent]->Write();
    h_mEpdEp1ShiftFull[iCent]->Write();
    h_mEpdEp1ShiftCorr[iCent]->Write();
  }
  file_OutPutShiftEp->Close();
}
