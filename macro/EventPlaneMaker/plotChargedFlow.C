#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void plotChargedFlow(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_ChargedFlow_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  TProfile *p_mZdcFullEpV1[mNumCentrality]; // v1 vs. eta
  TProfile *p_mEpdSubEpV1[mNumCentrality]; // v1 vs. eta
  // TProfile *p_mTpcSubEpV1[mNumCentrality]; // v1 vs. pT
  TProfile *p_mTpcSubEpV2[mNumCentrality]; // v2 vs. pT
  TProfile *p_mTpcSubEpV3[mNumCentrality]; // v3 vs. pT

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string proName = Form("p_mZdcFullEpV1Cent%d",iCent);
    p_mZdcFullEpV1[iCent] = (TProfile*)file_InPut->Get(proName.c_str());
    proName = Form("p_mEpdSubEpSideV1Cent%d",iCent);
    p_mEpdSubEpV1[iCent] = (TProfile*)file_InPut->Get(proName.c_str());
    // proName = Form("p_mTpcSubEpV1Cent%d",iCent);
    // p_mTpcSubEpV1[iCent] = (TProfile*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcSubEpV2Cent%d",iCent);
    p_mTpcSubEpV2[iCent] = (TProfile*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcSubEpV3Cent%d",iCent);
    p_mTpcSubEpV3[iCent] = (TProfile*)file_InPut->Get(proName.c_str());
  }

  TCanvas *c_ChargedFlow = new TCanvas("c_ChargedFlow","c_ChargedFlow",10,10,900,900);
  c_ChargedFlow->Divide(3,3);
  for(int iPad = 0; iPad < 9; ++iPad)
  {
    c_ChargedFlow->cd(iPad+1)->SetLeftMargin(0.15);
    c_ChargedFlow->cd(iPad+1)->SetRightMargin(0.15);
    c_ChargedFlow->cd(iPad+1)->SetBottomMargin(0.15);
    c_ChargedFlow->cd(iPad+1)->SetTicks(1,1);
    c_ChargedFlow->cd(iPad+1)->SetGrid(0,0);
  }

  std::string figName = Form("../../figures/EventPlaneMaker/%s/ChargedFlow_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_ChargedFlow->Print(figName.c_str());

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  { // plot v1 vs. eta from ZDC & EPD
    c_ChargedFlow->cd(iCent+1);
    p_mEpdSubEpV1[iCent]->SetMarkerStyle(24);
    p_mEpdSubEpV1[iCent]->SetMarkerSize(1.0);
    p_mEpdSubEpV1[iCent]->SetMarkerColor(2);
    p_mEpdSubEpV1[iCent]->GetYaxis()->SetNdivisions(505);
    p_mEpdSubEpV1[iCent]->GetYaxis()->SetRangeUser(-0.04,0.04);
    p_mEpdSubEpV1[iCent]->GetYaxis()->SetTitle("v_{1}");
    p_mEpdSubEpV1[iCent]->GetYaxis()->SetTitleSize(0.06);
    p_mEpdSubEpV1[iCent]->GetYaxis()->CenterTitle();
    p_mEpdSubEpV1[iCent]->GetXaxis()->SetTitle("#eta");
    p_mEpdSubEpV1[iCent]->GetXaxis()->SetTitleSize(0.06);
    p_mEpdSubEpV1[iCent]->Draw("pE");
    p_mZdcFullEpV1[iCent]->SetMarkerStyle(25);
    p_mZdcFullEpV1[iCent]->SetMarkerSize(1.0);
    p_mZdcFullEpV1[iCent]->SetMarkerColor(4);
    p_mZdcFullEpV1[iCent]->Draw("pE same");

    TLegend *leg = new TLegend(0.2,0.7,0.5,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(p_mEpdSubEpV1[iCent],"EPD","P");
    leg->AddEntry(p_mZdcFullEpV1[iCent],"ZDC","P");
    leg->Draw("same");
  }
  figName = Form("../../figures/EventPlaneMaker/%s/ChargedFlow_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_ChargedFlow->Update();
  c_ChargedFlow->Print(figName.c_str());

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  { // plot v2 & v3 vs. pT from TPC
    c_ChargedFlow->cd(iCent+1);
    p_mTpcSubEpV2[iCent]->SetMarkerStyle(24);
    p_mTpcSubEpV2[iCent]->SetMarkerSize(1.0);
    p_mTpcSubEpV2[iCent]->SetMarkerColor(2);
    p_mTpcSubEpV2[iCent]->GetYaxis()->SetNdivisions(505);
    p_mTpcSubEpV2[iCent]->GetYaxis()->SetRangeUser(-0.01,0.35);
    p_mTpcSubEpV2[iCent]->GetYaxis()->SetTitle("flow");
    p_mTpcSubEpV2[iCent]->GetYaxis()->SetTitleSize(0.06);
    p_mTpcSubEpV2[iCent]->GetYaxis()->CenterTitle();
    p_mTpcSubEpV2[iCent]->GetXaxis()->SetTitle("p_{T}");
    p_mTpcSubEpV2[iCent]->GetXaxis()->SetTitleSize(0.06);
    p_mTpcSubEpV2[iCent]->Draw("pE");
    p_mTpcSubEpV3[iCent]->SetMarkerStyle(25);
    p_mTpcSubEpV3[iCent]->SetMarkerSize(1.0);
    p_mTpcSubEpV3[iCent]->SetMarkerColor(4);
    p_mTpcSubEpV3[iCent]->Draw("pE same");

    TLegend *leg = new TLegend(0.2,0.7,0.5,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(p_mTpcSubEpV2[iCent],"TPC v_{2}","P");
    leg->AddEntry(p_mTpcSubEpV3[iCent],"TPC v_{3}","P");
    leg->Draw("same");
  }
  figName = Form("../../figures/EventPlaneMaker/%s/ChargedFlow_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_ChargedFlow->Update();
  c_ChargedFlow->Print(figName.c_str());


  figName = Form("../../figures/EventPlaneMaker/%s/ChargedFlow_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_ChargedFlow->Print(figName.c_str());
}
