#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getTpcResolution(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumShiftCorr = 20;
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_EpResolution_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  TProfile *p_mTpcSubEp2Res = (TProfile*)file_InPut->Get("p_mTpcSubEp2Res");
  TProfile *p_mTpcSubEp3Res = (TProfile*)file_InPut->Get("p_mTpcSubEp3Res");

  {
    TCanvas *c_TpcSubEp2Res = new TCanvas("c_TpcSubEp2Res","c_TpcSubEp2Res",10,10,800,800);
    c_TpcSubEp2Res->cd()->SetLeftMargin(0.15);
    c_TpcSubEp2Res->cd()->SetRightMargin(0.15);
    c_TpcSubEp2Res->cd()->SetBottomMargin(0.15);
    c_TpcSubEp2Res->cd()->SetTicks(1,1);
    c_TpcSubEp2Res->cd()->SetGrid(0,0);
    p_mTpcSubEp2Res->GetXaxis()->SetTitle("Centrality 9 Bins");
    p_mTpcSubEp2Res->GetYaxis()->SetTitle("cos(#Psi_{2}^{West}-#Psi_{2}^{East})");
    p_mTpcSubEp2Res->GetYaxis()->SetTitleOffset(1.25);
    p_mTpcSubEp2Res->GetYaxis()->SetNdivisions(505);
    p_mTpcSubEp2Res->GetYaxis()->SetLabelSize(0.03);
    p_mTpcSubEp2Res->Draw("PE");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/TpcSubEp2Resolution_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_TpcSubEp2Res->SaveAs(figName.c_str());
  }

  {
    TCanvas *c_TpcSubEp3Res = new TCanvas("c_TpcSubEp3Res","c_TpcSubEp3Res",10,10,800,800);
    c_TpcSubEp3Res->cd()->SetLeftMargin(0.15);
    c_TpcSubEp3Res->cd()->SetRightMargin(0.15);
    c_TpcSubEp3Res->cd()->SetBottomMargin(0.15);
    c_TpcSubEp3Res->cd()->SetTicks(1,1);
    c_TpcSubEp3Res->cd()->SetGrid(0,0);
    p_mTpcSubEp3Res->GetXaxis()->SetTitle("Centrality 9 Bins");
    p_mTpcSubEp3Res->GetYaxis()->SetTitle("cos(#Psi_{3}^{West}-#Psi_{3}^{East})");
    p_mTpcSubEp3Res->GetYaxis()->SetTitleOffset(1.25);
    p_mTpcSubEp3Res->GetYaxis()->SetNdivisions(505);
    p_mTpcSubEp3Res->GetYaxis()->SetLabelSize(0.03);
    p_mTpcSubEp3Res->Draw("PE");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/TpcSubEp3Resolution_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_TpcSubEp3Res->SaveAs(figName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/Resolution/file_TpcResolution_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  p_mTpcSubEp2Res->Write();
  p_mTpcSubEp3Res->Write();
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mTpcEp2ShiftEast[mNumCentrality]; // 2nd raw EP
  TH2F *h_mTpcEp2ShiftWest[mNumCentrality];
  TH2F *h_mTpcEp2ShiftFull[mNumCentrality];
  TH2F *h_mTpcEp2ShiftCorr[mNumCentrality]; // Psi2East vs Psi2West

  TH2F *h_mTpcEp3ShiftEast[mNumCentrality]; // 3rd raw EP
  TH2F *h_mTpcEp3ShiftWest[mNumCentrality];
  TH2F *h_mTpcEp3ShiftFull[mNumCentrality];
  TH2F *h_mTpcEp3ShiftCorr[mNumCentrality]; // Psi3East vs Psi3West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp2ShiftEastCent%d",iCent); // 2nd EP
    h_mTpcEp2ShiftEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2ShiftWestCent%d",iCent);
    h_mTpcEp2ShiftWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2ShiftFullCent%d",iCent);
    h_mTpcEp2ShiftFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2ShiftCorrCent%d",iCent);
    h_mTpcEp2ShiftCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mTpcEp3ShiftEastCent%d",iCent); // 3rd EP
    h_mTpcEp3ShiftEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3ShiftWestCent%d",iCent);
    h_mTpcEp3ShiftWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3ShiftFullCent%d",iCent);
    h_mTpcEp3ShiftFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3ShiftCorrCent%d",iCent);
    h_mTpcEp3ShiftCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  {
    TCanvas *c_TpcEpShiftDist = new TCanvas("c_TpcEpShiftDist","c_TpcEpShiftDist",10,10,800,800);
    c_TpcEpShiftDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_TpcEpShiftDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_TpcEpShiftDist->cd(iPad+1)->SetRightMargin(0.15);
      c_TpcEpShiftDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_TpcEpShiftDist->cd(iPad+1)->SetTicks(1,1);
      c_TpcEpShiftDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/%s/EventPlaneMaker/TpcShiftEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_TpcEpShiftDist->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      figName = Form("../../figures/%s/EventPlaneMaker/TpcShiftEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_TpcEpShiftDist->cd(1); h_mTpcEp2ShiftEast[iCent]->ProjectionY()->Draw();
      c_TpcEpShiftDist->cd(2); h_mTpcEp2ShiftWest[iCent]->ProjectionY()->Draw();
      c_TpcEpShiftDist->cd(3); h_mTpcEp2ShiftFull[iCent]->ProjectionY()->Draw();
      c_TpcEpShiftDist->cd(4); h_mTpcEp2ShiftCorr[iCent]->Draw("colz");
      c_TpcEpShiftDist->Update();
      c_TpcEpShiftDist->Print(figName.c_str());

      c_TpcEpShiftDist->cd(1); h_mTpcEp3ShiftEast[iCent]->ProjectionY()->Draw();
      c_TpcEpShiftDist->cd(2); h_mTpcEp3ShiftWest[iCent]->ProjectionY()->Draw();
      c_TpcEpShiftDist->cd(3); h_mTpcEp3ShiftFull[iCent]->ProjectionY()->Draw();
      c_TpcEpShiftDist->cd(4); h_mTpcEp3ShiftCorr[iCent]->Draw("colz");
      c_TpcEpShiftDist->Update();
      c_TpcEpShiftDist->Print(figName.c_str());
    }
    figName = Form("../../figures/%s/EventPlaneMaker/TpcShiftEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_TpcEpShiftDist->Print(figName.c_str());
  }

  string outputFileShiftEp = Form("../../data/%s/EventPlaneMaker/file_TpcShiftEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Shift EP: " << outputFileShiftEp.c_str() << endl;
  TFile *file_OutPutShiftEp = new TFile(outputFileShiftEp.c_str(),"RECREATE");
  file_OutPutShiftEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2ShiftEast[iCent]->Write();
    h_mTpcEp2ShiftWest[iCent]->Write();
    h_mTpcEp2ShiftFull[iCent]->Write();
    h_mTpcEp2ShiftCorr[iCent]->Write();

    h_mTpcEp3ShiftEast[iCent]->Write();
    h_mTpcEp3ShiftWest[iCent]->Write();
    h_mTpcEp3ShiftFull[iCent]->Write();
    h_mTpcEp3ShiftCorr[iCent]->Write();
  }
  file_OutPutShiftEp->Close();
}
