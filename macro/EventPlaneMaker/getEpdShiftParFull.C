#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getEpdShiftParFull(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumShiftCorr = 20;
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_ShiftParFull_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // Shift Correction for Full EP
  TProfile2D *p_mEpdQ1ShiftCosFull[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mEpdQ1ShiftSinFull[mNumVzBin][mNumShiftCorr];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1ShiftCos%dFullVz%d",iShift,iVz);
      p_mEpdQ1ShiftCosFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mEpdQ1ShiftSin%dFullVz%d",iShift,iVz);
      p_mEpdQ1ShiftSinFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
    }
  }

  {
    TCanvas *c_EpdQ1ShiftFull = new TCanvas("c_EpdQ1ShiftFull","c_EpdQ1ShiftFull",10,10,800,400);
    c_EpdQ1ShiftFull->Divide(2,1);
    for(int iPad = 0; iPad < 2; ++iPad)
    {
      c_EpdQ1ShiftFull->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdQ1ShiftFull->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdQ1ShiftFull->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdQ1ShiftFull->cd(iPad+1)->SetTicks(1,1);
      c_EpdQ1ShiftFull->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/%s/EventPlaneMaker/EpdQ1ShiftFull_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdQ1ShiftFull->Print(figName.c_str());
    for(int iVz = 0; iVz < mNumVzBin; ++iVz)
    {
      for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
      {
	figName = Form("../../figures/%s/EventPlaneMaker/EpdQ1ShiftFull_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
	c_EpdQ1ShiftFull->cd(1)->Clear(); c_EpdQ1ShiftFull->cd(1); p_mEpdQ1ShiftCosFull[iVz][iShift]->Draw("colz");
	c_EpdQ1ShiftFull->cd(2)->Clear(); c_EpdQ1ShiftFull->cd(2); p_mEpdQ1ShiftSinFull[iVz][iShift]->Draw("colz");
	c_EpdQ1ShiftFull->Update();
	c_EpdQ1ShiftFull->Print(figName.c_str());
      }
    }
    figName = Form("../../figures/%s/EventPlaneMaker/EpdQ1ShiftFull_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdQ1ShiftFull->Print(figName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ShiftPar/file_EpdShiftParFull_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mEpdQ1ShiftCosFull[iVz][iShift]->Write();
      p_mEpdQ1ShiftSinFull[iVz][iShift]->Write();
    }
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mEpdEp1ShiftEast[mNumCentrality]; // 1st weighted EP
  TH2F *h_mEpdEp1ShiftWest[mNumCentrality];
  TH2F *h_mEpdEp1ShiftFull[mNumCentrality];
  TH2F *h_mEpdEp1ShiftCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1ShiftEastCent%d",iCent); // 1st EP
    h_mEpdEp1ShiftEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ShiftWestCent%d",iCent);
    h_mEpdEp1ShiftWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ShiftFullCent%d",iCent);
    h_mEpdEp1ShiftFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ShiftCorrCent%d",iCent);
    h_mEpdEp1ShiftCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  {
    TCanvas *c_EpdEp1ShiftDist = new TCanvas("c_EpdEp1ShiftDist","c_EpdEp1ShiftDist",10,10,800,800);
    c_EpdEp1ShiftDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_EpdEp1ShiftDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetTicks(1,1);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/%s/EventPlaneMaker/EpdShiftEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1ShiftDist->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      figName = Form("../../figures/%s/EventPlaneMaker/EpdShiftEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_EpdEp1ShiftDist->cd(1)->Clear(); c_EpdEp1ShiftDist->cd(1); h_mEpdEp1ShiftEast[iCent]->ProjectionY()->Draw();
      c_EpdEp1ShiftDist->cd(2)->Clear(); c_EpdEp1ShiftDist->cd(2); h_mEpdEp1ShiftWest[iCent]->ProjectionY()->Draw();
      c_EpdEp1ShiftDist->cd(3)->Clear(); c_EpdEp1ShiftDist->cd(3); h_mEpdEp1ShiftFull[iCent]->ProjectionY()->Draw();
      c_EpdEp1ShiftDist->cd(4)->Clear(); c_EpdEp1ShiftDist->cd(4); h_mEpdEp1ShiftCorr[iCent]->Draw("colz");
      c_EpdEp1ShiftDist->Update();
      c_EpdEp1ShiftDist->Print(figName.c_str());
    }
    figName = Form("../../figures/%s/EventPlaneMaker/EpdShiftEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1ShiftDist->Print(figName.c_str());
  }

  string outputFileShiftEp = Form("../../data/%s/EventPlaneMaker/file_EpdShiftEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileShiftEp.c_str() << endl;
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
