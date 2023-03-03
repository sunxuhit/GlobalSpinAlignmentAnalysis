#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getZdcShiftParFull(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumShiftCorr = 20;
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_ShiftParFull_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // Shift Correction for Full EP
  TProfile2D *p_mZdcQ1ShiftCosFull[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
  TProfile2D *p_mZdcQ1ShiftSinFull[mNumVzBin][mNumShiftCorr];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQ1ShiftCos%dFullVz%d",iShift,iVz);
      p_mZdcQ1ShiftCosFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mZdcQ1ShiftSin%dFullVz%d",iShift,iVz);
      p_mZdcQ1ShiftSinFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
    }
  }

  {
    TCanvas *c_ZdcQ1ShiftFull = new TCanvas("c_ZdcQ1ShiftFull","c_ZdcQ1ShiftFull",10,10,800,400);
    c_ZdcQ1ShiftFull->Divide(2,1);
    for(int iPad = 0; iPad < 2; ++iPad)
    {
      c_ZdcQ1ShiftFull->cd(iPad+1)->SetLeftMargin(0.15);
      c_ZdcQ1ShiftFull->cd(iPad+1)->SetRightMargin(0.15);
      c_ZdcQ1ShiftFull->cd(iPad+1)->SetBottomMargin(0.15);
      c_ZdcQ1ShiftFull->cd(iPad+1)->SetTicks(1,1);
      c_ZdcQ1ShiftFull->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/ZdcQ1ShiftFull_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcQ1ShiftFull->Print(figName.c_str());
    for(int iVz = 0; iVz < mNumVzBin; ++iVz)
    {
      for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
      {
	figName = Form("../../figures/EventPlaneMaker/%s/ZdcQ1ShiftFull_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
	c_ZdcQ1ShiftFull->cd(1)->Clear(); c_ZdcQ1ShiftFull->cd(1); p_mZdcQ1ShiftCosFull[iVz][iShift]->Draw("colz");
	c_ZdcQ1ShiftFull->cd(2)->Clear(); c_ZdcQ1ShiftFull->cd(2); p_mZdcQ1ShiftSinFull[iVz][iShift]->Draw("colz");
	c_ZdcQ1ShiftFull->Update();
	c_ZdcQ1ShiftFull->Print(figName.c_str());
      }
    }
    figName = Form("../../figures/EventPlaneMaker/%s/ZdcQ1ShiftFull_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcQ1ShiftFull->Print(figName.c_str());
  }


  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ShiftPar/file_ZdcShiftParFull_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mZdcQ1ShiftCosFull[iVz][iShift]->Write();
      p_mZdcQ1ShiftSinFull[iVz][iShift]->Write();
    }
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mZdcEp1ShiftEast[mNumCentrality]; // shift EP
  TH2F *h_mZdcEp1ShiftWest[mNumCentrality];
  TH2F *h_mZdcEp1ShiftFull[mNumCentrality]; // Qwest-QEast
  TH2F *h_mZdcEp1ShiftCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEp1ShiftEastCent%d",iCent);
    h_mZdcEp1ShiftEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEp1ShiftWestCent%d",iCent);
    h_mZdcEp1ShiftWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEp1ShiftFullCent%d",iCent);
    h_mZdcEp1ShiftFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEp1ShiftCorrCent%d",iCent);
    h_mZdcEp1ShiftCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  {
    TCanvas *c_ZdcEp1ShiftDist = new TCanvas("c_ZdcEp1ShiftDist","c_ZdcEp1ShiftDist",10,10,800,800);
    c_ZdcEp1ShiftDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_ZdcEp1ShiftDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_ZdcEp1ShiftDist->cd(iPad+1)->SetRightMargin(0.15);
      c_ZdcEp1ShiftDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_ZdcEp1ShiftDist->cd(iPad+1)->SetTicks(1,1);
      c_ZdcEp1ShiftDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/ZdcShiftEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcEp1ShiftDist->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      figName = Form("../../figures/EventPlaneMaker/%s/ZdcShiftEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_ZdcEp1ShiftDist->cd(1)->Clear(); c_ZdcEp1ShiftDist->cd(1); h_mZdcEp1ShiftEast[iCent]->ProjectionY()->Draw();
      c_ZdcEp1ShiftDist->cd(2)->Clear(); c_ZdcEp1ShiftDist->cd(2); h_mZdcEp1ShiftWest[iCent]->ProjectionY()->Draw();
      c_ZdcEp1ShiftDist->cd(3)->Clear(); c_ZdcEp1ShiftDist->cd(3); h_mZdcEp1ShiftFull[iCent]->ProjectionY()->Draw();
      c_ZdcEp1ShiftDist->cd(4)->Clear(); c_ZdcEp1ShiftDist->cd(4); h_mZdcEp1ShiftCorr[iCent]->Draw("colz");
      c_ZdcEp1ShiftDist->Update();
      c_ZdcEp1ShiftDist->Print(figName.c_str());
    }
    figName = Form("../../figures/EventPlaneMaker/%s/ZdcShiftEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcEp1ShiftDist->Print(figName.c_str());
  }

  string outputFileShiftEp = Form("../../data/EventPlaneMaker/%s/file_ZdcShiftEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Shift EP: " << outputFileShiftEp.c_str() << endl;
  TFile *file_OutPutShiftEp = new TFile(outputFileShiftEp.c_str(),"RECREATE");
  file_OutPutShiftEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEp1ShiftEast[iCent]->Write();
    h_mZdcEp1ShiftWest[iCent]->Write();
    h_mZdcEp1ShiftFull[iCent]->Write();
    h_mZdcEp1ShiftCorr[iCent]->Write();
  }
  file_OutPutShiftEp->Close();
}
