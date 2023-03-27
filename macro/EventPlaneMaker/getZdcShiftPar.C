#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getZdcShiftPar(int beamType = 2)
{
  gStyle->SetOptStat(0);
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumShiftCorr = 20;
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_ShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // Shift Correction for East/West EP
  TProfile2D *p_mZdcQ1ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 0 = vertex neg/pos | 1 = shift correction harmonics
  TProfile2D *p_mZdcQ1ShiftSinEast[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mZdcQ1ShiftCosWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mZdcQ1ShiftSinWest[mNumVzBin][mNumShiftCorr];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQ1ShiftCos%dEastVz%d",iShift,iVz);
      p_mZdcQ1ShiftCosEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mZdcQ1ShiftSin%dEastVz%d",iShift,iVz);
      p_mZdcQ1ShiftSinEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mZdcQ1ShiftCos%dWestVz%d",iShift,iVz);
      p_mZdcQ1ShiftCosWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mZdcQ1ShiftSin%dWestVz%d",iShift,iVz);
      p_mZdcQ1ShiftSinWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
    }
  }

  {
    TCanvas *c_ZdcQ1Shift = new TCanvas("c_ZdcQ1Shift","c_ZdcQ1Shift",10,10,800,800);
    c_ZdcQ1Shift->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_ZdcQ1Shift->cd(iPad+1)->SetLeftMargin(0.15);
      c_ZdcQ1Shift->cd(iPad+1)->SetRightMargin(0.15);
      c_ZdcQ1Shift->cd(iPad+1)->SetBottomMargin(0.15);
      c_ZdcQ1Shift->cd(iPad+1)->SetTicks(1,1);
      c_ZdcQ1Shift->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/ZdcQ1Shift_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcQ1Shift->Print(figName.c_str());
    for(int iVz = 0; iVz < mNumVzBin; ++iVz)
    {
      for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
      {
	figName = Form("../../figures/EventPlaneMaker/%s/ZdcQ1Shift_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
	c_ZdcQ1Shift->cd(1)->Clear(); c_ZdcQ1Shift->cd(1); p_mZdcQ1ShiftCosEast[iVz][iShift]->Draw("colz");
	c_ZdcQ1Shift->cd(2)->Clear(); c_ZdcQ1Shift->cd(2); p_mZdcQ1ShiftSinEast[iVz][iShift]->Draw("colz");
	c_ZdcQ1Shift->cd(3)->Clear(); c_ZdcQ1Shift->cd(3); p_mZdcQ1ShiftCosWest[iVz][iShift]->Draw("colz");
	c_ZdcQ1Shift->cd(4)->Clear(); c_ZdcQ1Shift->cd(4); p_mZdcQ1ShiftSinWest[iVz][iShift]->Draw("colz");
	c_ZdcQ1Shift->Update();
	c_ZdcQ1Shift->Print(figName.c_str());
      }
    }
    figName = Form("../../figures/EventPlaneMaker/%s/ZdcQ1Shift_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcQ1Shift->Print(figName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ShiftPar/file_ZdcShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mZdcQ1ShiftCosEast[iVz][iShift]->Write();
      p_mZdcQ1ShiftSinEast[iVz][iShift]->Write();
      p_mZdcQ1ShiftCosWest[iVz][iShift]->Write();
      p_mZdcQ1ShiftSinWest[iVz][iShift]->Write();
    }
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mZdcEp1ReCtrEast[mNumCentrality]; // recenter EP
  TH2F *h_mZdcEp1ReCtrWest[mNumCentrality];
  TH2F *h_mZdcEp1ReCtrFull[mNumCentrality]; // Qwest-QEast
  TH2F *h_mZdcEp1ReCtrCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEp1ReCtrEastCent%d",iCent);
    h_mZdcEp1ReCtrEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEp1ReCtrWestCent%d",iCent);
    h_mZdcEp1ReCtrWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mZdcEp1ReCtrFullCent%d",iCent);
    h_mZdcEp1ReCtrFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEp1ReCtrCorrCent%d",iCent);
    h_mZdcEp1ReCtrCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  {
    TCanvas *c_ZdcEp1ReCtrDist = new TCanvas("c_ZdcEp1ReCtrDist","c_ZdcEp1ReCtrDist",10,10,800,800);
    c_ZdcEp1ReCtrDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_ZdcEp1ReCtrDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_ZdcEp1ReCtrDist->cd(iPad+1)->SetRightMargin(0.15);
      c_ZdcEp1ReCtrDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_ZdcEp1ReCtrDist->cd(iPad+1)->SetTicks(1,1);
      c_ZdcEp1ReCtrDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/ZdcReCtrEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcEp1ReCtrDist->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      figName = Form("../../figures/EventPlaneMaker/%s/ZdcReCtrEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_ZdcEp1ReCtrDist->cd(1)->Clear(); c_ZdcEp1ReCtrDist->cd(1); h_mZdcEp1ReCtrEast[iCent]->ProjectionY()->Draw();
      c_ZdcEp1ReCtrDist->cd(2)->Clear(); c_ZdcEp1ReCtrDist->cd(2); h_mZdcEp1ReCtrWest[iCent]->ProjectionY()->Draw();
      c_ZdcEp1ReCtrDist->cd(3)->Clear(); c_ZdcEp1ReCtrDist->cd(3); h_mZdcEp1ReCtrFull[iCent]->ProjectionY()->Draw();
      c_ZdcEp1ReCtrDist->cd(4)->Clear(); c_ZdcEp1ReCtrDist->cd(4); h_mZdcEp1ReCtrCorr[iCent]->Draw("colz");
      c_ZdcEp1ReCtrDist->Update();
      c_ZdcEp1ReCtrDist->Print(figName.c_str());
    }
    figName = Form("../../figures/EventPlaneMaker/%s/ZdcReCtrEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcEp1ReCtrDist->Print(figName.c_str());
  }

  string outputFileReCtrEp = Form("../../data/EventPlaneMaker/%s/file_ZdcReCtrEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of ReCtr EP: " << outputFileReCtrEp.c_str() << endl;
  TFile *file_OutPutReCtrEp = new TFile(outputFileReCtrEp.c_str(),"RECREATE");
  file_OutPutReCtrEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEp1ReCtrEast[iCent]->Write();
    h_mZdcEp1ReCtrWest[iCent]->Write();
    h_mZdcEp1ReCtrFull[iCent]->Write();
    h_mZdcEp1ReCtrCorr[iCent]->Write();
  }
  file_OutPutReCtrEp->Close();
}
