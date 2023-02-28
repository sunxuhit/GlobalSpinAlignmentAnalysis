#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getEpdShiftPar(int beamType = 0)
{
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumShiftCorr = 20;
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_ShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
  TProfile2D *p_mEpdQ1ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 1st EP
  TProfile2D *p_mEpdQ1ShiftSinEast[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mEpdQ1ShiftCosWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mEpdQ1ShiftSinWest[mNumVzBin][mNumShiftCorr];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mEpdQ1ShiftCos%dEastVz%d",iShift,iVz); // 1st EP
      p_mEpdQ1ShiftCosEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mEpdQ1ShiftSin%dEastVz%d",iShift,iVz);
      p_mEpdQ1ShiftSinEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mEpdQ1ShiftCos%dWestVz%d",iShift,iVz);
      p_mEpdQ1ShiftCosWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mEpdQ1ShiftSin%dWestVz%d",iShift,iVz);
      p_mEpdQ1ShiftSinWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
    }
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ShiftPar/file_EpdShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mEpdQ1ShiftCosEast[iVz][iShift]->Write();
      p_mEpdQ1ShiftSinEast[iVz][iShift]->Write();
      p_mEpdQ1ShiftCosWest[iVz][iShift]->Write();
      p_mEpdQ1ShiftSinWest[iVz][iShift]->Write();
    }
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mEpdEp1ReCtrEast[mNumCentrality]; // 1st weighted EP
  TH2F *h_mEpdEp1ReCtrWest[mNumCentrality];
  TH2F *h_mEpdEp1ReCtrFull[mNumCentrality];
  TH2F *h_mEpdEp1ReCtrCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1ReCtrEastCent%d",iCent); // 1st EP
    h_mEpdEp1ReCtrEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ReCtrWestCent%d",iCent);
    h_mEpdEp1ReCtrWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ReCtrFullCent%d",iCent);
    h_mEpdEp1ReCtrFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1ReCtrCorrCent%d",iCent);
    h_mEpdEp1ReCtrCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  TCanvas *c_EpdEp1ReCtrDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_EpdEp1ReCtrDistCent%d",iCent);
    c_EpdEp1ReCtrDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1600,400);
    c_EpdEp1ReCtrDist[iCent]->Divide(4,1);
    c_EpdEp1ReCtrDist[iCent]->cd(1); h_mEpdEp1ReCtrEast[iCent]->ProjectionY()->Draw();
    c_EpdEp1ReCtrDist[iCent]->cd(2); h_mEpdEp1ReCtrWest[iCent]->ProjectionY()->Draw();
    c_EpdEp1ReCtrDist[iCent]->cd(3); h_mEpdEp1ReCtrFull[iCent]->ProjectionY()->Draw();
    c_EpdEp1ReCtrDist[iCent]->cd(4); h_mEpdEp1ReCtrCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/EpdReCtrEpCent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1ReCtrDist[iCent]->SaveAs(figName.c_str());
  }

  string outputFileReCtrEp = Form("../../data/%s/EventPlaneMaker/file_EpdReCtrEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileReCtrEp.c_str() << endl;
  TFile *file_OutPutReCtrEp = new TFile(outputFileReCtrEp.c_str(),"RECREATE");
  file_OutPutReCtrEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1ReCtrEast[iCent]->Write();
    h_mEpdEp1ReCtrWest[iCent]->Write();
    h_mEpdEp1ReCtrFull[iCent]->Write();
    h_mEpdEp1ReCtrCorr[iCent]->Write();
  }
  file_OutPutReCtrEp->Close();
}
