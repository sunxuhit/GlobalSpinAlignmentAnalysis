#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getTpcShiftPar(int beamType = 0)
{
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumShiftCorr = 20;
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_ShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
  TProfile2D *p_mTpcQ1ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 1st EP
  TProfile2D *p_mTpcQ1ShiftSinEast[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ1ShiftCosWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ1ShiftSinWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ1ShiftCosFull[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ1ShiftSinFull[mNumVzBin][mNumShiftCorr];

  TProfile2D *p_mTpcQ2ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 2nd EP
  TProfile2D *p_mTpcQ2ShiftSinEast[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ2ShiftCosWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ2ShiftSinWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ2ShiftCosFull[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ2ShiftSinFull[mNumVzBin][mNumShiftCorr];

  TProfile2D *p_mTpcQ3ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 3rd EP
  TProfile2D *p_mTpcQ3ShiftSinEast[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ3ShiftCosWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ3ShiftSinWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ3ShiftCosFull[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ3ShiftSinFull[mNumVzBin][mNumShiftCorr];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mTpcQ1ShiftCos%dEastVz%d",iShift,iVz); // 1st EP
      p_mTpcQ1ShiftCosEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ1ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ1ShiftSinEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ1ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ1ShiftCosWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ1ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ1ShiftSinWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ1ShiftCos%dFullVz%d",iShift,iVz);
      p_mTpcQ1ShiftCosFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ1ShiftSin%dFullVz%d",iShift,iVz);
      p_mTpcQ1ShiftSinFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ2ShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mTpcQ2ShiftCosEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ2ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ2ShiftCos%dFullVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dFullVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dEastVz%d",iShift,iVz); // 3rd EP
      p_mTpcQ3ShiftCosEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftCosWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dFullVz%d",iShift,iVz);
      p_mTpcQ3ShiftCosFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dFullVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinFull[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
    }
  }

  {
    TCanvas *c_TpcQVecShift = new TCanvas("c_TpcQVecShift","c_TpcQVecShift",10,10,800,1200);
    c_TpcQVecShift->Divide(2,3);
    for(int iPad = 0; iPad < 6; ++iPad)
    {
      c_TpcQVecShift->cd(iPad+1)->SetLeftMargin(0.15);
      c_TpcQVecShift->cd(iPad+1)->SetRightMargin(0.15);
      c_TpcQVecShift->cd(iPad+1)->SetBottomMargin(0.15);
      c_TpcQVecShift->cd(iPad+1)->SetTicks(1,1);
      c_TpcQVecShift->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/TpcQVecShift_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_TpcQVecShift->Print(figName.c_str());
    for(int iVz = 0; iVz < mNumVzBin; ++iVz)
    {
      for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
      {
	figName = Form("../../figures/EventPlaneMaker/%s/TpcQVecShift_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
	c_TpcQVecShift->cd(1)->Clear(); c_TpcQVecShift->cd(1); p_mTpcQ1ShiftCosEast[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(2)->Clear(); c_TpcQVecShift->cd(2); p_mTpcQ1ShiftSinEast[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(3)->Clear(); c_TpcQVecShift->cd(3); p_mTpcQ1ShiftCosWest[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(4)->Clear(); c_TpcQVecShift->cd(4); p_mTpcQ1ShiftSinWest[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(5)->Clear(); c_TpcQVecShift->cd(5); p_mTpcQ1ShiftCosFull[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(6)->Clear(); c_TpcQVecShift->cd(6); p_mTpcQ1ShiftSinFull[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->Update();
	c_TpcQVecShift->Print(figName.c_str());

	c_TpcQVecShift->cd(1)->Clear(); c_TpcQVecShift->cd(1); p_mTpcQ2ShiftCosEast[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(2)->Clear(); c_TpcQVecShift->cd(2); p_mTpcQ2ShiftSinEast[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(3)->Clear(); c_TpcQVecShift->cd(3); p_mTpcQ2ShiftCosWest[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(4)->Clear(); c_TpcQVecShift->cd(4); p_mTpcQ2ShiftSinWest[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(5)->Clear(); c_TpcQVecShift->cd(5); p_mTpcQ2ShiftCosFull[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(6)->Clear(); c_TpcQVecShift->cd(6); p_mTpcQ2ShiftSinFull[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->Update();
	c_TpcQVecShift->Print(figName.c_str());

	c_TpcQVecShift->cd(1)->Clear(); c_TpcQVecShift->cd(1); p_mTpcQ3ShiftCosEast[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(2)->Clear(); c_TpcQVecShift->cd(2); p_mTpcQ3ShiftSinEast[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(3)->Clear(); c_TpcQVecShift->cd(3); p_mTpcQ3ShiftCosWest[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(4)->Clear(); c_TpcQVecShift->cd(4); p_mTpcQ3ShiftSinWest[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(5)->Clear(); c_TpcQVecShift->cd(5); p_mTpcQ3ShiftCosFull[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->cd(6)->Clear(); c_TpcQVecShift->cd(6); p_mTpcQ3ShiftSinFull[iVz][iShift]->DrawCopy("colz");
	c_TpcQVecShift->Update();
	c_TpcQVecShift->Print(figName.c_str());
      }
    }
    figName = Form("../../figures/EventPlaneMaker/%s/TpcQVecShift_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_TpcQVecShift->Print(figName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ShiftPar/file_TpcShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mTpcQ1ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ1ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ1ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ1ShiftSinWest[iVz][iShift]->Write();
      p_mTpcQ1ShiftCosFull[iVz][iShift]->Write();
      p_mTpcQ1ShiftSinFull[iVz][iShift]->Write();

      p_mTpcQ2ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ2ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinWest[iVz][iShift]->Write();
      p_mTpcQ2ShiftCosFull[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinFull[iVz][iShift]->Write();

      p_mTpcQ3ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ3ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinWest[iVz][iShift]->Write();
      p_mTpcQ3ShiftCosFull[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinFull[iVz][iShift]->Write();
    }
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mTpcEp1ReCtrEast[mNumCentrality]; // 2nd raw EP
  TH2F *h_mTpcEp1ReCtrWest[mNumCentrality];
  TH2F *h_mTpcEp1ReCtrFull[mNumCentrality];
  TH2F *h_mTpcEp1ReCtrCorr[mNumCentrality]; // Psi2East vs Psi2West

  TH2F *h_mTpcEp2ReCtrEast[mNumCentrality]; // 2nd raw EP
  TH2F *h_mTpcEp2ReCtrWest[mNumCentrality];
  TH2F *h_mTpcEp2ReCtrFull[mNumCentrality];
  TH2F *h_mTpcEp2ReCtrCorr[mNumCentrality]; // Psi2East vs Psi2West

  TH2F *h_mTpcEp3ReCtrEast[mNumCentrality]; // 3rd raw EP
  TH2F *h_mTpcEp3ReCtrWest[mNumCentrality];
  TH2F *h_mTpcEp3ReCtrFull[mNumCentrality];
  TH2F *h_mTpcEp3ReCtrCorr[mNumCentrality]; // Psi3East vs Psi3West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp1ReCtrEastCent%d",iCent); // 1st EP
    h_mTpcEp1ReCtrEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp1ReCtrWestCent%d",iCent);
    h_mTpcEp1ReCtrWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp1ReCtrFullCent%d",iCent);
    h_mTpcEp1ReCtrFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp1ReCtrCorrCent%d",iCent);
    h_mTpcEp1ReCtrCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mTpcEp2ReCtrEastCent%d",iCent); // 2nd EP
    h_mTpcEp2ReCtrEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2ReCtrWestCent%d",iCent);
    h_mTpcEp2ReCtrWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2ReCtrFullCent%d",iCent);
    h_mTpcEp2ReCtrFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2ReCtrCorrCent%d",iCent);
    h_mTpcEp2ReCtrCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mTpcEp3ReCtrEastCent%d",iCent); // 3rd EP
    h_mTpcEp3ReCtrEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3ReCtrWestCent%d",iCent);
    h_mTpcEp3ReCtrWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3ReCtrFullCent%d",iCent);
    h_mTpcEp3ReCtrFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3ReCtrCorrCent%d",iCent);
    h_mTpcEp3ReCtrCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  {
    TCanvas *c_TpcEpReCtrDist = new TCanvas("c_TpcEpReCtrDist","c_TpcEpReCtrDist",10,10,800,800);
    c_TpcEpReCtrDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_TpcEpReCtrDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_TpcEpReCtrDist->cd(iPad+1)->SetRightMargin(0.15);
      c_TpcEpReCtrDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_TpcEpReCtrDist->cd(iPad+1)->SetTicks(1,1);
      c_TpcEpReCtrDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/TpcReCtrEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_TpcEpReCtrDist->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      figName = Form("../../figures/EventPlaneMaker/%s/TpcReCtrEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_TpcEpReCtrDist->cd(1); h_mTpcEp1ReCtrEast[iCent]->ProjectionY()->DrawCopy();
      c_TpcEpReCtrDist->cd(2); h_mTpcEp1ReCtrWest[iCent]->ProjectionY()->DrawCopy();
      c_TpcEpReCtrDist->cd(3); h_mTpcEp1ReCtrFull[iCent]->ProjectionY()->DrawCopy();
      c_TpcEpReCtrDist->cd(4); h_mTpcEp1ReCtrCorr[iCent]->DrawCopy("colz");
      c_TpcEpReCtrDist->Update();
      c_TpcEpReCtrDist->Print(figName.c_str());

      c_TpcEpReCtrDist->cd(1); h_mTpcEp2ReCtrEast[iCent]->ProjectionY()->DrawCopy();
      c_TpcEpReCtrDist->cd(2); h_mTpcEp2ReCtrWest[iCent]->ProjectionY()->DrawCopy();
      c_TpcEpReCtrDist->cd(3); h_mTpcEp2ReCtrFull[iCent]->ProjectionY()->DrawCopy();
      c_TpcEpReCtrDist->cd(4); h_mTpcEp2ReCtrCorr[iCent]->DrawCopy("colz");
      c_TpcEpReCtrDist->Update();
      c_TpcEpReCtrDist->Print(figName.c_str());

      c_TpcEpReCtrDist->cd(1); h_mTpcEp3ReCtrEast[iCent]->ProjectionY()->DrawCopy();
      c_TpcEpReCtrDist->cd(2); h_mTpcEp3ReCtrWest[iCent]->ProjectionY()->DrawCopy();
      c_TpcEpReCtrDist->cd(3); h_mTpcEp3ReCtrFull[iCent]->ProjectionY()->DrawCopy();
      c_TpcEpReCtrDist->cd(4); h_mTpcEp3ReCtrCorr[iCent]->DrawCopy("colz");
      c_TpcEpReCtrDist->Update();
      c_TpcEpReCtrDist->Print(figName.c_str());
    }
    figName = Form("../../figures/EventPlaneMaker/%s/TpcReCtrEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_TpcEpReCtrDist->Print(figName.c_str());
  }

  string outputFileReCtrEp = Form("../../data/EventPlaneMaker/%s/file_TpcReCtrEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of ReCtr EP: " << outputFileReCtrEp.c_str() << endl;
  TFile *file_OutPutReCtrEp = new TFile(outputFileReCtrEp.c_str(),"RECREATE");
  file_OutPutReCtrEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp1ReCtrEast[iCent]->Write();
    h_mTpcEp1ReCtrWest[iCent]->Write();
    h_mTpcEp1ReCtrFull[iCent]->Write();
    h_mTpcEp1ReCtrCorr[iCent]->Write();

    h_mTpcEp2ReCtrEast[iCent]->Write();
    h_mTpcEp2ReCtrWest[iCent]->Write();
    h_mTpcEp2ReCtrFull[iCent]->Write();
    h_mTpcEp2ReCtrCorr[iCent]->Write();

    h_mTpcEp3ReCtrEast[iCent]->Write();
    h_mTpcEp3ReCtrWest[iCent]->Write();
    h_mTpcEp3ReCtrFull[iCent]->Write();
    h_mTpcEp3ReCtrCorr[iCent]->Write();
  }
  file_OutPutReCtrEp->Close();
  file_InPut->Close();
}
