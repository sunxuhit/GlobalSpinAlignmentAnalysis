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

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_ShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
  TProfile2D *p_mTpcQ2ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 2nd EP
  TProfile2D *p_mTpcQ2ShiftSinEast[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ2ShiftCosWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ2ShiftSinWest[mNumVzBin][mNumShiftCorr];

  TProfile2D *p_mTpcQ3ShiftCosEast[mNumVzBin][mNumShiftCorr]; // 3rd EP
  TProfile2D *p_mTpcQ3ShiftSinEast[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ3ShiftCosWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mTpcQ3ShiftSinWest[mNumVzBin][mNumShiftCorr];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mTpcQ2ShiftCos%dEastVz%d",iShift,iVz); // 2nd EP
      p_mTpcQ2ShiftCosEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ2ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftCosWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ2ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ2ShiftSinWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dEastVz%d",iShift,iVz); // 3rd EP
      p_mTpcQ3ShiftCosEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dEastVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mTpcQ3ShiftCos%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftCosWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mTpcQ3ShiftSin%dWestVz%d",iShift,iVz);
      p_mTpcQ3ShiftSinWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
    }
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ShiftPar/file_TpcShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mTpcQ2ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ2ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ2ShiftSinWest[iVz][iShift]->Write();

      p_mTpcQ3ShiftCosEast[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinEast[iVz][iShift]->Write();
      p_mTpcQ3ShiftCosWest[iVz][iShift]->Write();
      p_mTpcQ3ShiftSinWest[iVz][iShift]->Write();
    }
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mTpcEp2ReCtrEast[mNumCentrality]; // 2nd raw EP
  TH2F *h_mTpcEp2ReCtrWest[mNumCentrality];
  TH2F *h_mTpcEp2ReCtrCorr[mNumCentrality]; // Psi2East vs Psi2West

  TH2F *h_mTpcEp3ReCtrEast[mNumCentrality]; // 3rd raw EP
  TH2F *h_mTpcEp3ReCtrWest[mNumCentrality];
  TH2F *h_mTpcEp3ReCtrCorr[mNumCentrality]; // Psi3East vs Psi3West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp2ReCtrEastCent%d",iCent); // 2nd EP
    h_mTpcEp2ReCtrEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2ReCtrWestCent%d",iCent);
    h_mTpcEp2ReCtrWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2ReCtrCorrCent%d",iCent);
    h_mTpcEp2ReCtrCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mTpcEp3ReCtrEastCent%d",iCent); // 3rd EP
    h_mTpcEp3ReCtrEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3ReCtrWestCent%d",iCent);
    h_mTpcEp3ReCtrWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3ReCtrCorrCent%d",iCent);
    h_mTpcEp3ReCtrCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  TCanvas *c_TpcEpReCtrDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_TpcEpReCtrDistCent%d",iCent);
    c_TpcEpReCtrDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1500,1000);
    c_TpcEpReCtrDist[iCent]->Divide(3,2);
    c_TpcEpReCtrDist[iCent]->cd(1); h_mTpcEp2ReCtrEast[iCent]->ProjectionY()->Draw();
    c_TpcEpReCtrDist[iCent]->cd(2); h_mTpcEp2ReCtrWest[iCent]->ProjectionY()->Draw();
    c_TpcEpReCtrDist[iCent]->cd(3); h_mTpcEp2ReCtrCorr[iCent]->Draw("colz");

    c_TpcEpReCtrDist[iCent]->cd(4); h_mTpcEp3ReCtrEast[iCent]->ProjectionY()->Draw();
    c_TpcEpReCtrDist[iCent]->cd(5); h_mTpcEp3ReCtrWest[iCent]->ProjectionY()->Draw();
    c_TpcEpReCtrDist[iCent]->cd(6); h_mTpcEp3ReCtrCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/TpcReCtrEpCent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_TpcEpReCtrDist[iCent]->SaveAs(figName.c_str());
  }

  string outputFileReCtrEp = Form("../../data/%s/EventPlaneMaker/file_TpcReCtrEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of ReCtr EP: " << outputFileReCtrEp.c_str() << endl;
  TFile *file_OutPutReCtrEp = new TFile(outputFileReCtrEp.c_str(),"RECREATE");
  file_OutPutReCtrEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2ReCtrEast[iCent]->Write();
    h_mTpcEp2ReCtrWest[iCent]->Write();
    h_mTpcEp2ReCtrCorr[iCent]->Write();

    h_mTpcEp3ReCtrEast[iCent]->Write();
    h_mTpcEp3ReCtrWest[iCent]->Write();
    h_mTpcEp3ReCtrCorr[iCent]->Write();
  }
  file_OutPutReCtrEp->Close();
}
