#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

#define _RawEp_ 1
#define _WgtEp_ 0
#define _ReCtrEp_ 0
#define _ShiftEp_ 0
#define _ShiftFullEp_ 0

void getTpcReCtrPar(int beamType = 0)
{
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_%s_ReCenterPar.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // ReCenter Correction | x axis is runIndex, y axis is Centrality
  TProfile2D *p_mTpcQ2ReCenterXEast[mNumVzBin]; // 2nd EP
  TProfile2D *p_mTpcQ2ReCenterYEast[mNumVzBin];
  TProfile2D *p_mTpcQ2ReCenterXWest[mNumVzBin];
  TProfile2D *p_mTpcQ2ReCenterYWest[mNumVzBin];

  TProfile2D *p_mTpcQ3ReCenterXEast[mNumVzBin]; // 3rd EP
  TProfile2D *p_mTpcQ3ReCenterYEast[mNumVzBin];
  TProfile2D *p_mTpcQ3ReCenterXWest[mNumVzBin];
  TProfile2D *p_mTpcQ3ReCenterYWest[mNumVzBin];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQ2ReCenterXEastVz%d",iVz); // 2nd EP
    p_mTpcQ2ReCenterXEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCenterYEastVz%d",iVz);
    p_mTpcQ2ReCenterYEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mTpcQ2ReCenterXWestVz%d",iVz);
    p_mTpcQ2ReCenterXWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCenterYWestVz%d",iVz);
    p_mTpcQ2ReCenterYWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mTpcQ3ReCenterXEastVz%d",iVz); // 3rd EP
    p_mTpcQ3ReCenterXEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcQ3ReCenterYEastVz%d",iVz);
    p_mTpcQ3ReCenterYEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mTpcQ3ReCenterXWestVz%d",iVz);
    p_mTpcQ3ReCenterXWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcQ3ReCenterYWestVz%d",iVz);
    p_mTpcQ3ReCenterYWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_TpcReCenterPar.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mTpcQ2ReCenterXEast[iVz]->Write();
    p_mTpcQ2ReCenterYEast[iVz]->Write();
    p_mTpcQ2ReCenterXWest[iVz]->Write();
    p_mTpcQ2ReCenterYWest[iVz]->Write();

    p_mTpcQ3ReCenterXEast[iVz]->Write();
    p_mTpcQ3ReCenterYEast[iVz]->Write();
    p_mTpcQ3ReCenterXWest[iVz]->Write();
    p_mTpcQ3ReCenterYWest[iVz]->Write();
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mTpcEp2RawEast[mNumCentrality]; // 2nd raw EP
  TH2F *h_mTpcEp2RawWest[mNumCentrality];
  TH2F *h_mTpcEp2RawCorr[mNumCentrality]; // Psi2East vs Psi2West

  TH2F *h_mTpcEp3RawEast[mNumCentrality]; // 3rd raw EP
  TH2F *h_mTpcEp3RawWest[mNumCentrality];
  TH2F *h_mTpcEp3RawCorr[mNumCentrality]; // Psi3East vs Psi3West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp2RawEastCent%d",iCent); // 2nd EP
    h_mTpcEp2RawEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2RawWestCent%d",iCent);
    h_mTpcEp2RawWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2RawCorrCent%d",iCent);
    h_mTpcEp2RawCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mTpcEp3RawEastCent%d",iCent); // 3rd EP
    h_mTpcEp3RawEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3RawWestCent%d",iCent);
    h_mTpcEp3RawWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3RawCorrCent%d",iCent);
    h_mTpcEp3RawCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  TCanvas *c_TpcEp2RawDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_TpcEp2RawDistCent%d",iCent);
    c_TpcEp2RawDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1500,500);
    c_TpcEp2RawDist[iCent]->Divide(3,1);
    c_TpcEp2RawDist[iCent]->cd(1); h_mTpcEp2RawEast[iCent]->ProjectionY()->Draw();
    c_TpcEp2RawDist[iCent]->cd(2); h_mTpcEp2RawWest[iCent]->ProjectionY()->Draw();
    c_TpcEp2RawDist[iCent]->cd(3); h_mTpcEp2RawCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/TpcRawEp2Cent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_TpcEp2RawDist[iCent]->SaveAs(figName.c_str());
  }

  TCanvas *c_TpcEp3RawDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_TpcEp3RawDistCent%d",iCent);
    c_TpcEp3RawDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1500,500);
    c_TpcEp3RawDist[iCent]->Divide(3,1);
    c_TpcEp3RawDist[iCent]->cd(1); h_mTpcEp3RawEast[iCent]->ProjectionY()->Draw();
    c_TpcEp3RawDist[iCent]->cd(2); h_mTpcEp3RawWest[iCent]->ProjectionY()->Draw();
    c_TpcEp3RawDist[iCent]->cd(3); h_mTpcEp3RawCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/TpcRawEp3Cent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_TpcEp3RawDist[iCent]->SaveAs(figName.c_str());
  }

  string outputFileWgtEp = Form("../../data/%s/EventPlaneMaker/file_%s_TpcRawEpDist.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileWgtEp.c_str() << endl;
  TFile *file_OutPutWgtEp = new TFile(outputFileWgtEp.c_str(),"RECREATE");
  file_OutPutWgtEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2RawEast[iCent]->Write();
    h_mTpcEp2RawWest[iCent]->Write();
    h_mTpcEp2RawCorr[iCent]->Write();

    h_mTpcEp3RawEast[iCent]->Write();
    h_mTpcEp3RawWest[iCent]->Write();
    h_mTpcEp3RawCorr[iCent]->Write();
  }
  file_OutPutWgtEp->Close();
}
