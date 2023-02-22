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

void getZdcReCtrPar(int beamType = 0)
{
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_%s_ReCenterPar.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // ReCenter Correction | x axis is runIndex, y axis is Centrality
  TProfile2D *p_mZdcQReCenterVertEast[mNumVzBin];
  TProfile2D *p_mZdcQReCenterHoriEast[mNumVzBin];
  TProfile2D *p_mZdcQReCenterVertWest[mNumVzBin];
  TProfile2D *p_mZdcQReCenterHoriWest[mNumVzBin];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mZdcQReCenterVertEastVz%d",iVz);
    p_mZdcQReCenterVertEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mZdcQReCenterHoriEastVz%d",iVz);
    p_mZdcQReCenterHoriEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mZdcQReCenterVertWestVz%d",iVz);
    p_mZdcQReCenterVertWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mZdcQReCenterHoriWestVz%d",iVz);
    p_mZdcQReCenterHoriWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_ZdcReCenterPar.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mZdcQReCenterVertEast[iVz]->Write();
    p_mZdcQReCenterHoriEast[iVz]->Write();
    p_mZdcQReCenterVertWest[iVz]->Write();
    p_mZdcQReCenterHoriWest[iVz]->Write();
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mZdcEpRawEast[mNumCentrality]; // raw EP
  TH2F *h_mZdcEpRawWest[mNumCentrality];
  TH2F *h_mZdcEpRawFull[mNumCentrality]; // Qwest-QEast
  TH2F *h_mZdcEpRawCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEpRawEastCent%d",iCent);
    h_mZdcEpRawEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEpRawWestCent%d",iCent);
    h_mZdcEpRawWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mZdcEpRawFullCent%d",iCent);
    h_mZdcEpRawFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEpRawCorrCent%d",iCent);
    h_mZdcEpRawCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  TCanvas *c_ZdcEp1RawDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_EpdEp1RawDistCent%d",iCent);
    c_ZdcEp1RawDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1600,400);
    c_ZdcEp1RawDist[iCent]->Divide(4,1);
    c_ZdcEp1RawDist[iCent]->cd(1); h_mZdcEpRawEast[iCent]->ProjectionY()->Draw();
    c_ZdcEp1RawDist[iCent]->cd(2); h_mZdcEpRawWest[iCent]->ProjectionY()->Draw();
    c_ZdcEp1RawDist[iCent]->cd(3); h_mZdcEpRawFull[iCent]->ProjectionY()->Draw();
    c_ZdcEp1RawDist[iCent]->cd(4); h_mZdcEpRawCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/ZdcRawEpCent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_ZdcEp1RawDist[iCent]->SaveAs(figName.c_str());
  }

  string outputFileRawEp = Form("../../data/%s/EventPlaneMaker/file_%s_ZdcRawEpDist.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileRawEp.c_str() << endl;
  TFile *file_OutPutRawEp = new TFile(outputFileRawEp.c_str(),"RECREATE");
  file_OutPutRawEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEpRawEast[iCent]->Write();
    h_mZdcEpRawWest[iCent]->Write();
    h_mZdcEpRawFull[iCent]->Write();
    h_mZdcEpRawCorr[iCent]->Write();
  }
  file_OutPutRawEp->Close();
}
