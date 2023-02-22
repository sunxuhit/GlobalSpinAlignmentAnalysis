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

void getEpdReCtrPar(int beamType = 0)
{
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_%s_ReCenterPar.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // ReCtr Correction | x axis is runIndex, y axis is Centrality
  TProfile2D *p_mEpdQ1ReCtrXEast[mNumVzBin]; // 1st EP
  TProfile2D *p_mEpdQ1ReCtrYEast[mNumVzBin];
  TProfile2D *p_mEpdQ1ReCtrXWest[mNumVzBin];
  TProfile2D *p_mEpdQ1ReCtrYWest[mNumVzBin];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mEpdQ1ReCtrXEastVz%d",iVz); // 1st EP
    p_mEpdQ1ReCtrXEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mEpdQ1ReCtrYEastVz%d",iVz);
    p_mEpdQ1ReCtrYEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mEpdQ1ReCtrXWestVz%d",iVz);
    p_mEpdQ1ReCtrXWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mEpdQ1ReCtrYWestVz%d",iVz);
    p_mEpdQ1ReCtrYWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_EpdReCenterPar.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mEpdQ1ReCtrXEast[iVz]->Write();
    p_mEpdQ1ReCtrYEast[iVz]->Write();
    p_mEpdQ1ReCtrXWest[iVz]->Write();
    p_mEpdQ1ReCtrYWest[iVz]->Write();
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mEpdEp1WgtEast[mNumCentrality]; // 1st weighted EP
  TH2F *h_mEpdEp1WgtWest[mNumCentrality];
  TH2F *h_mEpdEp1WgtFull[mNumCentrality];
  TH2F *h_mEpdEp1WgtCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1WgtEastCent%d",iCent); // 1st EP
    h_mEpdEp1WgtEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1WgtWestCent%d",iCent);
    h_mEpdEp1WgtWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1WgtFullCent%d",iCent);
    h_mEpdEp1WgtFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1WgtCorrCent%d",iCent);
    h_mEpdEp1WgtCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  TCanvas *c_EpdEp1RawDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_EpdEp1RawDistCent%d",iCent);
    c_EpdEp1RawDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1600,400);
    c_EpdEp1RawDist[iCent]->Divide(4,1);
    c_EpdEp1RawDist[iCent]->cd(1); h_mEpdEp1WgtEast[iCent]->ProjectionY()->Draw();
    c_EpdEp1RawDist[iCent]->cd(2); h_mEpdEp1WgtWest[iCent]->ProjectionY()->Draw();
    c_EpdEp1RawDist[iCent]->cd(3); h_mEpdEp1WgtFull[iCent]->ProjectionY()->Draw();
    c_EpdEp1RawDist[iCent]->cd(4); h_mEpdEp1WgtCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/EpdWgtEpCent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1RawDist[iCent]->SaveAs(figName.c_str());
  }

  string outputFileWgtEp = Form("../../data/%s/EventPlaneMaker/file_%s_EpdWgtEpDist.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileWgtEp.c_str() << endl;
  TFile *file_OutPutWgtEp = new TFile(outputFileWgtEp.c_str(),"RECREATE");
  file_OutPutWgtEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1WgtEast[iCent]->Write();
    h_mEpdEp1WgtWest[iCent]->Write();
    h_mEpdEp1WgtFull[iCent]->Write();
    h_mEpdEp1WgtCorr[iCent]->Write();
  }
  file_OutPutWgtEp->Close();
}
