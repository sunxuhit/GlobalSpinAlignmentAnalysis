#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void calcEpdPhiWgtCorr(int beamType = 0)
{
  string inputFile = Form("../../data/%s/EventPlaneMaker/file_%s_GainCorr.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  const int mNumCentrality = 9;
  TH2F *h_mEpdPhiWgtEast[mNumCentrality];
  TH2F *h_mEpdPhiAveEast[mNumCentrality];
  TH2F *h_mEpdPhiWgtWest[mNumCentrality];
  TH2F *h_mEpdPhiAveWest[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdPhiWgtEastCent%d",iCent);
    h_mEpdPhiWgtEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    h_mEpdPhiWgtEast[iCent]->Sumw2();
    histName = Form("h_mEpdPhiAveEastCent%d",iCent);
    h_mEpdPhiAveEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    h_mEpdPhiAveEast[iCent]->Sumw2();

    histName = Form("h_mEpdPhiWgtWestCent%d",iCent);
    h_mEpdPhiWgtWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    h_mEpdPhiWgtWest[iCent]->Sumw2();
    histName = Form("h_mEpdPhiAveWestCent%d",iCent);
    h_mEpdPhiAveWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    h_mEpdPhiAveWest[iCent]->Sumw2();
  }

  TCanvas *c_EpdPhiWgtEast = new TCanvas("c_EpdPhiWgtEast","c_EpdPhiWgtEast",10,10,900,900);
  c_EpdPhiWgtEast->Divide(3,3);
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    c_EpdPhiWgtEast->cd(iCent+1);
    h_mEpdPhiWgtEast[iCent]->SetStats(0);
    h_mEpdPhiWgtEast[iCent]->Divide(h_mEpdPhiAveEast[iCent]);
    h_mEpdPhiWgtEast[iCent]->Draw("colz");
  }
  std::string figName = Form("../../figures/%s/EventPlaneMaker/EpdPhiWgtEast_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_EpdPhiWgtEast->SaveAs(figName.c_str());

  TCanvas *c_EpdPhiWgtWest = new TCanvas("c_EpdPhiWgtWest","c_EpdPhiWgtWest",10,10,900,900);
  c_EpdPhiWgtWest->Divide(3,3);
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    c_EpdPhiWgtWest->cd(iCent+1);
    h_mEpdPhiWgtWest[iCent]->SetStats(0);
    h_mEpdPhiWgtWest[iCent]->Divide(h_mEpdPhiAveWest[iCent]);
    h_mEpdPhiWgtWest[iCent]->Draw("colz");
  }
  figName = Form("../../figures/%s/EventPlaneMaker/EpdPhiWgtWest_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_EpdPhiWgtEast->SaveAs(figName.c_str());

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/GainCorrPar/file_%s_EpdPhiWgtPar.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdPhiWgtEast[iCent]->Write();
    h_mEpdPhiWgtWest[iCent]->Write();
  }
  file_OutPut->Close();
}
