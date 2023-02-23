#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void calcEpdPhiWgtCorr(int beamType = 0)
{
  string inputFile = Form("../../data/%s/EventPlaneMaker/file_GainCorr_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  const int mNumCentrality = 9;
  TH2F *h_mEpdPhiWgtEast[mNumCentrality]; // phi wgt
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

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/GainCorrPar/file_EpdPhiWgtPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdPhiWgtEast[iCent]->Write();
    h_mEpdPhiWgtWest[iCent]->Write();
  }
  file_OutPut->Close();

  //----------------------------------------------------------------

  TH2F *h_mEpdEp1RawEast[mNumCentrality]; // 1st raw EP
  TH2F *h_mEpdEp1RawWest[mNumCentrality];
  TH2F *h_mEpdEp1RawFull[mNumCentrality];
  TH2F *h_mEpdEp1RawCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1RawEastCent%d",iCent); // 2nd EP
    h_mEpdEp1RawEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1RawWestCent%d",iCent);
    h_mEpdEp1RawWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1RawFullCent%d",iCent);
    h_mEpdEp1RawFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1RawCorrCent%d",iCent);
    h_mEpdEp1RawCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  TCanvas *c_EpdEp1RawDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_EpdEp1RawDistCent%d",iCent);
    c_EpdEp1RawDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1600,400);
    c_EpdEp1RawDist[iCent]->Divide(4,1);
    c_EpdEp1RawDist[iCent]->cd(1); h_mEpdEp1RawEast[iCent]->ProjectionY()->Draw();
    c_EpdEp1RawDist[iCent]->cd(2); h_mEpdEp1RawWest[iCent]->ProjectionY()->Draw();
    c_EpdEp1RawDist[iCent]->cd(3); h_mEpdEp1RawFull[iCent]->ProjectionY()->Draw();
    c_EpdEp1RawDist[iCent]->cd(4); h_mEpdEp1RawCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/EpdRawEpCent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1RawDist[iCent]->SaveAs(figName.c_str());
  }

  string outputFileRawEp = Form("../../data/%s/EventPlaneMaker/file_EpdRawEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileRawEp.c_str() << endl;
  TFile *file_OutPutRawEp = new TFile(outputFileRawEp.c_str(),"RECREATE");
  file_OutPutRawEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1RawEast[iCent]->Write();
    h_mEpdEp1RawWest[iCent]->Write();
    h_mEpdEp1RawFull[iCent]->Write();
    h_mEpdEp1RawCorr[iCent]->Write();
  }
  file_OutPutRawEp->Close();
}
