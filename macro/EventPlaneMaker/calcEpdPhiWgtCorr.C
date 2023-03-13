#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void calcEpdPhiWgtCorr(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_GainCorr_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

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

  {
    TCanvas *c_EpdPhiWgt = new TCanvas("c_EpdPhiWgt","c_EpdPhiWgt",10,10,900,900);
    c_EpdPhiWgt->Divide(3,3);
    for(int iPad = 0; iPad < 9; ++iPad)
    {
      c_EpdPhiWgt->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdPhiWgt->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdPhiWgt->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdPhiWgt->cd(iPad+1)->SetTicks(1,1);
      c_EpdPhiWgt->cd(iPad+1)->SetGrid(0,0);
    }
    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdPhiWgt_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdPhiWgt->Print(figName.c_str());

    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      c_EpdPhiWgt->cd(iCent+1);
      h_mEpdPhiWgtEast[iCent]->SetStats(0);
      h_mEpdPhiWgtEast[iCent]->Divide(h_mEpdPhiAveEast[iCent]);
      h_mEpdPhiWgtEast[iCent]->Draw("colz");
    }
    figName = Form("../../figures/EventPlaneMaker/%s/EpdPhiWgt_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdPhiWgt->Update();
    c_EpdPhiWgt->Print(figName.c_str());

    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      c_EpdPhiWgt->cd(iCent+1);
      h_mEpdPhiWgtWest[iCent]->SetStats(0);
      h_mEpdPhiWgtWest[iCent]->Divide(h_mEpdPhiAveWest[iCent]);
      h_mEpdPhiWgtWest[iCent]->Draw("colz");
    }
    c_EpdPhiWgt->Update();
    c_EpdPhiWgt->Print(figName.c_str());

    figName = Form("../../figures/EventPlaneMaker/%s/EpdPhiWgt_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdPhiWgt->Print(figName.c_str());
  }

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

  {
    TCanvas *c_EpdEp1RawDist = new TCanvas("c_EpdEp1RawDist","c_EpdEp1RawDist",10,10,800,800);
    c_EpdEp1RawDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_EpdEp1RawDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdEp1RawDist->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdEp1RawDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdEp1RawDist->cd(iPad+1)->SetTicks(1,1);
      c_EpdEp1RawDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdRawEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1RawDist->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      figName = Form("../../figures/EventPlaneMaker/%s/EpdRawEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_EpdEp1RawDist->cd(1)->Clear(); c_EpdEp1RawDist->cd(1); h_mEpdEp1RawEast[iCent]->ProjectionY()->Draw();
      c_EpdEp1RawDist->cd(2)->Clear(); c_EpdEp1RawDist->cd(2); h_mEpdEp1RawWest[iCent]->ProjectionY()->Draw();
      c_EpdEp1RawDist->cd(3)->Clear(); c_EpdEp1RawDist->cd(3); h_mEpdEp1RawFull[iCent]->ProjectionY()->Draw();
      c_EpdEp1RawDist->cd(4)->Clear(); c_EpdEp1RawDist->cd(4); h_mEpdEp1RawCorr[iCent]->Draw("colz");
      c_EpdEp1RawDist->Update();
      c_EpdEp1RawDist->Print(figName.c_str());
    }
    figName = Form("../../figures/EventPlaneMaker/%s/EpdRawEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1RawDist->Print(figName.c_str());
  }

  string outputFileRawEp = Form("../../data/EventPlaneMaker/%s/file_EpdRawEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
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
