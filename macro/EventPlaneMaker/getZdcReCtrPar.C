#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getZdcReCtrPar(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_ReCenterPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // ReCenter Correction | x axis is runIndex, y axis is Centrality
  TProfile2D *p_mZdcQ1ReCtrVertEast[mNumVzBin];
  TProfile2D *p_mZdcQ1ReCtrHoriEast[mNumVzBin];
  TProfile2D *p_mZdcQ1ReCtrVertWest[mNumVzBin];
  TProfile2D *p_mZdcQ1ReCtrHoriWest[mNumVzBin];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mZdcQ1ReCtrVertEastVz%d",iVz);
    p_mZdcQ1ReCtrVertEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mZdcQ1ReCtrHoriEastVz%d",iVz);
    p_mZdcQ1ReCtrHoriEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mZdcQ1ReCtrVertWestVz%d",iVz);
    p_mZdcQ1ReCtrVertWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mZdcQ1ReCtrHoriWestVz%d",iVz);
    p_mZdcQ1ReCtrHoriWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
  }

  {
    TCanvas *c_ZdcQ1ReCtr = new TCanvas("c_ZdcQ1ReCtr","c_ZdcQ1ReCtr",10,10,800,800);
    c_ZdcQ1ReCtr->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_ZdcQ1ReCtr->cd(iPad+1)->SetLeftMargin(0.15);
      c_ZdcQ1ReCtr->cd(iPad+1)->SetRightMargin(0.15);
      c_ZdcQ1ReCtr->cd(iPad+1)->SetBottomMargin(0.15);
      c_ZdcQ1ReCtr->cd(iPad+1)->SetTicks(1,1);
      c_ZdcQ1ReCtr->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/ZdcQ1ReCtr_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcQ1ReCtr->Print(figName.c_str());
    for(int iVz = 0; iVz < mNumVzBin; ++iVz)
    {
      figName = Form("../../figures/EventPlaneMaker/%s/ZdcQ1ReCtr_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_ZdcQ1ReCtr->cd(1)->Clear(); c_ZdcQ1ReCtr->cd(1); p_mZdcQ1ReCtrVertEast[iVz]->Draw("colz");
      c_ZdcQ1ReCtr->cd(2)->Clear(); c_ZdcQ1ReCtr->cd(2); p_mZdcQ1ReCtrHoriEast[iVz]->Draw("colz");
      c_ZdcQ1ReCtr->cd(3)->Clear(); c_ZdcQ1ReCtr->cd(3); p_mZdcQ1ReCtrVertWest[iVz]->Draw("colz");
      c_ZdcQ1ReCtr->cd(4)->Clear(); c_ZdcQ1ReCtr->cd(4); p_mZdcQ1ReCtrHoriWest[iVz]->Draw("colz");
      c_ZdcQ1ReCtr->Update();
      c_ZdcQ1ReCtr->Print(figName.c_str());
    }
    figName = Form("../../figures/EventPlaneMaker/%s/ZdcQ1ReCtr_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcQ1ReCtr->Print(figName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ReCenterPar/file_ZdcReCenterPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mZdcQ1ReCtrVertEast[iVz]->Write();
    p_mZdcQ1ReCtrHoriEast[iVz]->Write();
    p_mZdcQ1ReCtrVertWest[iVz]->Write();
    p_mZdcQ1ReCtrHoriWest[iVz]->Write();
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mZdcEp1RawEast[mNumCentrality]; // raw EP
  TH2F *h_mZdcEp1RawWest[mNumCentrality];
  TH2F *h_mZdcEp1RawFull[mNumCentrality]; // Qwest-QEast
  TH2F *h_mZdcEp1RawCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEp1RawEastCent%d",iCent);
    h_mZdcEp1RawEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEp1RawWestCent%d",iCent);
    h_mZdcEp1RawWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mZdcEp1RawFullCent%d",iCent);
    h_mZdcEp1RawFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEp1RawCorrCent%d",iCent);
    h_mZdcEp1RawCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  {
    TCanvas *c_ZdcEp1RawDist = new TCanvas("c_ZdcEp1RawDist","c_ZdcEp1RawDist",10,10,800,800);
    c_ZdcEp1RawDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_ZdcEp1RawDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_ZdcEp1RawDist->cd(iPad+1)->SetRightMargin(0.15);
      c_ZdcEp1RawDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_ZdcEp1RawDist->cd(iPad+1)->SetTicks(1,1);
      c_ZdcEp1RawDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/ZdcRawEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcEp1RawDist->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      figName = Form("../../figures/EventPlaneMaker/%s/ZdcRawEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_ZdcEp1RawDist->cd(1)->Clear(); c_ZdcEp1RawDist->cd(1); h_mZdcEp1RawEast[iCent]->ProjectionY()->Draw();
      c_ZdcEp1RawDist->cd(2)->Clear(); c_ZdcEp1RawDist->cd(2); h_mZdcEp1RawWest[iCent]->ProjectionY()->Draw();
      c_ZdcEp1RawDist->cd(3)->Clear(); c_ZdcEp1RawDist->cd(3); h_mZdcEp1RawFull[iCent]->ProjectionY()->Draw();
      c_ZdcEp1RawDist->cd(4)->Clear(); c_ZdcEp1RawDist->cd(4); h_mZdcEp1RawCorr[iCent]->Draw("colz");
      c_ZdcEp1RawDist->Update();
      c_ZdcEp1RawDist->Print(figName.c_str());
    }
    figName = Form("../../figures/EventPlaneMaker/%s/ZdcRawEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_ZdcEp1RawDist->Print(figName.c_str());
  }

  string outputFileRawEp = Form("../../data/EventPlaneMaker/%s/file_ZdcRawEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileRawEp.c_str() << endl;
  TFile *file_OutPutRawEp = new TFile(outputFileRawEp.c_str(),"RECREATE");
  file_OutPutRawEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEp1RawEast[iCent]->Write();
    h_mZdcEp1RawWest[iCent]->Write();
    h_mZdcEp1RawFull[iCent]->Write();
    h_mZdcEp1RawCorr[iCent]->Write();
  }
  file_OutPutRawEp->Close();
}
