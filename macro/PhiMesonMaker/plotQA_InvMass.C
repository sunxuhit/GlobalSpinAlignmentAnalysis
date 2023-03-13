#include <string>
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void plotQA_InvMass(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumCentrality = 9;

  string inputFileSE = Form("../../data/PhiMesonMaker/%s/file_RecoPhiSE_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPutSE = TFile::Open(inputFileSE.c_str());
  if(!file_InPutSE->IsOpen()) cout << "inputFileSE: " << inputFileSE.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFileSE.c_str() << endl;

  TH2F *h_mInvMassPhiSE[mNumCentrality];
  TH1F *h_mInvMassPhiSEDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    // std::string histName = Form("h_mInvMassPhiCent%d",iCent);
    std::string histName = Form("h_mPhiMass2Cent%d",iCent);
    h_mInvMassPhiSE[iCent] = (TH2F*)file_InPutSE->Get(histName.c_str());
    histName = Form("h_mInvMassPhiSECent%d",iCent);
    h_mInvMassPhiSEDist[iCent] = (TH1F*)h_mInvMassPhiSE[iCent]->ProjectionY()->Clone(histName.c_str());
    h_mInvMassPhiSEDist[iCent]->Sumw2();
  }

  string inputFileME = Form("../../data/PhiMesonMaker/%s/file_RecoPhiME_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPutME = TFile::Open(inputFileME.c_str());
  if(!file_InPutME->IsOpen()) cout << "inputFileME: " << inputFileME.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFileME.c_str() << endl;

  TH2F *h_mInvMassPhiME[mNumCentrality]; // phi wgt
  TH1F *h_mInvMassPhiMEDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    // std::string histName = Form("h_mInvMassPhiCent%d",iCent);
    std::string histName = Form("h_mPhiMass2Cent%d",iCent);
    h_mInvMassPhiME[iCent] = (TH2F*)file_InPutME->Get(histName.c_str());
    histName = Form("h_mInvMassPhiMECent%d",iCent);
    h_mInvMassPhiMEDist[iCent] = (TH1F*)h_mInvMassPhiME[iCent]->ProjectionY()->Clone(histName.c_str());
    h_mInvMassPhiMEDist[iCent]->Sumw2();
  }

  TCanvas *c_InvMassPhi = new TCanvas("c_InvMassPhi","c_InvMassPhi",10,10,900,900);
  c_InvMassPhi->Divide(3,3);
  for(int iPad = 0; iPad < 9; ++iPad)
  {
    c_InvMassPhi->cd(iPad+1)->SetLeftMargin(0.15);
    c_InvMassPhi->cd(iPad+1)->SetRightMargin(0.15);
    c_InvMassPhi->cd(iPad+1)->SetBottomMargin(0.15);
    c_InvMassPhi->cd(iPad+1)->SetTicks(1,1);
    c_InvMassPhi->cd(iPad+1)->SetGrid(0,0);
  }

  for(int iCent = 0; iCent< 9; ++iCent)
  {
    c_InvMassPhi->cd(iCent+1);

    int binInteLo = h_mInvMassPhiSEDist[iCent]->FindBin(1.04);
    int binInteHi = h_mInvMassPhiSEDist[iCent]->FindBin(1.08);

    float yieldInteSE = h_mInvMassPhiSEDist[iCent]->Integral(binInteLo,binInteHi);
    float yieldInteME = h_mInvMassPhiMEDist[iCent]->Integral(binInteLo,binInteHi);

    h_mInvMassPhiMEDist[iCent]->Scale(yieldInteSE/yieldInteME);

    h_mInvMassPhiSEDist[iCent]->SetSt
    h_mInvMassPhiSEDist[iCent]->SetMarkerStyle(24);
    h_mInvMassPhiSEDist[iCent]->SetMarkerSize(1.5);
    h_mInvMassPhiSEDist[iCent]->SetMarkerColor(kGray+3);
    h_mInvMassPhiSEDist[iCent]->Draw("pE");
    h_mInvMassPhiMEDist[iCent]->SetFillColor(2);
    h_mInvMassPhiMEDist[iCent]->SetFillStyle(3002);
    h_mInvMassPhiMEDist[iCent]->Draw("hE same");
  }

  std::string figName = Form("../../figures/PhiMesonMaker/%s/invMassPhiQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_InvMassPhi->SaveAs(figName.c_str());
}
