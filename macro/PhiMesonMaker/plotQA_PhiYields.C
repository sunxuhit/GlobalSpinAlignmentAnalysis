#include <string>
#include "TString.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void plotQA_PhiYields(int beamType = 2, int jobId = 14)
{
  const int mNumCentrality = 9;

  string inputFileSE = Form("../../data/PhiMesonMaker/%s/file_RecoPhiSE_%s_%d.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str(),jobId);
  TFile *file_InPutSE = TFile::Open(inputFileSE.c_str());
  if(!file_InPutSE->IsOpen()) cout << "inputFileSE: " << inputFileSE.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFileSE.c_str() << endl;

  TH2F *h_mInvMassPhiSE[mNumCentrality];
  TH1F *h_mInvMassPhiSEDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    // std::string histName = Form("h_mInvMassPhiCent%d",iCent);
    std::string histName = Form("h_mInvMassPhiCent%d",iCent);
    h_mInvMassPhiSE[iCent] = (TH2F*)file_InPutSE->Get(histName.c_str());
    histName = Form("h_mInvMassPhiSECent%d",iCent);
    h_mInvMassPhiSEDist[iCent] = (TH1F*)h_mInvMassPhiSE[iCent]->ProjectionY()->Clone(histName.c_str());
    h_mInvMassPhiSEDist[iCent]->Sumw2();
  }

  TH2F *h_mInvMassPhiMinBiasSE;
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    if(iCent == 0) h_mInvMassPhiMinBiasSE = (TH2F*)h_mInvMassPhiSE[iCent]->Clone("h_mInvMassPhiCent9");
    else h_mInvMassPhiMinBiasSE->Add(h_mInvMassPhiSE[iCent],1.0);
  }
  TH1F *h_mInvMassPhiMinBiasSEDist = (TH1F*)h_mInvMassPhiMinBiasSE->ProjectionY()->Clone("h_mInvMassPhiSECent9");
  h_mInvMassPhiMinBiasSEDist->Draw("h");

  TCanvas *c_InvMassPhi = new TCanvas("c_InvMassPhi","c_InvMassPhi",10,10,2000,800);
  c_InvMassPhi->Divide(5,2);
  for(int iPad = 0; iPad < 10; ++iPad)
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
    std::string histName = Form("h_mInvMassPhiCent%d",iCent);
    h_mInvMassPhiSEDist[iCent]->SetTitle(histName.c_str());
    h_mInvMassPhiSEDist[iCent]->GetXaxis()->SetNdivisions(505);
    h_mInvMassPhiSEDist[iCent]->Draw("h");
  }
  c_InvMassPhi->cd(10);
  h_mInvMassPhiMinBiasSEDist->SetTitle("0-80%");
  h_mInvMassPhiMinBiasSEDist->GetXaxis()->SetNdivisions(505);
  h_mInvMassPhiMinBiasSEDist->Draw("h");


  std::string figName = Form("../../figures/PhiMesonMaker/%s/invMassPhiYieldsQA_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_InvMassPhi->SaveAs(figName.c_str());
}
