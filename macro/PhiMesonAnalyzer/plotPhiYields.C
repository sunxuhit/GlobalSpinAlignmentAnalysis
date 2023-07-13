#include <string>
#include <map>

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void plotPhiYields(int beamType = 2)
{
  const int mNumCentBinFxtPhiYileds = 9;  // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
  const int mNumPtBinFxtPhiFlow     = 10; // 0.4 - 3.0 GeV/c
  const double ptLo[mNumPtBinFxtPhiFlow] = {0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
  const double ptHi[mNumPtBinFxtPhiFlow] = {0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,3.0};

  // FXT phi Yields histograms
  // 0 = centrality: 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%
  // 1 = pt bin
  std::map<string,TH1F*> h_mInvMassFxtPhiYieldsSE;
  std::map<string,TH1F*> h_mInvMassFxtPhiYieldsME;

  string inputFileSE = Form("../../data/PhiMesonAnalyzer/%s/file_FlowPhiSE_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPutSE = TFile::Open(inputFileSE.c_str());
  if(!file_InPutSE->IsOpen()) cout << "inputFileSE: " << inputFileSE.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFileSE.c_str() << endl;
  for(int iCent = 0; iCent < mNumCentBinFxtPhiYileds; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinFxtPhiFlow; ++iPt)
    {
      std::string fxtPhiYieldsKey = Form("h_mInvMassFxtPhiYieldsSECent%dPt%d",iCent,iPt);
      h_mInvMassFxtPhiYieldsSE[fxtPhiYieldsKey] = (TH1F*)file_InPutSE->Get(fxtPhiYieldsKey.c_str())->Clone();
      h_mInvMassFxtPhiYieldsSE[fxtPhiYieldsKey]->Sumw2();
      if(iCent == 2)
      {
	std::string fxtPhiYieldsMinBiasKey = Form("h_mInvMassFxtPhiYieldsSECent9Pt%d",iPt); // 0-60%
	h_mInvMassFxtPhiYieldsSE[fxtPhiYieldsMinBiasKey] = (TH1F*)h_mInvMassFxtPhiYieldsSE[fxtPhiYieldsKey]->Clone(fxtPhiYieldsMinBiasKey.c_str());
      }
      if(iCent > 2 && iCent <=8)
      {
	std::string fxtPhiYieldsMinBiasKey = Form("h_mInvMassFxtPhiYieldsSECent9Pt%d",iPt);
	h_mInvMassFxtPhiYieldsSE[fxtPhiYieldsMinBiasKey]->Add(h_mInvMassFxtPhiYieldsSE[fxtPhiYieldsKey],1.0);
      }
    }
  }

  string inputFileME = Form("../../data/PhiMesonAnalyzer/%s/file_FlowPhiME_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPutME = TFile::Open(inputFileME.c_str());
  if(!file_InPutME->IsOpen()) cout << "inputFileME: " << inputFileME.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFileME.c_str() << endl;
  for(int iCent = 0; iCent < mNumCentBinFxtPhiYileds; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinFxtPhiFlow; ++iPt)
    {
      std::string fxtPhiYieldsKey = Form("h_mInvMassFxtPhiYieldsMECent%dPt%d",iCent,iPt);
      h_mInvMassFxtPhiYieldsME[fxtPhiYieldsKey] = (TH1F*)file_InPutME->Get(fxtPhiYieldsKey.c_str())->Clone();
    }
  }

  {
    TCanvas *c_phiYields = new TCanvas("c_phiYields","c_phiYields",10,10,1000,400);
    c_phiYields->Divide(5,2);
    for(int iPad = 0; iPad < mNumCentBinFxtPhiYileds; ++iPad)
    {
      c_phiYields->cd(iPad+1)->SetLeftMargin(0.1);
      c_phiYields->cd(iPad+1)->SetRightMargin(0.1);
      c_phiYields->cd(iPad+1)->SetBottomMargin(0.1);
      c_phiYields->cd(iPad+1)->SetTicks(1,1);
      c_phiYields->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/PhiMesonAnalyzer/%s/phiYields_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_phiYields->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentBinFxtPhiYileds+1; ++iCent)
    {
      figName = Form("../../figures/PhiMesonAnalyzer/%s/phiYields_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      for(int iPt = 0; iPt < mNumPtBinFxtPhiFlow; ++iPt)
      {
	c_phiYields->cd(iPt+1)->Clear();
	c_phiYields->cd(iPt+1); 
	std::string fxtPhiYieldsKey = Form("h_mInvMassFxtPhiYieldsSECent%dPt%d",iCent,iPt);
	h_mInvMassFxtPhiYieldsSE[fxtPhiYieldsKey]->Draw("hE");
      }
      c_phiYields->Update();
      c_phiYields->Print(figName.c_str());
    }
    figName = Form("../../figures/PhiMesonAnalyzer/%s/phiYields_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_phiYields->Print(figName.c_str());
  }

}
