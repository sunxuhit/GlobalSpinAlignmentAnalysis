#include <string>
#include <map>

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void plotPhiAcpt(int beamType = 2)
{
  gStyle->SetOptStat(0);
  const int mNumCentBinQA = 10; // 0: 70-80%, 1: 60-70%, 2: 50-60%, 3: 40-50%, 4: 30-40%, 5: 20-30%, 6: 10-20%, 7: 5-10%, 8: 0-5%, 9: 20-60%
  const int mNumPtBinQA   = 25; // 25 bins from 0 to 5.0 GeV/c
  const int mNumRapBinQA  = 25; // 25 bins from -1.25 to 1.25

  // Acceptance Histograms
  // 0 = centrality: 9 = 20%-60%, 0-8 from StRefMultCorr
  // 1 = pt bin: 25 bins from 0 to 5.0 GeV/c
  // 2 = rapidity bin: 25 bins from -1.25 to 1.25
  std::map<string,TH1F*> h_mAcptPhiSELab;
  std::map<string,TH1F*> h_mAcptPhiSECms;
  std::map<string,TH1F*> h_mAcptPhiMELab;
  std::map<string,TH1F*> h_mAcptPhiMECms;

  string inputFileSE = Form("../../data/PhiMesonAnalyzer/%s/file_QaPhiSE_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPutSE = TFile::Open(inputFileSE.c_str());
  if(!file_InPutSE->IsOpen()) cout << "inputFileSE: " << inputFileSE.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFileSE.c_str() << endl;
  for(int iCent = 0; iCent < mNumCentBinQA; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinQA; ++iPt)
    {
      for(int iRap = 0; iRap < mNumRapBinQA; ++iRap)
      {
	std::string acptPhiLabKey = Form("h_mAcptPhiSELabCent%dPt%dRap%d",iCent,iPt,iRap);
	h_mAcptPhiSELab[acptPhiLabKey] = (TH1F*)file_InPutSE->Get(acptPhiLabKey.c_str())->Clone();
	std::string acptPhiCmsKey = Form("h_mAcptPhiSECmsCent%dPt%dRap%d",iCent,iPt,iRap);
	h_mAcptPhiSECms[acptPhiCmsKey] = (TH1F*)file_InPutSE->Get(acptPhiCmsKey.c_str())->Clone();
      }
    }
  }

  string inputFileME = Form("../../data/PhiMesonAnalyzer/%s/file_QaPhiME_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPutME = TFile::Open(inputFileME.c_str());
  if(!file_InPutME->IsOpen()) cout << "inputFileME: " << inputFileME.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFileME.c_str() << endl;
  for(int iCent = 0; iCent < mNumCentBinQA; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinQA; ++iPt)
    {
      for(int iRap = 0; iRap < mNumRapBinQA; ++iRap)
      {
	std::string acptPhiLabKey = Form("h_mAcptPhiMELabCent%dPt%dRap%d",iCent,iPt,iRap);
	h_mAcptPhiMELab[acptPhiLabKey] = (TH1F*)file_InPutME->Get(acptPhiLabKey.c_str())->Clone();
	std::string acptPhiCmsKey = Form("h_mAcptPhiMECmsCent%dPt%dRap%d",iCent,iPt,iRap);
	h_mAcptPhiMECms[acptPhiCmsKey] = (TH1F*)file_InPutME->Get(acptPhiCmsKey.c_str())->Clone();
      }
    }
  }

  {
    TCanvas *c_phiAcptLab = new TCanvas("c_phiAcptLab","c_phiAcptLab",10,10,mNumRapBinQA*100,mNumPtBinQA*100);
    c_phiAcptLab->Divide(mNumRapBinQA,mNumPtBinQA,0,0);
    for(int iPad = 0; iPad < mNumRapBinQA*mNumPtBinQA; ++iPad)
    {
      c_phiAcptLab->cd(iPad+1)->SetLeftMargin(0.0);
      c_phiAcptLab->cd(iPad+1)->SetRightMargin(0.0);
      c_phiAcptLab->cd(iPad+1)->SetTopMargin(0.0);
      c_phiAcptLab->cd(iPad+1)->SetBottomMargin(0.0);
      c_phiAcptLab->cd(iPad+1)->SetTicks(1,1);
      c_phiAcptLab->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/PhiMesonAnalyzer/%s/phiAcptLab_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_phiAcptLab->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentBinQA; ++iCent)
    {
      figName = Form("../../figures/PhiMesonAnalyzer/%s/phiAcptLab_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      for(int iPt = 0; iPt < mNumPtBinQA; ++iPt)
      {
	for(int iRap = 0; iRap < mNumRapBinQA; ++iRap)
	{
	  int padId = (mNumPtBinQA-(iPt+1))*mNumRapBinQA + (iRap+1);
	  c_phiAcptLab->cd(padId)->Clear();
	  c_phiAcptLab->cd(padId); 
	  std::string acptPhiLabKey = Form("h_mAcptPhiSELabCent%dPt%dRap%d",iCent,iPt,iRap);
	  h_mAcptPhiSELab[acptPhiLabKey]->DrawCopy("hE");
	}
      }
      c_phiAcptLab->Update();
      c_phiAcptLab->Print(figName.c_str());
    }
    figName = Form("../../figures/PhiMesonAnalyzer/%s/phiAcptLab_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_phiAcptLab->Print(figName.c_str());
  }

  {
    TCanvas *c_phiAcptCms = new TCanvas("c_phiAcptCms","c_phiAcptCms",10,10,mNumRapBinQA*100,mNumPtBinQA*100);
    c_phiAcptCms->Divide(mNumRapBinQA,mNumPtBinQA,0,0);
    for(int iPad = 0; iPad < mNumRapBinQA*mNumPtBinQA; ++iPad)
    {
      c_phiAcptCms->cd(iPad+1)->SetLeftMargin(0.0);
      c_phiAcptCms->cd(iPad+1)->SetRightMargin(0.0);
      c_phiAcptCms->cd(iPad+1)->SetTopMargin(0.0);
      c_phiAcptCms->cd(iPad+1)->SetBottomMargin(0.0);
      c_phiAcptCms->cd(iPad+1)->SetTicks(1,1);
      c_phiAcptCms->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/PhiMesonAnalyzer/%s/phiAcptCms_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_phiAcptCms->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentBinQA; ++iCent)
    {
      figName = Form("../../figures/PhiMesonAnalyzer/%s/phiAcptCms_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      for(int iPt = 0; iPt < mNumPtBinQA; ++iPt)
      {
	for(int iRap = 0; iRap < mNumRapBinQA; ++iRap)
	{
	  int padId = (mNumPtBinQA-(iPt+1))*mNumRapBinQA + (iRap+1);
	  c_phiAcptCms->cd(padId)->Clear();
	  c_phiAcptCms->cd(padId); 
	  std::string acptPhiCmsKey = Form("h_mAcptPhiSECmsCent%dPt%dRap%d",iCent,iPt,iRap);
	  h_mAcptPhiSECms[acptPhiCmsKey]->DrawCopy("hE");
	}
      }
      c_phiAcptCms->Update();
      c_phiAcptCms->Print(figName.c_str());
    }
    figName = Form("../../figures/PhiMesonAnalyzer/%s/phiAcptCms_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_phiAcptCms->Print(figName.c_str());
  }
}
