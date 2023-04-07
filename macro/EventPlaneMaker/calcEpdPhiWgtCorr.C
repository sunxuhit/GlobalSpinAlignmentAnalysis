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
  const int mNumRingsGrps  = 2;

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
      h_mEpdPhiWgtEast[iCent]->DrawCopy("colz");
    }
    figName = Form("../../figures/EventPlaneMaker/%s/EpdPhiWgt_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdPhiWgt->Update();
    c_EpdPhiWgt->Print(figName.c_str());

    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      c_EpdPhiWgt->cd(iCent+1);
      h_mEpdPhiWgtWest[iCent]->SetStats(0);
      h_mEpdPhiWgtWest[iCent]->Divide(h_mEpdPhiAveWest[iCent]);
      h_mEpdPhiWgtWest[iCent]->DrawCopy("colz");
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
  // plot Event Plane Distribution
  TH2F *h_mEpdEp1SideRawEast[mNumCentrality]; // 1st raw EP
  TH2F *h_mEpdEp1SideRawWest[mNumCentrality];
  TH2F *h_mEpdEp1SideRawFull[mNumCentrality];
  TH2F *h_mEpdEp1SideRawCorr[mNumCentrality]; // Psi1East vs Psi1West
  TH2F *h_mEpdEp1GrpRawEast[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpRawWest[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpRawFull[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpRawCorr[mNumCentrality][mNumRingsGrps];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    if(beamType == 0 || beamType == 1) // ZrZr200GeV_2018 & RuRu200GeV_2018
    {
      std::string histName = Form("h_mEpdEp1SideRawEastCent%d",iCent); // 2nd EP
      h_mEpdEp1SideRawEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideRawWestCent%d",iCent);
      h_mEpdEp1SideRawWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideRawFullCent%d",iCent);
      h_mEpdEp1SideRawFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideRawCorrCent%d",iCent);
      h_mEpdEp1SideRawCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    }
    if(beamType == 2) // Fxt3p85GeV_2018
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	std::string histName = Form("h_mEpdEp1Grp%dRawEastCent%d",iGrp,iCent);
	h_mEpdEp1GrpRawEast[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dRawWestCent%d",iGrp,iCent);
	h_mEpdEp1GrpRawWest[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dRawFullCent%d",iGrp,iCent);
	h_mEpdEp1GrpRawFull[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dRawCorrCent%d",iGrp,iCent);
	h_mEpdEp1GrpRawCorr[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
      }
    }
  }

  {
    TCanvas *c_EpdEp1SideRawDist = new TCanvas("c_EpdEp1SideRawDist","c_EpdEp1SideRawDist",10,10,800,800);
    c_EpdEp1SideRawDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_EpdEp1SideRawDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdEp1SideRawDist->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdEp1SideRawDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdEp1SideRawDist->cd(iPad+1)->SetTicks(1,1);
      c_EpdEp1SideRawDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdRawEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1SideRawDist->Print(figName.c_str());
    figName = Form("../../figures/EventPlaneMaker/%s/EpdRawEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      if(beamType == 0 || beamType == 1)
      {
	c_EpdEp1SideRawDist->cd(1)->Clear(); c_EpdEp1SideRawDist->cd(1); h_mEpdEp1SideRawEast[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1SideRawDist->cd(2)->Clear(); c_EpdEp1SideRawDist->cd(2); h_mEpdEp1SideRawWest[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1SideRawDist->cd(3)->Clear(); c_EpdEp1SideRawDist->cd(3); h_mEpdEp1SideRawFull[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1SideRawDist->cd(4)->Clear(); c_EpdEp1SideRawDist->cd(4); h_mEpdEp1SideRawCorr[iCent]->DrawCopy("colz");
	c_EpdEp1SideRawDist->Update();
	c_EpdEp1SideRawDist->Print(figName.c_str());
      }
      if(beamType == 2)
      {
	for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	{
	  c_EpdEp1SideRawDist->cd(1)->Clear(); c_EpdEp1SideRawDist->cd(1); h_mEpdEp1GrpRawEast[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1SideRawDist->cd(2)->Clear(); c_EpdEp1SideRawDist->cd(2); h_mEpdEp1GrpRawWest[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1SideRawDist->cd(3)->Clear(); c_EpdEp1SideRawDist->cd(3); h_mEpdEp1GrpRawFull[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1SideRawDist->cd(4)->Clear(); c_EpdEp1SideRawDist->cd(4); h_mEpdEp1GrpRawCorr[iCent][iGrp]->DrawCopy("colz");
	  c_EpdEp1SideRawDist->Update();
	  c_EpdEp1SideRawDist->Print(figName.c_str());
	}
      }
    }
    figName = Form("../../figures/EventPlaneMaker/%s/EpdRawEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1SideRawDist->Print(figName.c_str());
  }

  string outputFileRawEp = Form("../../data/EventPlaneMaker/%s/file_EpdRawEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileRawEp.c_str() << endl;
  TFile *file_OutPutRawEp = new TFile(outputFileRawEp.c_str(),"RECREATE");
  file_OutPutRawEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    if(beamType == 0 || beamType == 1)
    {
      h_mEpdEp1SideRawEast[iCent]->Write();
      h_mEpdEp1SideRawWest[iCent]->Write();
      h_mEpdEp1SideRawFull[iCent]->Write();
      h_mEpdEp1SideRawCorr[iCent]->Write();
    }
    if(beamType == 2)
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	h_mEpdEp1GrpRawEast[iCent][iGrp]->Write();
	h_mEpdEp1GrpRawWest[iCent][iGrp]->Write();
	h_mEpdEp1GrpRawFull[iCent][iGrp]->Write();
	h_mEpdEp1GrpRawCorr[iCent][iGrp]->Write();
      }
    }
  }
  file_OutPutRawEp->Close();
  file_InPut->Close();
}
