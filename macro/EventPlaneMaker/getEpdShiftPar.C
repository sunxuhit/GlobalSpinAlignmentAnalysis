#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getEpdShiftPar(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumShiftCorr = 20;
  const int mNumCentrality = 9;
  const int mNumRingsGrps  = 2;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_ShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // Shift Correction for East/West | 0 = vertex neg/pos | 1 = shift correction harmonics
  TProfile2D *p_mEpdQ1SideShiftCosEast[mNumVzBin][mNumShiftCorr]; // 1st EP
  TProfile2D *p_mEpdQ1SideShiftSinEast[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mEpdQ1SideShiftCosWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mEpdQ1SideShiftSinWest[mNumVzBin][mNumShiftCorr];
  TProfile2D *p_mEpdQ1GrpShiftCosEast[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
  TProfile2D *p_mEpdQ1GrpShiftSinEast[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
  TProfile2D *p_mEpdQ1GrpShiftCosWest[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
  TProfile2D *p_mEpdQ1GrpShiftSinWest[mNumVzBin][mNumShiftCorr][mNumRingsGrps];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      if(beamType == 0 || beamType == 1)
      {
	std::string proName = Form("p_mEpdQ1SideShiftCos%dEastVz%d",iShift,iVz); // 1st EP
	p_mEpdQ1SideShiftCosEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
	proName = Form("p_mEpdQ1SideShiftSin%dEastVz%d",iShift,iVz);
	p_mEpdQ1SideShiftSinEast[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());

	proName = Form("p_mEpdQ1SideShiftCos%dWestVz%d",iShift,iVz);
	p_mEpdQ1SideShiftCosWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
	proName = Form("p_mEpdQ1SideShiftSin%dWestVz%d",iShift,iVz);
	p_mEpdQ1SideShiftSinWest[iVz][iShift] = (TProfile2D*)file_InPut->Get(proName.c_str());
      }
      if(beamType == 2)
      {
	for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	{
	  std::string proName = Form("p_mEpdQ1Grp%dShiftCos%dEastVz%d",iGrp,iShift,iVz); // 1st EP
	  p_mEpdQ1GrpShiftCosEast[iVz][iShift][iGrp] = (TProfile2D*)file_InPut->Get(proName.c_str());
	  proName = Form("p_mEpdQ1Grp%dShiftSin%dEastVz%d",iGrp,iShift,iVz);
	  p_mEpdQ1GrpShiftSinEast[iVz][iShift][iGrp] = (TProfile2D*)file_InPut->Get(proName.c_str());

	  proName = Form("p_mEpdQ1Grp%dShiftCos%dWestVz%d",iGrp,iShift,iVz);
	  p_mEpdQ1GrpShiftCosWest[iVz][iShift][iGrp] = (TProfile2D*)file_InPut->Get(proName.c_str());
	  proName = Form("p_mEpdQ1Grp%dShiftSin%dWestVz%d",iGrp,iShift,iVz);
	  p_mEpdQ1GrpShiftSinWest[iVz][iShift][iGrp] = (TProfile2D*)file_InPut->Get(proName.c_str());
	}
      }
    }
  }

  {
    TCanvas *c_EpdQ1Shift = new TCanvas("c_EpdQ1Shift","c_EpdQ1Shift",10,10,800,800);
    c_EpdQ1Shift->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_EpdQ1Shift->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdQ1Shift->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdQ1Shift->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdQ1Shift->cd(iPad+1)->SetTicks(1,1);
      c_EpdQ1Shift->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdQ1Shift_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdQ1Shift->Print(figName.c_str());
    figName = Form("../../figures/EventPlaneMaker/%s/EpdQ1Shift_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    for(int iVz = 0; iVz < mNumVzBin; ++iVz)
    {
      for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
      {
	if(beamType == 0 || beamType == 1)
	{
	  c_EpdQ1Shift->cd(1)->Clear(); c_EpdQ1Shift->cd(1); p_mEpdQ1SideShiftCosEast[iVz][iShift]->DrawCopy("colz");
	  c_EpdQ1Shift->cd(2)->Clear(); c_EpdQ1Shift->cd(2); p_mEpdQ1SideShiftSinEast[iVz][iShift]->DrawCopy("colz");
	  c_EpdQ1Shift->cd(3)->Clear(); c_EpdQ1Shift->cd(3); p_mEpdQ1SideShiftCosWest[iVz][iShift]->DrawCopy("colz");
	  c_EpdQ1Shift->cd(4)->Clear(); c_EpdQ1Shift->cd(4); p_mEpdQ1SideShiftSinWest[iVz][iShift]->DrawCopy("colz");
	  c_EpdQ1Shift->Update();
	  c_EpdQ1Shift->Print(figName.c_str());
	}
	if(beamType == 2)
	{
	  for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	  {
	    c_EpdQ1Shift->cd(1)->Clear(); c_EpdQ1Shift->cd(1); p_mEpdQ1GrpShiftCosEast[iVz][iShift][iGrp]->DrawCopy("colz");
	    c_EpdQ1Shift->cd(2)->Clear(); c_EpdQ1Shift->cd(2); p_mEpdQ1GrpShiftSinEast[iVz][iShift][iGrp]->DrawCopy("colz");
	    c_EpdQ1Shift->cd(3)->Clear(); c_EpdQ1Shift->cd(3); p_mEpdQ1GrpShiftCosWest[iVz][iShift][iGrp]->DrawCopy("colz");
	    c_EpdQ1Shift->cd(4)->Clear(); c_EpdQ1Shift->cd(4); p_mEpdQ1GrpShiftSinWest[iVz][iShift][iGrp]->DrawCopy("colz");
	    c_EpdQ1Shift->Update();
	    c_EpdQ1Shift->Print(figName.c_str());
	  }
	}
      }
    }
    figName = Form("../../figures/EventPlaneMaker/%s/EpdQ1Shift_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdQ1Shift->Print(figName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ShiftPar/file_EpdShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      if(beamType == 0 || beamType == 1)
      {
	p_mEpdQ1SideShiftCosEast[iVz][iShift]->Write();
	p_mEpdQ1SideShiftSinEast[iVz][iShift]->Write();
	p_mEpdQ1SideShiftCosWest[iVz][iShift]->Write();
	p_mEpdQ1SideShiftSinWest[iVz][iShift]->Write();
      }
      if(beamType == 2)
      {
	for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	{
	  p_mEpdQ1GrpShiftCosEast[iVz][iShift][iGrp]->Write();
	  p_mEpdQ1GrpShiftSinEast[iVz][iShift][iGrp]->Write();
	  p_mEpdQ1GrpShiftCosWest[iVz][iShift][iGrp]->Write();
	  p_mEpdQ1GrpShiftSinWest[iVz][iShift][iGrp]->Write();
	}
      }
    }
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mEpdEp1SideReCtrEast[mNumCentrality]; // 1st weighted EP
  TH2F *h_mEpdEp1SideReCtrWest[mNumCentrality];
  TH2F *h_mEpdEp1SideReCtrFull[mNumCentrality];
  TH2F *h_mEpdEp1SideReCtrCorr[mNumCentrality]; // Psi1East vs Psi1West
  TH2F *h_mEpdEp1GrpReCtrEast[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpReCtrWest[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpReCtrFull[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpReCtrCorr[mNumCentrality][mNumRingsGrps];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    if(beamType == 0 || beamType == 1)
    {
      std::string histName = Form("h_mEpdEp1SideReCtrEastCent%d",iCent); // 1st EP
      h_mEpdEp1SideReCtrEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideReCtrWestCent%d",iCent);
      h_mEpdEp1SideReCtrWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideReCtrFullCent%d",iCent);
      h_mEpdEp1SideReCtrFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideReCtrCorrCent%d",iCent);
      h_mEpdEp1SideReCtrCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    }
    if(beamType == 2)
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	std::string histName = Form("h_mEpdEp1Grp%dReCtrEastCent%d",iGrp,iCent);
	h_mEpdEp1GrpReCtrEast[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dReCtrWestCent%d",iGrp,iCent);
	h_mEpdEp1GrpReCtrWest[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dReCtrFullCent%d",iGrp,iCent);
	h_mEpdEp1GrpReCtrFull[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dReCtrCorrCent%d",iGrp,iCent);
	h_mEpdEp1GrpReCtrCorr[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
      }
    }
  }

  {
    TCanvas *c_EpdEp1SideReCtrDist = new TCanvas("c_EpdEp1SideReCtrDist","c_EpdEp1SideReCtrDist",10,10,800,800);
    c_EpdEp1SideReCtrDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_EpdEp1SideReCtrDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdEp1SideReCtrDist->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdEp1SideReCtrDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdEp1SideReCtrDist->cd(iPad+1)->SetTicks(1,1);
      c_EpdEp1SideReCtrDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdReCtrEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1SideReCtrDist->Print(figName.c_str());
    figName = Form("../../figures/EventPlaneMaker/%s/EpdReCtrEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      if(beamType == 0 || beamType == 1)
      {
	c_EpdEp1SideReCtrDist->cd(1)->Clear(); c_EpdEp1SideReCtrDist->cd(1); h_mEpdEp1SideReCtrEast[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1SideReCtrDist->cd(2)->Clear(); c_EpdEp1SideReCtrDist->cd(2); h_mEpdEp1SideReCtrWest[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1SideReCtrDist->cd(3)->Clear(); c_EpdEp1SideReCtrDist->cd(3); h_mEpdEp1SideReCtrFull[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1SideReCtrDist->cd(4)->Clear(); c_EpdEp1SideReCtrDist->cd(4); h_mEpdEp1SideReCtrCorr[iCent]->DrawCopy("colz");
	c_EpdEp1SideReCtrDist->Update();
	c_EpdEp1SideReCtrDist->Print(figName.c_str());
      }
      if(beamType == 2)
      {
	for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	{
	  c_EpdEp1SideReCtrDist->cd(1)->Clear(); c_EpdEp1SideReCtrDist->cd(1); h_mEpdEp1GrpReCtrEast[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1SideReCtrDist->cd(2)->Clear(); c_EpdEp1SideReCtrDist->cd(2); h_mEpdEp1GrpReCtrWest[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1SideReCtrDist->cd(3)->Clear(); c_EpdEp1SideReCtrDist->cd(3); h_mEpdEp1GrpReCtrFull[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1SideReCtrDist->cd(4)->Clear(); c_EpdEp1SideReCtrDist->cd(4); h_mEpdEp1GrpReCtrCorr[iCent][iGrp]->DrawCopy("colz");
	  c_EpdEp1SideReCtrDist->Update();
	  c_EpdEp1SideReCtrDist->Print(figName.c_str());
	}
      }
    }
    figName = Form("../../figures/EventPlaneMaker/%s/EpdReCtrEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1SideReCtrDist->Print(figName.c_str());
  }

  string outputFileReCtrEp = Form("../../data/EventPlaneMaker/%s/file_EpdReCtrEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileReCtrEp.c_str() << endl;
  TFile *file_OutPutReCtrEp = new TFile(outputFileReCtrEp.c_str(),"RECREATE");
  file_OutPutReCtrEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    if(beamType == 0 || beamType == 1)
    {
      h_mEpdEp1SideReCtrEast[iCent]->Write();
      h_mEpdEp1SideReCtrWest[iCent]->Write();
      h_mEpdEp1SideReCtrFull[iCent]->Write();
      h_mEpdEp1SideReCtrCorr[iCent]->Write();
    }
    if(beamType == 2)
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	h_mEpdEp1GrpReCtrEast[iCent][iGrp]->Write();
	h_mEpdEp1GrpReCtrWest[iCent][iGrp]->Write();
	h_mEpdEp1GrpReCtrFull[iCent][iGrp]->Write();
	h_mEpdEp1GrpReCtrCorr[iCent][iGrp]->Write();
      }
    }
  }
  file_OutPutReCtrEp->Close();
  file_InPut->Close();
}
