#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getEpdResolution(int beamType = 0)
{
  gStyle->SetOptStat(0);
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumShiftCorr = 20;
  const int mNumCentrality = 9;
  const int mNumRingsGrps = 2;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_EpResolution_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // Shift Correction for Full EP
  TProfile *p_mEpdSubEp1SideRes;
  TProfile *p_mEpdSubEp1GrpRes[mNumRingsGrps]; // resolution of same group
  if(beamType == 0 || beamType == 1)
  {
    p_mEpdSubEp1SideRes = (TProfile*)file_InPut->Get("p_mEpdSubEp1SideRes");
  }
  if(beamType == 2)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      std::string proName = Form("p_mEpdSubEp1Grp%dRes",iGrp);
      p_mEpdSubEp1GrpRes[iGrp] = (TProfile*)file_InPut->Get(proName.c_str());
    }
  }

  if(beamType == 0 || beamType == 1)
  {
    TCanvas *c_EpdSubEp1Res = new TCanvas("c_EpdSubEp1Res","c_EpdSubEp1Res",10,10,800,800);
    c_EpdSubEp1Res->cd()->SetLeftMargin(0.15);
    c_EpdSubEp1Res->cd()->SetRightMargin(0.15);
    c_EpdSubEp1Res->cd()->SetBottomMargin(0.15);
    c_EpdSubEp1Res->cd()->SetTicks(1,1);
    c_EpdSubEp1Res->cd()->SetGrid(0,0);
    p_mEpdSubEp1SideRes->GetXaxis()->SetTitle("Centrality 9 Bins");
    p_mEpdSubEp1SideRes->GetYaxis()->SetTitle("cos(#Psi_{1}^{West}-#Psi_{1}^{East})");
    p_mEpdSubEp1SideRes->GetYaxis()->SetTitleOffset(1.25);
    p_mEpdSubEp1SideRes->GetYaxis()->SetNdivisions(505);
    p_mEpdSubEp1SideRes->GetYaxis()->SetLabelSize(0.03);
    p_mEpdSubEp1SideRes->DrawCopy("PE");

    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdSubEp1Resolution_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdSubEp1Res->SaveAs(figName.c_str());
  }
  if(beamType == 2)
  {
    TCanvas *c_EpdSubEp1Res = new TCanvas("c_EpdSubEp1Res","c_EpdSubEp1Res",10,10,800,400);
    c_EpdSubEp1Res->Divide(2,1);
    for(int iGrp = 0; iGrp < 2; ++iGrp)
    {
      c_EpdSubEp1Res->cd(iGrp+1)->SetLeftMargin(0.15);
      c_EpdSubEp1Res->cd(iGrp+1)->SetRightMargin(0.15);
      c_EpdSubEp1Res->cd(iGrp+1)->SetBottomMargin(0.15);
      c_EpdSubEp1Res->cd(iGrp+1)->SetTicks(1,1);
      c_EpdSubEp1Res->cd(iGrp+1)->SetGrid(0,0);
      p_mEpdSubEp1GrpRes[iGrp]->GetXaxis()->SetTitle("Centrality 9 Bins");
      p_mEpdSubEp1GrpRes[iGrp]->GetYaxis()->SetTitle("cos(#Psi_{1}^{West}-#Psi_{1}^{East})");
      p_mEpdSubEp1GrpRes[iGrp]->GetYaxis()->SetTitleOffset(1.25);
      p_mEpdSubEp1GrpRes[iGrp]->GetYaxis()->SetNdivisions(505);
      p_mEpdSubEp1GrpRes[iGrp]->GetYaxis()->SetLabelSize(0.03);
      p_mEpdSubEp1GrpRes[iGrp]->DrawCopy("PE");
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdSubEp1Resolution_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdSubEp1Res->SaveAs(figName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/Resolution/file_EpdEpResolution_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  if(beamType == 0 || beamType == 1)
  {
    p_mEpdSubEp1SideRes->Write();
  }
  if(beamType == 2)
  {
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      p_mEpdSubEp1GrpRes[iGrp]->Write();
    }
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mEpdEp1SideShiftEast[mNumCentrality]; // shift EP
  TH2F *h_mEpdEp1SideShiftWest[mNumCentrality];
  TH2F *h_mEpdEp1SideShiftFull[mNumCentrality]; // Qwest-QEast
  TH2F *h_mEpdEp1SideShiftCorr[mNumCentrality]; // Psi1East vs Psi1West
  TH2F *h_mEpdEp1SideShiftFullCorr[mNumCentrality];
  TH2F *h_mEpdEp1GrpShiftTrkAveEast[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpShiftTrkAveWest[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpShiftTrkAveFull[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpShiftTrkAveCorr[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpShiftEvtAveEast[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpShiftEvtAveWest[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpShiftEvtAveFull[mNumCentrality][mNumRingsGrps];
  TH2F *h_mEpdEp1GrpShiftEvtAveCorr[mNumCentrality][mNumRingsGrps];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    if(beamType == 0 || beamType == 1)
    {
      std::string histName = Form("h_mEpdEp1SideShiftEastCent%d",iCent);
      h_mEpdEp1SideShiftEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideShiftWestCent%d",iCent);
      h_mEpdEp1SideShiftWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideShiftFullCent%d",iCent);
      h_mEpdEp1SideShiftFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideShiftCorrCent%d",iCent);
      h_mEpdEp1SideShiftCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
      histName = Form("h_mEpdEp1SideShiftFullCorrCent%d",iCent);
      h_mEpdEp1SideShiftFullCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    }
    if(beamType == 2)
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	std::string histName = Form("h_mEpdEp1Grp%dShiftTrkAveEastCent%d",iGrp,iCent);
	h_mEpdEp1GrpShiftTrkAveEast[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dShiftTrkAveWestCent%d",iGrp,iCent);
	h_mEpdEp1GrpShiftTrkAveWest[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dShiftTrkAveFullCent%d",iGrp,iCent);
	h_mEpdEp1GrpShiftTrkAveFull[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dShiftTrkAveCorrCent%d",iGrp,iCent);
	h_mEpdEp1GrpShiftTrkAveCorr[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());

	histName = Form("h_mEpdEp1Grp%dShiftEvtAveEastCent%d",iGrp,iCent);
	h_mEpdEp1GrpShiftEvtAveEast[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dShiftEvtAveWestCent%d",iGrp,iCent);
	h_mEpdEp1GrpShiftEvtAveWest[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dShiftEvtAveFullCent%d",iGrp,iCent);
	h_mEpdEp1GrpShiftEvtAveFull[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
	histName = Form("h_mEpdEp1Grp%dShiftEvtAveCorrCent%d",iGrp,iCent);
	h_mEpdEp1GrpShiftEvtAveCorr[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
      }
    }
  }

  {
    TCanvas *c_EpdEp1ShiftDist = new TCanvas("c_EpdEp1ShiftDist","c_EpdEp1ShiftDist",10,10,800,1200);
    c_EpdEp1ShiftDist->Divide(2,3);
    for(int iPad = 0; iPad < 6; ++iPad)
    {
      c_EpdEp1ShiftDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetTicks(1,1);
      c_EpdEp1ShiftDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdShiftFullEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1ShiftDist->Print(figName.c_str());
    figName = Form("../../figures/EventPlaneMaker/%s/EpdShiftFullEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      if(beamType == 0 || beamType == 1)
      {
	c_EpdEp1ShiftDist->cd(1)->Clear(); c_EpdEp1ShiftDist->cd(1); h_mEpdEp1SideShiftEast[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1ShiftDist->cd(2)->Clear(); c_EpdEp1ShiftDist->cd(2); h_mEpdEp1SideShiftWest[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1ShiftDist->cd(3)->Clear(); c_EpdEp1ShiftDist->cd(3); h_mEpdEp1SideShiftFull[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1ShiftDist->cd(4)->Clear(); c_EpdEp1ShiftDist->cd(4); h_mEpdEp1SideShiftFullCorr[iCent]->ProjectionY()->DrawCopy();
	c_EpdEp1ShiftDist->cd(5)->Clear(); c_EpdEp1ShiftDist->cd(5); h_mEpdEp1SideShiftCorr[iCent]->DrawCopy("colz");
	c_EpdEp1ShiftDist->Update();
	c_EpdEp1ShiftDist->Print(figName.c_str());
      }
      if(beamType == 2)
      {
	for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	{
	  c_EpdEp1ShiftDist->cd(1)->Clear(); c_EpdEp1ShiftDist->cd(1); h_mEpdEp1GrpShiftTrkAveEast[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1ShiftDist->cd(2)->Clear(); c_EpdEp1ShiftDist->cd(2); h_mEpdEp1GrpShiftTrkAveWest[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1ShiftDist->cd(3)->Clear(); c_EpdEp1ShiftDist->cd(3); h_mEpdEp1GrpShiftTrkAveFull[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1ShiftDist->cd(5)->Clear(); c_EpdEp1ShiftDist->cd(5); h_mEpdEp1GrpShiftTrkAveCorr[iCent][iGrp]->DrawCopy("colz");
	  c_EpdEp1ShiftDist->Update();
	  c_EpdEp1ShiftDist->Print(figName.c_str());
	}
	for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
	{
	  c_EpdEp1ShiftDist->cd(1)->Clear(); c_EpdEp1ShiftDist->cd(1); h_mEpdEp1GrpShiftEvtAveEast[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1ShiftDist->cd(2)->Clear(); c_EpdEp1ShiftDist->cd(2); h_mEpdEp1GrpShiftEvtAveWest[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1ShiftDist->cd(3)->Clear(); c_EpdEp1ShiftDist->cd(3); h_mEpdEp1GrpShiftEvtAveFull[iCent][iGrp]->ProjectionY()->DrawCopy();
	  c_EpdEp1ShiftDist->cd(5)->Clear(); c_EpdEp1ShiftDist->cd(5); h_mEpdEp1GrpShiftEvtAveCorr[iCent][iGrp]->DrawCopy("colz");
	  c_EpdEp1ShiftDist->Update();
	  c_EpdEp1ShiftDist->Print(figName.c_str());
	}
      }
    }
    figName = Form("../../figures/EventPlaneMaker/%s/EpdShiftFullEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1ShiftDist->Print(figName.c_str());
  }

  string outputFileShiftEp = Form("../../data/EventPlaneMaker/%s/file_EpdShiftFullEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Shift EP: " << outputFileShiftEp.c_str() << endl;
  TFile *file_OutPutShiftEp = new TFile(outputFileShiftEp.c_str(),"RECREATE");
  file_OutPutShiftEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    if(beamType == 0 || beamType == 1)
    {
      h_mEpdEp1SideShiftEast[iCent]->Write();
      h_mEpdEp1SideShiftWest[iCent]->Write();
      h_mEpdEp1SideShiftFull[iCent]->Write();
      h_mEpdEp1SideShiftCorr[iCent]->Write();
      h_mEpdEp1SideShiftFullCorr[iCent]->Write();
    }
    if(beamType == 2)
    {
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	h_mEpdEp1GrpShiftTrkAveEast[iCent][iGrp]->Write();
	h_mEpdEp1GrpShiftTrkAveWest[iCent][iGrp]->Write();
	h_mEpdEp1GrpShiftTrkAveFull[iCent][iGrp]->Write();
	h_mEpdEp1GrpShiftTrkAveCorr[iCent][iGrp]->Write();

	h_mEpdEp1GrpShiftEvtAveEast[iCent][iGrp]->Write();
	h_mEpdEp1GrpShiftEvtAveWest[iCent][iGrp]->Write();
	h_mEpdEp1GrpShiftEvtAveFull[iCent][iGrp]->Write();
	h_mEpdEp1GrpShiftEvtAveCorr[iCent][iGrp]->Write();
      }
    }
  }
  file_OutPutShiftEp->Close();
  file_InPut->Close();
}
