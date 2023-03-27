#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getEpdReCtrPar(int beamType = 1)
{
  gStyle->SetOptStat(0);
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumCentrality = 9;
  const int mNumRingsGrps  = 2;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_ReCenterPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // ReCtr Correction | x axis is runIndex, y axis is Centrality
  TProfile2D *p_mEpdQ1SideReCtrXEast[mNumVzBin]; // 1st EP
  TProfile2D *p_mEpdQ1SideReCtrYEast[mNumVzBin];
  TProfile2D *p_mEpdQ1SideReCtrXWest[mNumVzBin];
  TProfile2D *p_mEpdQ1SideReCtrYWest[mNumVzBin];
  TProfile2D *p_mEpdQ1GrpReCtrXEast[mNumVzBin][mNumRingsGrps];
  TProfile2D *p_mEpdQ1GrpReCtrYEast[mNumVzBin][mNumRingsGrps];
  TProfile2D *p_mEpdQ1GrpReCtrXWest[mNumVzBin][mNumRingsGrps];
  TProfile2D *p_mEpdQ1GrpReCtrYWest[mNumVzBin][mNumRingsGrps];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mEpdQ1SideReCtrXEastVz%d",iVz); // 1st EP
    p_mEpdQ1SideReCtrXEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mEpdQ1SideReCtrYEastVz%d",iVz);
    p_mEpdQ1SideReCtrYEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mEpdQ1SideReCtrXWestVz%d",iVz);
    p_mEpdQ1SideReCtrXWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mEpdQ1SideReCtrYWestVz%d",iVz);
    p_mEpdQ1SideReCtrYWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      proName = Form("p_mEpdQ1Grp%dReCtrXEastVz%d",iGrp,iVz); // 1st EP
      p_mEpdQ1GrpReCtrXEast[iVz][iGrp] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mEpdQ1Grp%dReCtrYEastVz%d",iGrp,iVz);
      p_mEpdQ1GrpReCtrYEast[iVz][iGrp] = (TProfile2D*)file_InPut->Get(proName.c_str());

      proName = Form("p_mEpdQ1Grp%dReCtrXWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpReCtrXWest[iVz][iGrp] = (TProfile2D*)file_InPut->Get(proName.c_str());
      proName = Form("p_mEpdQ1Grp%dReCtrYWestVz%d",iGrp,iVz);
      p_mEpdQ1GrpReCtrYWest[iVz][iGrp] = (TProfile2D*)file_InPut->Get(proName.c_str());
    }
  }

  {
    TCanvas *c_EpdQ1ReCtr = new TCanvas("c_EpdQ1ReCtr","c_EpdQ1ReCtr",10,10,800,800);
    c_EpdQ1ReCtr->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_EpdQ1ReCtr->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdQ1ReCtr->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdQ1ReCtr->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdQ1ReCtr->cd(iPad+1)->SetTicks(1,1);
      c_EpdQ1ReCtr->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdQ1ReCtr_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdQ1ReCtr->Print(figName.c_str());
    for(int iVz = 0; iVz < mNumVzBin; ++iVz)
    {
      figName = Form("../../figures/EventPlaneMaker/%s/EpdQ1ReCtr_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_EpdQ1ReCtr->cd(1)->Clear(); c_EpdQ1ReCtr->cd(1); p_mEpdQ1SideReCtrXEast[iVz]->Draw("colz");
      c_EpdQ1ReCtr->cd(2)->Clear(); c_EpdQ1ReCtr->cd(2); p_mEpdQ1SideReCtrYEast[iVz]->Draw("colz");
      c_EpdQ1ReCtr->cd(3)->Clear(); c_EpdQ1ReCtr->cd(3); p_mEpdQ1SideReCtrXWest[iVz]->Draw("colz");
      c_EpdQ1ReCtr->cd(4)->Clear(); c_EpdQ1ReCtr->cd(4); p_mEpdQ1SideReCtrYWest[iVz]->Draw("colz");
      c_EpdQ1ReCtr->Update();
      c_EpdQ1ReCtr->Print(figName.c_str());
      for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
      {
	c_EpdQ1ReCtr->cd(1)->Clear(); c_EpdQ1ReCtr->cd(1); p_mEpdQ1GrpReCtrXEast[iVz][iGrp]->Draw("colz");
	c_EpdQ1ReCtr->cd(2)->Clear(); c_EpdQ1ReCtr->cd(2); p_mEpdQ1GrpReCtrYEast[iVz][iGrp]->Draw("colz");
	c_EpdQ1ReCtr->cd(3)->Clear(); c_EpdQ1ReCtr->cd(3); p_mEpdQ1GrpReCtrXWest[iVz][iGrp]->Draw("colz");
	c_EpdQ1ReCtr->cd(4)->Clear(); c_EpdQ1ReCtr->cd(4); p_mEpdQ1GrpReCtrYWest[iVz][iGrp]->Draw("colz");
	c_EpdQ1ReCtr->Update();
	c_EpdQ1ReCtr->Print(figName.c_str());
      }
    }
    figName = Form("../../figures/EventPlaneMaker/%s/EpdQ1ReCtr_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdQ1ReCtr->Print(figName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ReCenterPar/file_EpdReCenterPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mEpdQ1SideReCtrXEast[iVz]->Write();
    p_mEpdQ1SideReCtrYEast[iVz]->Write();
    p_mEpdQ1SideReCtrXWest[iVz]->Write();
    p_mEpdQ1SideReCtrYWest[iVz]->Write();
    for(int iGrp = 0; iGrp < mNumRingsGrps; ++iGrp)
    {
      p_mEpdQ1GrpReCtrXEast[iVz][iGrp]->Write();
      p_mEpdQ1GrpReCtrYEast[iVz][iGrp]->Write();
      p_mEpdQ1GrpReCtrXWest[iVz][iGrp]->Write();
      p_mEpdQ1GrpReCtrYWest[iVz][iGrp]->Write();
    }
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mEpdEp1SideWgtEast[mNumCentrality]; // 1st weighted EP
  TH2F *h_mEpdEp1SideWgtWest[mNumCentrality];
  TH2F *h_mEpdEp1SideWgtFull[mNumCentrality];
  TH2F *h_mEpdEp1SideWgtCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1SideWgtEastCent%d",iCent); // 1st EP
    h_mEpdEp1SideWgtEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1SideWgtWestCent%d",iCent);
    h_mEpdEp1SideWgtWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1SideWgtFullCent%d",iCent);
    h_mEpdEp1SideWgtFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mEpdEp1SideWgtCorrCent%d",iCent);
    h_mEpdEp1SideWgtCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  {
    TCanvas *c_EpdEp1SideWgtDist = new TCanvas("c_EpdEp1SideWgtDist","c_EpdEp1SideWgtDist",10,10,800,800);
    c_EpdEp1SideWgtDist->Divide(2,2);
    for(int iPad = 0; iPad < 4; ++iPad)
    {
      c_EpdEp1SideWgtDist->cd(iPad+1)->SetLeftMargin(0.15);
      c_EpdEp1SideWgtDist->cd(iPad+1)->SetRightMargin(0.15);
      c_EpdEp1SideWgtDist->cd(iPad+1)->SetBottomMargin(0.15);
      c_EpdEp1SideWgtDist->cd(iPad+1)->SetTicks(1,1);
      c_EpdEp1SideWgtDist->cd(iPad+1)->SetGrid(0,0);
    }

    std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdWgtEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1SideWgtDist->Print(figName.c_str());
    for(int iCent = 0; iCent < mNumCentrality; ++iCent)
    {
      figName = Form("../../figures/EventPlaneMaker/%s/EpdWgtEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
      c_EpdEp1SideWgtDist->cd(1)->Clear(); c_EpdEp1SideWgtDist->cd(1); h_mEpdEp1SideWgtEast[iCent]->ProjectionY()->Draw();
      c_EpdEp1SideWgtDist->cd(2)->Clear(); c_EpdEp1SideWgtDist->cd(2); h_mEpdEp1SideWgtWest[iCent]->ProjectionY()->Draw();
      c_EpdEp1SideWgtDist->cd(3)->Clear(); c_EpdEp1SideWgtDist->cd(3); h_mEpdEp1SideWgtFull[iCent]->ProjectionY()->Draw();
      c_EpdEp1SideWgtDist->cd(4)->Clear(); c_EpdEp1SideWgtDist->cd(4); h_mEpdEp1SideWgtCorr[iCent]->Draw("colz");
      c_EpdEp1SideWgtDist->Update();
      c_EpdEp1SideWgtDist->Print(figName.c_str());
    }
    figName = Form("../../figures/EventPlaneMaker/%s/EpdWgtEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1SideWgtDist->Print(figName.c_str());
  }

  string outputFileWgtEp = Form("../../data/EventPlaneMaker/%s/file_EpdWgtEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Wgt EP: " << outputFileWgtEp.c_str() << endl;
  TFile *file_OutPutWgtEp = new TFile(outputFileWgtEp.c_str(),"RECREATE");
  file_OutPutWgtEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1SideWgtEast[iCent]->Write();
    h_mEpdEp1SideWgtWest[iCent]->Write();
    h_mEpdEp1SideWgtFull[iCent]->Write();
    h_mEpdEp1SideWgtCorr[iCent]->Write();
  }
  file_OutPutWgtEp->Close();
}
