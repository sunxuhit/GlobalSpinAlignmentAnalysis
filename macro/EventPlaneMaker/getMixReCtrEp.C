#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getMixReCtrEp(int beamType = 2)
{
  gStyle->SetOptStat(0);
  const int mNumCentrality = 9;
  const int mNumEpGroup    = 6;

  const std::string str_LabelX[mNumEpGroup] = {"#Psi_{1,Epd}^{Grp0}","#Psi_{1,Epd}^{Grp0}","#Psi_{1,Epd}^{Grp1}","#Psi_{1,Epd}^{Grp1}","#Psi_{1,Epd}^{Grp0}","#Psi_{1,Tpc}^{East}"};
  const std::string str_LabelY[mNumEpGroup] = {"#Psi_{1,Tpc}^{East}","#Psi_{1,Tpc}^{West}","#Psi_{1,Tpc}^{East}","#Psi_{1,Tpc}^{West}","#Psi_{1,Epd}^{Grp1}","#Psi_{1,Tpc}^{West}"};

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_ShiftPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  TH2F *h_mMixEp1ReCtrCorr[mNumCentrality][mNumEpGroup];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      std::string histName = Form("h_mMixEp1Grp%dReCtrCent%d",iGrp,iCent);
      h_mMixEp1ReCtrCorr[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
      h_mMixEp1ReCtrCorr[iCent][iGrp]->GetXaxis()->SetTitle(str_LabelX[iGrp].c_str());
      h_mMixEp1ReCtrCorr[iCent][iGrp]->GetYaxis()->SetTitle(str_LabelY[iGrp].c_str());
    }
  }

  TCanvas *c_MixEpReCtr = new TCanvas("c_MixEpReCtr","c_MixEpReCtr",10,10,800,1200);
  c_MixEpReCtr->Divide(2,3);
  for(int iPad = 0; iPad < 6; ++iPad)
  {
    c_MixEpReCtr->cd(iPad+1)->SetLeftMargin(0.15);
    c_MixEpReCtr->cd(iPad+1)->SetRightMargin(0.15);
    c_MixEpReCtr->cd(iPad+1)->SetBottomMargin(0.15);
    c_MixEpReCtr->cd(iPad+1)->SetTicks(1,1);
    c_MixEpReCtr->cd(iPad+1)->SetGrid(0,0);
  }
  std::string figName = Form("../../figures/EventPlaneMaker/%s/MixReCtrEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_MixEpReCtr->Print(figName.c_str());
  figName = Form("../../figures/EventPlaneMaker/%s/MixReCtrEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      c_MixEpReCtr->cd(iGrp+1)->Clear(); 
      c_MixEpReCtr->cd(iGrp+1);
      h_mMixEp1ReCtrCorr[iCent][iGrp]->DrawCopy("colz");
    }
    c_MixEpReCtr->Update();
    c_MixEpReCtr->Print(figName.c_str());
  }
  figName = Form("../../figures/EventPlaneMaker/%s/MixReCtrEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_MixEpReCtr->Print(figName.c_str());


  string outputFileReCtrEp = Form("../../data/EventPlaneMaker/%s/file_MixReCtrEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of ReCtr EP: " << outputFileReCtrEp.c_str() << endl;
  TFile *file_OutPutReCtrEp = new TFile(outputFileReCtrEp.c_str(),"RECREATE");
  file_OutPutReCtrEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      h_mMixEp1ReCtrCorr[iCent][iGrp]->Write();
    }
  }
  file_OutPutReCtrEp->Close();
  file_InPut->Close();
}

