#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getMixShiftEp(int beamType = 2)
{
  gStyle->SetOptStat(0);
  const int mNumCentrality = 9;
  const int mNumEpGroup    = 6;

  const std::string str_LabelX[mNumEpGroup] = {"#Psi_{1,Epd}^{Grp0}","#Psi_{1,Epd}^{Grp0}","#Psi_{1,Epd}^{Grp1}","#Psi_{1,Epd}^{Grp1}","#Psi_{1,Epd}^{Grp0}","#Psi_{1,Tpc}^{East}"};
  const std::string str_LabelY[mNumEpGroup] = {"#Psi_{1,Tpc}^{East}","#Psi_{1,Tpc}^{West}","#Psi_{1,Tpc}^{East}","#Psi_{1,Tpc}^{West}","#Psi_{1,Epd}^{Grp1}","#Psi_{1,Tpc}^{West}"};

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_EpResolution_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  TH2F *h_mMixEp1ShiftCorr[mNumCentrality][mNumEpGroup];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      std::string histName = Form("h_mMixEp1Grp%dShiftCent%d",iGrp,iCent);
      h_mMixEp1ShiftCorr[iCent][iGrp] = (TH2F*)file_InPut->Get(histName.c_str());
      h_mMixEp1ShiftCorr[iCent][iGrp]->GetXaxis()->SetTitle(str_LabelX[iGrp].c_str());
      h_mMixEp1ShiftCorr[iCent][iGrp]->GetYaxis()->SetTitle(str_LabelY[iGrp].c_str());
    }
  }

  TCanvas *c_MixEpShift = new TCanvas("c_MixEpShift","c_MixEpShift",10,10,800,1200);
  c_MixEpShift->Divide(2,3);
  for(int iPad = 0; iPad < 6; ++iPad)
  {
    c_MixEpShift->cd(iPad+1)->SetLeftMargin(0.15);
    c_MixEpShift->cd(iPad+1)->SetRightMargin(0.15);
    c_MixEpShift->cd(iPad+1)->SetBottomMargin(0.15);
    c_MixEpShift->cd(iPad+1)->SetTicks(1,1);
    c_MixEpShift->cd(iPad+1)->SetGrid(0,0);
  }
  std::string figName = Form("../../figures/EventPlaneMaker/%s/MixShiftEp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_MixEpShift->Print(figName.c_str());
  figName = Form("../../figures/EventPlaneMaker/%s/MixShiftEp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      c_MixEpShift->cd(iGrp+1)->Clear(); 
      c_MixEpShift->cd(iGrp+1);
      h_mMixEp1ShiftCorr[iCent][iGrp]->DrawCopy("colz");
    }
    c_MixEpShift->Update();
    c_MixEpShift->Print(figName.c_str());
  }
  figName = Form("../../figures/EventPlaneMaker/%s/MixShiftEp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_MixEpShift->Print(figName.c_str());


  string outputFileShiftEp = Form("../../data/EventPlaneMaker/%s/file_MixShiftEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Shift EP: " << outputFileShiftEp.c_str() << endl;
  TFile *file_OutPutShiftEp = new TFile(outputFileShiftEp.c_str(),"RECREATE");
  file_OutPutShiftEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      h_mMixEp1ShiftCorr[iCent][iGrp]->Write();
    }
  }
  file_OutPutShiftEp->Close();
  file_InPut->Close();
}


