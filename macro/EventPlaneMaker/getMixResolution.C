#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

double propMixEpResErr(double valA, double sigA, double valB, double sigB, double valC, double sigC);

void getMixResolution(int beamType = 2)
{
  gStyle->SetOptStat(0);
  const int mNumCentrality = 9;
  const int mNumEpGroup    = 6;

  string inputFile = Form("../../data/EventPlaneMaker/%s/file_EpResolution_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;
  // 0: EpdEpGrp0 vs. TpcEpEast | 1: EpdEpGrp0 vs. TpcEpWest | 2: EpdEpGrp1 vs. TpcEpEast
  // 3: EpdEpGrp1 vs. TpcEpWest | 4: EpdEpGrp0 vs. EpdEpGrp1 | 5: TpcEpEast vs. TpcEpWest
  TProfile *p_mMixSubEp1Res[mNumEpGroup]; // 1st EP
  for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
  {
    std::string proName = Form("p_mMixSubEp1ResGrp%d",iGrp); // 1st EP
    p_mMixSubEp1Res[iGrp] = (TProfile*)file_InPut->Get(proName.c_str());
  }

  //--------------------------------------------------------------------------------
  // calculate 3-sub EP resolution
  double valMixEp[mNumCentrality][mNumEpGroup];
  double errMixEp[mNumCentrality][mNumEpGroup];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      valMixEp[iCent][iGrp] = p_mMixSubEp1Res[iGrp]->GetBinContent(iCent+1);
      errMixEp[iCent][iGrp] = p_mMixSubEp1Res[iGrp]->GetBinError(iCent+1);
    }
  }

  // double subRes1Grp0 = TMath::Cos(Psi1EpdGrp0-Psi1EpdGrp1)*TMath::Cos(Psi1EpdGrp0-Psi1TpcWest)/TMath::Cos(Psi1EpdGrp1-Psi1TpcWest);
  // double subRes1Grp1 = TMath::Cos(Psi1EpdGrp0-Psi1EpdGrp1)*TMath::Cos(Psi1EpdGrp0-Psi1TpcEast)/TMath::Cos(Psi1EpdGrp1-Psi1TpcEast);
  // double subRes1Grp2 = TMath::Cos(Psi1EpdGrp0-Psi1TpcEast)*TMath::Cos(Psi1EpdGrp0-Psi1TpcWest)/TMath::Cos(Psi1TpcEast-Psi1TpcWest);
  // double subRes1Grp3 = TMath::Cos(Psi1EpdGrp1-Psi1EpdGrp0)*TMath::Cos(Psi1EpdGrp1-Psi1TpcWest)/TMath::Cos(Psi1EpdGrp0-Psi1TpcWest);
  // double subRes1Grp4 = TMath::Cos(Psi1EpdGrp1-Psi1EpdGrp0)*TMath::Cos(Psi1EpdGrp1-Psi1TpcEast)/TMath::Cos(Psi1EpdGrp0-Psi1TpcEast);
  // double subRes1Grp5 = TMath::Cos(Psi1EpdGrp1-Psi1TpcEast)*TMath::Cos(Psi1EpdGrp1-Psi1TpcWest)/TMath::Cos(Psi1TpcEast-Psi1TpcWest);

  // 0: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpWest (default) | 1: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpEast | 2: EpdEpGrp0 vs. TpcEpEast && TpcEpWest
  // 3: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpWest (mainSys) | 4: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpEast | 5: EpdEpGrp1 vs. TpcEpEast && TpcEpWest
  double valRes1Temp[mNumCentrality][mNumEpGroup];
  double errRes1Temp[mNumCentrality][mNumEpGroup];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      valRes1Temp[iCent][0] = valMixEp[iCent][4]*valMixEp[iCent][1]/valMixEp[iCent][3];
      valRes1Temp[iCent][1] = valMixEp[iCent][4]*valMixEp[iCent][0]/valMixEp[iCent][2];
      valRes1Temp[iCent][2] = valMixEp[iCent][0]*valMixEp[iCent][1]/valMixEp[iCent][5];
      valRes1Temp[iCent][3] = valMixEp[iCent][4]*valMixEp[iCent][3]/valMixEp[iCent][1];
      valRes1Temp[iCent][4] = valMixEp[iCent][4]*valMixEp[iCent][2]/valMixEp[iCent][0];
      valRes1Temp[iCent][5] = valMixEp[iCent][2]*valMixEp[iCent][3]/valMixEp[iCent][5];
      errRes1Temp[iCent][0] = propMixEpResErr(valMixEp[iCent][4],errMixEp[iCent][4],valMixEp[iCent][1],errMixEp[iCent][1],valMixEp[iCent][3],errMixEp[iCent][3]);
      errRes1Temp[iCent][1] = propMixEpResErr(valMixEp[iCent][4],errMixEp[iCent][4],valMixEp[iCent][0],errMixEp[iCent][0],valMixEp[iCent][2],errMixEp[iCent][2]);
      errRes1Temp[iCent][2] = propMixEpResErr(valMixEp[iCent][0],errMixEp[iCent][0],valMixEp[iCent][1],errMixEp[iCent][1],valMixEp[iCent][5],errMixEp[iCent][5]);
      errRes1Temp[iCent][3] = propMixEpResErr(valMixEp[iCent][4],errMixEp[iCent][4],valMixEp[iCent][3],errMixEp[iCent][3],valMixEp[iCent][1],errMixEp[iCent][1]);
      errRes1Temp[iCent][4] = propMixEpResErr(valMixEp[iCent][4],errMixEp[iCent][4],valMixEp[iCent][2],errMixEp[iCent][2],valMixEp[iCent][0],errMixEp[iCent][0]);
      errRes1Temp[iCent][5] = propMixEpResErr(valMixEp[iCent][2],errMixEp[iCent][2],valMixEp[iCent][3],errMixEp[iCent][3],valMixEp[iCent][5],errMixEp[iCent][5]);
    }
  }

  // 0: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpWest (default) | 1: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpEast | 2: EpdEpGrp0 vs. TpcEpEast && TpcEpWest
  // 3: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpWest (mainSys) | 4: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpEast | 5: EpdEpGrp1 vs. TpcEpEast && TpcEpWest
  TH1F *h_mMixSubEp1Res[mNumEpGroup]; // 1st EP
  for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
  {
    std::string histName = Form("h_mMixSubEp1ResGrp%d",iGrp); // 1st EP
    h_mMixSubEp1Res[iGrp] = new TH1F(histName.c_str(),histName.c_str(),mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  }
  double mMixSubEp1ResVal[mNumCentrality][mNumEpGroup];
  double mMixSubEp1ResErr[mNumCentrality][mNumEpGroup];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      mMixSubEp1ResVal[iCent][iGrp]  = 0.0;
      mMixSubEp1ResErr[iCent][iGrp]  = 0.0;
    }
  }
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      double valRes1Sub  = -999.9;
      double errRes1Sub  = 1.0;
      double valRes1Raw = valRes1Temp[iCent][iGrp];
      double errRes1Raw = errRes1Temp[iCent][iGrp];

      if(valRes1Raw > 0)
      {
	valRes1Sub = TMath::Sqrt(valRes1Raw);
	errRes1Sub = errRes1Raw/(2.0*valRes1Sub);
      }
      mMixSubEp1ResVal[iCent][iGrp] = valRes1Sub;
      mMixSubEp1ResErr[iCent][iGrp] = errRes1Sub;
      h_mMixSubEp1Res[iGrp]->SetBinContent(iCent+1,mMixSubEp1ResVal[iCent][iGrp]);
      h_mMixSubEp1Res[iGrp]->SetBinError(iCent+1,mMixSubEp1ResErr[iCent][iGrp]);
    }
  }
  //--------------------------------------------------------------------------------

  TCanvas *c_MixSubEp1Res = new TCanvas("c_MixSubEp1Res","c_MixSubEp1Res",10,10,1200,800);
  c_MixSubEp1Res->Divide(3,2);
  for(int iPad = 0; iPad < 6; ++iPad)
  {
    c_MixSubEp1Res->cd(iPad+1)->SetLeftMargin(0.15);
    c_MixSubEp1Res->cd(iPad+1)->SetRightMargin(0.15);
    c_MixSubEp1Res->cd(iPad+1)->SetBottomMargin(0.15);
    c_MixSubEp1Res->cd(iPad+1)->SetTicks(1,1);
    c_MixSubEp1Res->cd(iPad+1)->SetGrid(0,0);
  }
  std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdSubEp1Resolution_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_MixSubEp1Res->Print(figName.c_str());
  figName = Form("../../figures/EventPlaneMaker/%s/EpdSubEp1Resolution_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
  {
    c_MixSubEp1Res->cd(iGrp+1)->Clear(); 
    c_MixSubEp1Res->cd(iGrp+1); 
    // 0: EpdEpGrp0 vs. TpcEpEast | 1: EpdEpGrp0 vs. TpcEpWest | 2: EpdEpGrp1 vs. TpcEpEast
    // 3: EpdEpGrp1 vs. TpcEpWest | 4: EpdEpGrp0 vs. EpdEpGrp1 | 5: TpcEpEast vs. TpcEpWest
    p_mMixSubEp1Res[iGrp]->GetXaxis()->SetTitle("Centrality 9 Bins");
    if(iGrp == 0) p_mMixSubEp1Res[iGrp]->GetYaxis()->SetTitle("<cos(#Psi_{1}^{Epd0AB}-#Psi_{1}^{TpcEast})>");
    if(iGrp == 1) p_mMixSubEp1Res[iGrp]->GetYaxis()->SetTitle("<cos(#Psi_{1}^{Epd0AB}-#Psi_{1}^{TpcWest})>");
    if(iGrp == 2) p_mMixSubEp1Res[iGrp]->GetYaxis()->SetTitle("<cos(#Psi_{1}^{Epd1CD}-#Psi_{1}^{TpcEast})>");
    if(iGrp == 3) p_mMixSubEp1Res[iGrp]->GetYaxis()->SetTitle("<cos(#Psi_{1}^{Epd1CD}-#Psi_{1}^{TpcWest})>");
    if(iGrp == 4) p_mMixSubEp1Res[iGrp]->GetYaxis()->SetTitle("<cos(#Psi_{1}^{Epd0AB}-#Psi_{1}^{Epd1CD})>");
    if(iGrp == 5) p_mMixSubEp1Res[iGrp]->GetYaxis()->SetTitle("<cos(#Psi_{1}^{TpcEast}-#Psi_{1}^{TpcWest})>");
    p_mMixSubEp1Res[iGrp]->DrawCopy();
  }
  c_MixSubEp1Res->Update();
  c_MixSubEp1Res->Print(figName.c_str());

  for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
  {
    c_MixSubEp1Res->cd(iGrp+1)->Clear(); 
    c_MixSubEp1Res->cd(iGrp+1); 
    // 0: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpWest (default) | 1: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpEast | 2: EpdEpGrp0 vs. TpcEpEast && TpcEpWest
    // 3: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpWest (mainSys) | 4: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpEast | 5: EpdEpGrp1 vs. TpcEpEast && TpcEpWest
    if(iGrp == 0) h_mMixSubEp1Res[iGrp]->SetTitle("Epd0AB vs. Epd1CD & TpcWest");
    if(iGrp == 1) h_mMixSubEp1Res[iGrp]->SetTitle("Epd0AB vs. Epd1CD & TpcEast");
    if(iGrp == 2) h_mMixSubEp1Res[iGrp]->SetTitle("Epd0AB vs. TpcEast& TpcWest");
    if(iGrp == 3) h_mMixSubEp1Res[iGrp]->SetTitle("Epd1CD vs. Epd0AB & TpcWest");
    if(iGrp == 4) h_mMixSubEp1Res[iGrp]->SetTitle("Epd1CD vs. Epd0AB & TpcEast");
    if(iGrp == 5) h_mMixSubEp1Res[iGrp]->SetTitle("Epd1CD vs. TpcEast & TpcWest");
    h_mMixSubEp1Res[iGrp]->GetXaxis()->SetTitle("Centrality 9 Bins");
    h_mMixSubEp1Res[iGrp]->DrawCopy();
  }
  c_MixSubEp1Res->Update();
  c_MixSubEp1Res->Print(figName.c_str());

  figName = Form("../../figures/EventPlaneMaker/%s/EpdSubEp1Resolution_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_MixSubEp1Res->Print(figName.c_str());

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/Resolution/file_MixEpResolution_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
  {
    p_mMixSubEp1Res[iGrp]->Write();
  }
  file_OutPut->Close();
  file_InPut->Close();
}

double propMixEpResErr(double valA, double sigA, double valB, double sigB, double valC, double sigC)
{ // return the error of valA*valB/valC
  double errA = sigA/valA;
  double errB = sigB/valB;
  double errC = sigC/valC;
  double valAxB = valA*valB;
  double sigAxB = valAxB*TMath::Sqrt(errA*errA+errB*errB); // valA*valB*sqrt((sigA/valA)^2+(sigB/valB)^2)
  double errAB = sigAxB/valAxB;

  double sigABdivC = valAxB/valC*TMath::Sqrt(errAB*errAB+errC*errC); // valAB/valC*sqrt((sigAB/valAB)^2+(sigC/valC)^2);

  return sigABdivC;
}

