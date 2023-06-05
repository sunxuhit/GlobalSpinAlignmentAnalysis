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

double funcMixEp1Res1(double *x_val, double *par)
{
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  double res1Full = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return res1Full;
}

double funcMixEp1Res2(double *x_val, double *par)
{
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi())/(2.0*TMath::Sqrt(2.0));
  double besselOneHalf = TMath::Sqrt(2.0*arg/TMath::Pi()) * TMath::SinH(arg)/arg;
  double besselThreeHalf = TMath::Sqrt(2.0*arg/TMath::Pi()) * (TMath::CosH(arg)/arg - TMath::SinH(arg)/(arg*arg) );

  double res12Sub = norm * chi * TMath::Exp(-1.0*arg) * (besselOneHalf + besselThreeHalf);

  return res12Sub;
}

void plotMixResolution(int beamType = 2)
{
  gStyle->SetOptStat(0);
  const int mNumCentrality = 9;
  const int mNumEpGroup    = 6;
  const double mCentrality[mNumCentrality] = {75.0, 65.0, 55.0, 45.0, 35.0, 25.0, 15.0, 7.5, 2.5};

  string inputFileXH = "../../data/EventPlaneMaker/Fxt3p85GeV_2018/res_3GeV.root";
  TFile *file_InPutXH = TFile::Open(inputFileXH.c_str());
  TGraphErrors *g_resXH = (TGraphErrors*)file_InPutXH->Get("Graph");

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
  double mMixSubEp1Res1Val[mNumCentrality][mNumEpGroup];
  double mMixSubEp1Res1Err[mNumCentrality][mNumEpGroup];
  double mMixSubEp1Res2Val[mNumCentrality][mNumEpGroup];
  double mMixSubEp1Res2Err[mNumCentrality][mNumEpGroup];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      mMixSubEp1Res1Val[iCent][iGrp] = 0.0;
      mMixSubEp1Res1Err[iCent][iGrp] = 0.0;
      mMixSubEp1Res2Val[iCent][iGrp] = 0.0;
      mMixSubEp1Res2Err[iCent][iGrp] = 0.0;
    }
  }

  TGraphErrors *g_mMixSubEp1Res1[mNumEpGroup]; // 1st EP Res
  TGraphErrors *g_mMixSubEp1Res2[mNumEpGroup]; // 2nd EP Res
  for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
  {
    std::string grapName = Form("g_mMixSubEp1ResGrp%d",iGrp);
    g_mMixSubEp1Res1[iGrp] = new TGraphErrors();
    g_mMixSubEp1Res1[iGrp]->SetName(grapName.c_str());

    grapName = Form("g_mMixSubEp1Res2Grp%d",iGrp);
    g_mMixSubEp1Res2[iGrp] = new TGraphErrors();
    g_mMixSubEp1Res2[iGrp]->SetName(grapName.c_str());
  }

  TF1 *f_MixEp1Res1 = new TF1("f_MixEp1Res1",funcMixEp1Res1,0,10,0);
  TF1 *f_MixEp1Res2 = new TF1("f_MixEp1Res2",funcMixEp1Res2,0,10,0);
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      double valRes1Sub  = -999.9;
      double errRes1Sub  = 1.0;
      double valRes1Raw  = valRes1Temp[iCent][iGrp];
      double errRes1Raw  = errRes1Temp[iCent][iGrp];

      double valRes12Sub = -999.9;
      double errRes12Sub = 1.0;

      if(valRes1Raw > 0 && valRes1Raw < 1.0)
      {
	valRes1Sub = TMath::Sqrt(valRes1Raw);
	errRes1Sub = errRes1Raw/(2.0*valRes1Sub);

	double chi1Sub    = f_MixEp1Res1->GetX(valRes1Sub); // calculate 2nd EP Res
	valRes12Sub       = f_MixEp1Res2->Eval(chi1Sub);
	double errChi1Sub = errRes1Sub/f_MixEp1Res1->Derivative(chi1Sub); // error propagation
	errRes12Sub       = f_MixEp1Res2->Derivative(chi1Sub)*errChi1Sub;
      }
      mMixSubEp1Res1Val[iCent][iGrp] = valRes1Sub;
      mMixSubEp1Res1Err[iCent][iGrp] = errRes1Sub;
      g_mMixSubEp1Res1[iGrp]->SetPoint(iCent, mCentrality[iCent], mMixSubEp1Res1Val[iCent][iGrp]);
      g_mMixSubEp1Res1[iGrp]->SetPointError(iCent, 0.0, mMixSubEp1Res1Err[iCent][iGrp]);

      mMixSubEp1Res2Val[iCent][iGrp] = valRes12Sub;
      mMixSubEp1Res2Err[iCent][iGrp] = errRes12Sub;
      g_mMixSubEp1Res2[iGrp]->SetPoint(iCent, mCentrality[iCent], mMixSubEp1Res2Val[iCent][iGrp]);
      g_mMixSubEp1Res2[iGrp]->SetPointError(iCent, 0.0, mMixSubEp1Res2Err[iCent][iGrp]);
    }
  }
  //--------------------------------------------------------------------------------

  TCanvas *c_MixSubEp1Res = new TCanvas("c_MixSubEp1Res","c_MixSubEp1Res",10,10,800,800);
  c_MixSubEp1Res->cd()->SetLeftMargin(0.15);
  c_MixSubEp1Res->cd()->SetRightMargin(0.15);
  c_MixSubEp1Res->cd()->SetBottomMargin(0.15);
  c_MixSubEp1Res->cd()->SetTicks(1,1);
  c_MixSubEp1Res->cd()->SetGrid(0,0);
  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.5,99.5);
  for(int iBin = 0; iBin < 100; ++iBin)
  {
    h_frame->SetBinContent(iBin+1,-10.0);
    h_frame->SetBinError(iBin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(-0.5,80.0);
  h_frame->GetYaxis()->SetRangeUser(-0.05,0.8);
  h_frame->GetXaxis()->SetNdivisions(505);
  h_frame->GetYaxis()->SetNdivisions(505);
  h_frame->GetXaxis()->SetTitle("Centrality (%)");
  h_frame->GetYaxis()->SetTitle("1st EP Resolution");
  h_frame->GetXaxis()->CenterTitle();
  h_frame->GetYaxis()->CenterTitle();

  std::string figName = Form("../../figures/EventPlaneMaker/%s/EpdMixEp1ResComp_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_MixSubEp1Res->Print(figName.c_str());
  figName = Form("../../figures/EventPlaneMaker/%s/EpdMixEp1ResComp_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());

  h_frame->DrawCopy("pE");
  g_resXH->SetMarkerStyle(24);
  g_resXH->SetMarkerSize(1.6);
  g_resXH->SetMarkerColor(2);
  g_resXH->Draw("pE same");

  g_mMixSubEp1Res1[0]->SetMarkerStyle(24);
  g_mMixSubEp1Res1[0]->SetMarkerSize(1.6);
  g_mMixSubEp1Res1[0]->SetMarkerColor(4);
  g_mMixSubEp1Res1[0]->Draw("pE same");
  g_mMixSubEp1Res1[3]->SetMarkerStyle(25);
  g_mMixSubEp1Res1[3]->SetMarkerSize(1.4);
  g_mMixSubEp1Res1[3]->SetMarkerColor(4);
  g_mMixSubEp1Res1[3]->Draw("pE same");

  g_mMixSubEp1Res2[0]->SetMarkerStyle(20);
  g_mMixSubEp1Res2[0]->SetMarkerSize(1.6);
  g_mMixSubEp1Res2[0]->SetMarkerColor(4);
  g_mMixSubEp1Res2[0]->Draw("pE same");
  g_mMixSubEp1Res2[3]->SetMarkerStyle(21);
  g_mMixSubEp1Res2[3]->SetMarkerSize(1.4);
  g_mMixSubEp1Res2[3]->SetMarkerColor(4);
  g_mMixSubEp1Res2[3]->Draw("pE same");

  TLegend *leg = new TLegend(0.20,0.2,0.65,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0,0.0);
  leg->AddEntry(g_resXH,"R1: Xionghong's Results","P");
  leg->AddEntry(g_mMixSubEp1Res1[0],"R1: Epd0AB vs. Epd1CD && TpcWest","P");
  leg->AddEntry(g_mMixSubEp1Res1[3],"R1: Epd1CD vs. Epd0AB && TpcWest","P");
  leg->Draw("same");

  c_MixSubEp1Res->Update();
  c_MixSubEp1Res->Print(figName.c_str());

  c_MixSubEp1Res->Update();
  c_MixSubEp1Res->Print(figName.c_str());

  figName = Form("../../figures/EventPlaneMaker/%s/EpdMixEp1ResComp_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_MixSubEp1Res->Print(figName.c_str());

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

