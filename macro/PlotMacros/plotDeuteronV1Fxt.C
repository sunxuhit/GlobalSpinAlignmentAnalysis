#include <string>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphErrors.h"

#include "../../Utility/include/StSpinAlignmentCons.h"
#include "../../Utility/include/draw.h"

void plotDeuteronV1Fxt()
{
  gStyle->SetOptStat(0);

  string inputFile = "../../data/EventPlaneMaker/Fxt3p85GeV_2018/file_ChargedFlow_Fxt3p85GeV_2018.root";
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;
  TProfile *p_mMixSubEpDeuV1 = (TProfile*)file_InPut->Get("p_mMixSubEpDeuV1");

  std::vector<double> vecYCms, vecV1Val, vecV1Err;
  vecYCms.clear();
  vecV1Val.clear();
  vecV1Err.clear();
  for(int iPoints = 0; iPoints < p_mMixSubEpDeuV1->GetNbinsX(); ++iPoints)
  {
    double yCms = p_mMixSubEpDeuV1->GetBinCenter(iPoints+1);
    double v1Val = p_mMixSubEpDeuV1->GetBinContent(iPoints+1);
    double v1Err = p_mMixSubEpDeuV1->GetBinError(iPoints+1);
    cout << "yCms = " << yCms << ", v1Val = " << v1Val << ", v1Err = " << v1Err << endl;
    if(v1Err > 0)
    {
      vecYCms.push_back(-1.0*yCms);
      vecV1Val.push_back(-1.0*v1Val);
      vecV1Err.push_back(v1Err);
    }
  }

  TGraphErrors *g_mMixSubEpDeuV1 = new TGraphErrors();
  for(int iPoints = 0; iPoints < vecYCms.size(); ++iPoints)
  {
    g_mMixSubEpDeuV1->SetPoint(iPoints,vecYCms[iPoints],vecV1Val[iPoints]);
    g_mMixSubEpDeuV1->SetPointError(iPoints,0.0,vecV1Err[iPoints]);
  }

  inputFile = "../../data/EventPlaneMaker/Fxt3p85GeV_2018/v1y.root";
  TFile *file_InPutXH = TFile::Open(inputFile.c_str());
  TGraphErrors *g_DeuV1XH = (TGraphErrors*)file_InPutXH->Get("v1y_d_1040");

  TCanvas *c_DeuteronV1 = new TCanvas("c_DeuteronV1","c_DeuteronV1",10,10,800,800);
  c_DeuteronV1->cd()->SetLeftMargin(0.15);
  c_DeuteronV1->cd()->SetRightMargin(0.15);
  c_DeuteronV1->cd()->SetBottomMargin(0.15);
  c_DeuteronV1->cd()->SetTicks(1,1);
  c_DeuteronV1->cd()->SetGrid(0,0);
  p_mMixSubEpDeuV1->SetStats(0);
  p_mMixSubEpDeuV1->GetYaxis()->SetRangeUser(-0.8,0.8);
  p_mMixSubEpDeuV1->GetXaxis()->SetNdivisions(505);
  p_mMixSubEpDeuV1->GetYaxis()->SetNdivisions(505);
  p_mMixSubEpDeuV1->GetXaxis()->SetTitle("y_{c.m.}");
  p_mMixSubEpDeuV1->GetYaxis()->SetTitle("v_{1}");
  p_mMixSubEpDeuV1->GetXaxis()->CenterTitle();
  p_mMixSubEpDeuV1->GetYaxis()->CenterTitle();
  p_mMixSubEpDeuV1->SetMarkerStyle(20);
  p_mMixSubEpDeuV1->SetMarkerSize(1.4);
  p_mMixSubEpDeuV1->SetMarkerColor(2);
  p_mMixSubEpDeuV1->Draw("pEX0");

  g_mMixSubEpDeuV1->SetMarkerStyle(24);
  g_mMixSubEpDeuV1->SetMarkerSize(1.4);
  g_mMixSubEpDeuV1->SetMarkerColor(2);
  g_mMixSubEpDeuV1->Draw("pE same");

  g_DeuV1XH->SetMarkerStyle(25);
  g_DeuV1XH->SetMarkerSize(1.4);
  g_DeuV1XH->SetMarkerColor(4);
  g_DeuV1XH->Draw("pE same");
 
  PlotLine(0.0, 0.0, -0.8, 0.8, 1, 2, 2);
  PlotLine(-1.0, 1.0, 0.0, 0.0, 1, 2, 2);

  TLegend *leg = new TLegend(0.15,0.7,0.5,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0,0.0);
  leg->AddEntry(p_mMixSubEpDeuV1,"d (no Eff. Corr.)","P");
  leg->AddEntry(g_mMixSubEpDeuV1,"d flipped","P");
  leg->AddEntry(g_DeuV1XH,"Xionghong's Results","P");
  leg->Draw("same");

  plotTopLegend((char*)"Au+Au 10\%-40\%",0.55,0.25,0.04,1,0.0,42,1);
  plotTopLegend((char*)"#sqrt{s_{NN}} = 3.0 GeV",0.55,0.20,0.04,1,0.0,42,1);

  c_DeuteronV1->SaveAs("../../figures/EventPlaneMaker/Fxt3p85GeV_2018/DeuteronV1_Fxt3p85GeV_2018.pdf");
}
