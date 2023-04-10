#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TF1.h>

#include "../../Utility/include/draw.h"

using namespace std;

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color);
void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color);

double rho00_theory(double *x_var, double *par)
{
  double s12 = x_var[0];
  double scurrent = par[0]; // strangeness current
  double ms = par[1]; // mass of s-quark MeV

  double mphi = 1020.0; // MeV
  double c1 = 300.0/(pow(200.0,1.0/3.0) * pow((-0.4+0.39*log(200.0*200.0)),1.0/3.0)); // 300.0/Teff[200.0]
  double gphi = 2.0*sqrt(2.0);

  double Teff = pow(s12,1.0/3.0) * pow((-0.4+0.39*log(s12*s12)),1.0/3.0);

  double denom_phifield = 27.0*pow(ms,4.0)*pow(mphi,4.0)*pow(c1,2.0)*pow(Teff,2.0);
  double numer_phifield = scurrent*pow(197.0,8.0)*1.8*1.0e+5; // <p^2>_phi = 0.18 GeV^2

  double rho00 = 1.0/3.0 + numer_phifield/denom_phifield;

  return rho00;
}


void plotRho00BESI()
{
  gStyle->SetOptDate(0);
  const int style_phi_1st = 21;
  const int color_phi_1st = kGray+2;
  const int style_phi_2nd = 29;
  const int color_phi_2nd = kRed-4;
  const int style_phi_ALICE = 30;
  const int color_phi_ALICE = kGray+1;
  const int colorDiff_phi = 0;

  const int style_Kstr = 20;
  const int color_Kstr = kAzure-9;
  const int style_Kstr_ALICE = 24;
  const int color_Kstr_ALICE = kGray+1;
  const int colorDiff_Kstr = 2;

  const float size_marker = 1.4;
  const float size_font = 0.035;
  
  //----------------------------------------------------------
  // phi-meson STAR
  //beam-energy dependence of phi-meson rho00 from STAR, pT: 1.2 - 5.4 GeV/c, 20-60%
  TFile *File_InputPhi = TFile::Open("/Users/xusun/WorkSpace/STAR/SpinAlignment/VecMesonSpinAlignment/data/PaperDraft/Nature/Phi/NewF_JHChen/rho00_stat_sys_Laxis.root");
  TGraphAsymmErrors *g_rhoPhi_2nd_stat_Laxis = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_stat");
  TGraphAsymmErrors *g_rhoPhi_2nd_sys_Laxis = (TGraphAsymmErrors*)File_InputPhi->Get("rho00_2ndEP_energy_sys");
  TGraphAsymmErrors *g_rho_2nd_fit_Laxis  = new TGraphAsymmErrors();
  for(int i_energy = 0; i_energy < g_rhoPhi_2nd_stat_Laxis->GetN(); ++i_energy) // combine stat & sys for fit
  {
    double energy, rho;
    g_rhoPhi_2nd_stat_Laxis->GetPoint(i_energy,energy,rho);
    double err_stat = g_rhoPhi_2nd_stat_Laxis->GetErrorYhigh(i_energy);
    double err_sys = g_rhoPhi_2nd_sys_Laxis->GetErrorYhigh(i_energy);
    double err_fit = TMath::Sqrt(err_stat*err_stat+err_sys*err_sys);

    g_rho_2nd_fit_Laxis->SetPoint(i_energy,energy,rho);
    g_rho_2nd_fit_Laxis->SetPointError(i_energy,0.0,0.0,err_fit,err_fit);
  }

  //----------------------------------------------------------
  TFile *File_InputKstar = TFile::Open("/Users/xusun/WorkSpace/STAR/SpinAlignment/VecMesonSpinAlignment/data/PaperDraft/Nature/Kstar/data_Kstar_rho00_sNN_Oct12_2022.root");
  // K* STAR
  //beam-energy dependence of kstar rho00 from STAR, pT: 1.0 - 1.5 GeV/c, 20-60%
  TGraphAsymmErrors *g_rhoKstar_stat       = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_stat");
  TGraphAsymmErrors *g_rhoKstar_sys        = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_sys");
  //data points from ALICE
  //K* ALICE pT: 0.8 - 1.2 GeV/c
  TGraphAsymmErrors *g_rhoKstar_ALICE_stat = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_ALICE_stat");
  TGraphAsymmErrors *g_rhoKstar_ALICE_sys  = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoKstar_ALICE_sys");
  //Phi ALICE pT: 1.0 - 5.0 GeV/c
  TGraphAsymmErrors *g_rhoPhi_ALICE_stat   = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoPhi_ALICE_stat");
  TGraphAsymmErrors *g_rhoPhi_ALICE_sys    = (TGraphAsymmErrors*)File_InputKstar->Get("g_rhoPhi_ALICE_sys");
  //----------------------------------------------------------

  TH1F *h_frame = new TH1F("h_frame","h_frame",5000,0,5000);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10.0);
    h_frame->SetBinError(i_bin+1,-10.0);
  }
  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  c_rho00->cd()->SetLogx();
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(1.0,4996.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.23,0.45);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  // h_frame->GetYaxis()->SetTitle("#rho_{00} (Out-of-Plane)");
  h_frame->GetYaxis()->SetTitle("#rho_{00}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(1.0,4996.0,1.0/3.0,1.0/3.0,1,2,2);

  // K* STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstar_stat,style_Kstr,color_Kstr,colorDiff_Kstr,size_marker-0.4);
  plotSysErrorsBox(g_rhoKstar_sys,color_Kstr+2);

  // K* ALICE
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoKstar_ALICE_stat,style_Kstr_ALICE,color_Kstr_ALICE,0,size_marker-0.4);
  plotSysErrorsBox(g_rhoKstar_ALICE_sys,color_Kstr_ALICE);

  // phi-meson STAR
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_2nd_stat_Laxis,style_phi_2nd,color_phi_2nd,colorDiff_phi,size_marker+0.2);
  plotSysErrorsBox(g_rhoPhi_2nd_sys_Laxis,color_phi_2nd);

  // phi-meson ALICE
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_rhoPhi_ALICE_stat,style_phi_ALICE,color_phi_ALICE,0,size_marker);
  plotSysErrorsBox(g_rhoPhi_ALICE_sys,color_phi_ALICE);

  // plot theory curve
  double ms = 450.0;
  int Style[3] = {kDashed,kDotted,kSolid};
  int Color[3] = {kBlue,kMagenta,kRed};
  TF1 *f_rho00_Laxis = new TF1("f_rho00_Laxis",rho00_theory,1,201,2);
  f_rho00_Laxis->SetParameter(0,1000.0);
  f_rho00_Laxis->FixParameter(1,ms);
  f_rho00_Laxis->SetRange(19.0,200.0);
  g_rho_2nd_fit_Laxis->Fit(f_rho00_Laxis,"MNR");
  double chi2_Laxis = f_rho00_Laxis->GetChisquare();
  int ndf_Laxis = f_rho00_Laxis->GetNDF();
  double chi2_ndf_Laxis = chi2_Laxis/(double)ndf_Laxis;
  double p_Laxis = TMath::Prob(chi2_Laxis,ndf_Laxis);
  cout << "chi2 for Laxis: " << chi2_Laxis << endl;
  cout << "ndf for Laxis: " << ndf_Laxis << endl;
  cout << "chi2/ndf for Laxis: " << chi2_ndf_Laxis  << ", p_Laxis: " << p_Laxis << endl;
  cout << "C^{y}_{s} = " << f_rho00_Laxis->GetParameter(0) << " +/- " << f_rho00_Laxis->GetParError(0) << endl;

  f_rho00_Laxis->SetLineColor(kRed);
  f_rho00_Laxis->SetLineStyle(1);
  f_rho00_Laxis->SetLineWidth(4);
  f_rho00_Laxis->Draw("l same");

  TF1 *f_rho00_plot = new TF1("f_rho00_plot",rho00_theory,1,4000,2);
  f_rho00_plot->FixParameter(0,f_rho00_Laxis->GetParameter(0));
  f_rho00_plot->FixParameter(1,ms);
  f_rho00_plot->SetLineStyle(2);
  f_rho00_plot->SetLineWidth(4);
  f_rho00_plot->SetRange(200.0,3000.0);
  f_rho00_plot->Draw("l same");

  // plot Legend
  // plotTopLegend((char*)"filled: STAR (Au+Au & 20\% - 60\% Centrality)",0.185,0.25,size_font,1,0.0,42,1);
  // plotTopLegend((char*)"open: ALICE (Pb+Pb & 10\% - 50\% Centrality)",0.185,0.20,size_font,1,0.0,42,1);

  plotTopLegend((char*)"STAR (Au+Au & 20\% - 60\%)",1.8,0.435,size_font,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(2,0.425,0.0,0.0,0.0,0.0,style_phi_2nd,color_phi_2nd,size_marker+0.2);
  plotTopLegend((char*)"#phi    (|y| < 1.0 & 1.2 < p_{T} < 5.4 GeV/c)",2.5,0.4225,size_font,1,0.0,42,0);

  Draw_TGAE_Point_new_Symbol(2,0.410,0.0,0.0,0.0,0.0,style_Kstr,color_Kstr,size_marker-0.4);
  plotTopLegend((char*)"K^{*0} (|y| < 1.0 & 1.0 < p_{T} < 5.0 GeV/c)",2.5,0.4075,size_font,1,0.0,42,0);

  plotTopLegend((char*)"open symbols:",0.385,0.25,size_font,1,0.0,42,1);
  plotTopLegend((char*)"ALICE (Pb+Pb & 10\% - 50\%)",0.385,0.20,size_font,1,0.0,42,1);

  // theory
  const float Gs = f_rho00_Laxis->GetParameter(0)*1.8e5*pow(197.0,8)/(pow(450.0,2)*pow(1020.0,4)*pow(138.05,4));
  const float GsErr = f_rho00_Laxis->GetParError(0)*1.8e5*pow(197.0,8)/(pow(450.0,2)*pow(1020.0,4)*pow(138.05,4));
  string leg_current_Laxis = Form("G^{(y)}_{s} = %1.2f #pm %1.2f m_{#pi}^{4}", Gs, GsErr);
  // string leg_current_Laxis = Form("C^{(y)}_{s} = %1.0f #pm %1.0f fm^{-8}",f_rho00_Laxis->GetParameter(0),f_rho00_Laxis->GetParError(0));
  // string leg_chi2_Laxis = Form("#chi^{2}/ndf: %1.1f",chi2_ndf_Laxis);
  // string leg_p_Laxis = Form("p-value: %1.3f", p_Laxis);
  // string leg_stat_Laxis = Form("#chi^{2}/ndf: %1.1f & p-value: %1.3f",chi2_ndf_Laxis,p_Laxis);
  // TLegend *leg = new TLegend(0.335,0.18,0.675,0.23);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(10);
  // leg->SetFillStyle(0);
  // leg->AddEntry(f_rho00_Laxis,leg_current_Laxis.c_str(),"l");
  // leg->AddEntry((TObject*)0,leg_stat_Laxis.c_str(),"");
  // leg->Draw("same");
  PlotLine(1.7,2.2,0.395,0.395,2,3,1);
  plotTopLegend((char*)leg_current_Laxis.c_str(),2.5,0.3925,size_font,1,0.0,42,0);

  // future measurement
  TBox *bHighMuB= new TBox(2.0, 0.23, 7.5, 0.38);
  bHighMuB->SetFillColor(kYellow);
  bHighMuB->SetFillStyle(3001);
  bHighMuB->Draw("same");

  // NICA @ JINR: 3 - 11 GeV
  TLatex *lNica;
  lNica = new TLatex(4.5, 0.31,"NICA");
  lNica->SetTextFont(42);
  lNica->SetTextSize(0.023);
  lNica->SetTextColor(kGreen+3);
  lNica->SetLineWidth(2);
  lNica->Draw();
  TArrow *aNica;
  aNica = new TArrow(3.0, 0.308, 11.0, 0.308, 0.02,"<>");
  aNica->SetFillColor(1);
  aNica->SetFillStyle(1001);
  aNica->SetLineColor(kGreen+3);
  aNica->SetLineWidth(1);
  aNica->SetLineStyle(1);
  aNica->SetAngle(29);
  aNica->Draw();
  
  //FXT : 3.0 - 7.7
  TLatex *lFxt;
  lFxt = new TLatex(2.75, 0.29,"STAR-FXT");
  lFxt->SetTextFont(42);
  lFxt->SetTextSize(0.023);
  lFxt->SetTextColor(kRed-2);
  lFxt->SetLineWidth(2);
  lFxt->Draw();
  TArrow *aFxt;
  aFxt = new TArrow(3.0, 0.288, 7.7, 0.288, 0.02,"<>");
  aFxt->SetFillColor(kRed-7);
  aFxt->SetFillStyle(1002);
  aFxt->SetLineColor(kRed-7);
  aFxt->SetLineWidth(1);
  aFxt->SetAngle(30);
  aFxt->Draw();
  
  //BES-II: 3.0 - 7.7
  TLatex *lBESII;
  lBESII = new TLatex(7.75, 0.275,"STAR-BESII");
  lBESII->SetTextFont(42);
  lBESII->SetTextSize(0.023);
  lBESII->SetTextColor(kAzure+2);
  lBESII->SetLineWidth(2);
  lBESII->Draw();
  TArrow *aBESII;
  aBESII = new TArrow(7.7, 0.282, 27.0, 0.282, 0.02,"<>");
  aBESII->SetFillColor(kAzure+2);
  aBESII->SetFillStyle(1002);
  aBESII->SetLineColor(kAzure+2);
  aBESII->SetLineWidth(1);
  aBESII->SetAngle(30);
  aBESII->Draw();
  
  // FAIR @ GSI: 2 - 4.9 GeV
  TLatex *lFair;
  lFair = new TLatex(2.40, 0.27,"FAIR");
  lFair->SetTextFont(42);
  lFair->SetTextSize(0.02283105);
  lFair->SetLineWidth(1);
  lFair->Draw();
  TArrow *aFair;
  aFair = new TArrow(2.0, 0.268, 4.9, 0.268, 0.02,"<>");
  aFair->SetFillColor(1);
  aFair->SetFillStyle(1001);
  aFair->SetLineColor(17);
  aFair->SetLineWidth(1);
  aFair->SetAngle(30);
  aFair->Draw();
  
  // HIAF @ IMP: 2 - 4 GeV
  TLatex *lHiaf;
  lHiaf = new TLatex(1.7, 0.25,"CEE@HIAF");
  lHiaf->SetTextFont(42);
  lHiaf->SetTextSize(0.022);
  lHiaf->SetTextColor(kRed-2);
  lHiaf->SetLineWidth(1);
  lHiaf->Draw();
  TArrow *aHiaf;
  aHiaf = new TArrow(2.0, 0.248, 4.0, 0.248, 0.02,"<>");
  aHiaf->SetFillColor(1);
  aHiaf->SetFillStyle(1001);
  aHiaf->SetLineColor(kRed-2);
  aHiaf->SetLineWidth(1);
  aHiaf->SetAngle(30);
  aHiaf->Draw();

  c_rho00->SaveAs("../../figures/AnalysisNote/fig_rho00FutureMeasurement.eps");
}

void plotSysErrors(TGraphAsymmErrors *g_rho, int plot_color)
{
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double err = g_rho->GetErrorYhigh(i_energy);

    PlotLine(energy*0.95,energy*1.05,rho+err,rho+err,plot_color,2,1);
    PlotLine(energy*0.95,energy*0.95,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy*1.05,energy*1.05,rho+err-0.001,rho+err,plot_color,2,1);
    PlotLine(energy*0.95,energy*1.05,rho-err,rho-err,plot_color,2,1);
    PlotLine(energy*0.95,energy*0.95,rho-err+0.001,rho-err,plot_color,2,1);
    PlotLine(energy*1.05,energy*1.05,rho-err+0.001,rho-err,plot_color,2,1);
  }
}

void plotSysErrorsBox(TGraphAsymmErrors *g_rho, int plot_color)
{
  const int nEnergy = g_rho->GetN();
  TBox *bSys[nEnergy];
  for(int i_energy = 0; i_energy < g_rho->GetN(); ++i_energy) // plot sys errors
  {
    double energy, rho;
    g_rho->GetPoint(i_energy,energy,rho);
    double errHigh = rho + g_rho->GetErrorYhigh(i_energy);
    double errLow  = rho - g_rho->GetErrorYlow(i_energy);
    if(errLow < 0.23) errLow = 0.23; // only for h_frame->GetYaxis()->SetRangeUser(0.23,0.42);

    bSys[i_energy] = new TBox(energy/1.08,errLow,energy*1.08,errHigh);
    // bSys[i_energy] = new TBox(energy-1.5,rho-err,energy+1.5,rho+err);
    bSys[i_energy]->SetFillColor(0);
    bSys[i_energy]->SetFillStyle(0);
    bSys[i_energy]->SetLineStyle(1);
    bSys[i_energy]->SetLineWidth(1);
    bSys[i_energy]->SetLineColor(plot_color);
    bSys[i_energy]->Draw("b Same");
  }
}
