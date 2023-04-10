#include <string>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include "../../Utility/include/draw.h"

int estStatError()
{
  gStyle->SetOptDate(0);

  float beamEnergy[3] = {11.5, 14.6, 19.6};
  float beamShiftPhi = -0.5;
  float beamShiftKstar = 0.5;
  float numEventBesI[3]  = {11.7,  12.6,  36.0};
  float numEventBesII[3] = {235.0, 324.0, 582.0};
  float statErrPhi[3]  = {0.018250128, 0.0, 0.0075676721};
  float statErrKstar[3]  = {0.0241712,  0.0207129,  0.0136966};
  float numEvent7 = 100.0;

  TCanvas *c_rho00 = new TCanvas("c_rho00","c_rho00",10,10,800,800);
  c_rho00->cd();
  c_rho00->cd()->SetLeftMargin(0.15);
  c_rho00->cd()->SetBottomMargin(0.15);
  c_rho00->cd()->SetTicks(1,1);
  c_rho00->cd()->SetGrid(0,0);
  // c_rho00->cd()->SetLogx();

  TH1F *h_frame = new TH1F("h_frame","h_frame",100,0,100.0);
  h_frame->SetTitle("");
  h_frame->SetStats(0);
  h_frame->GetXaxis()->SetRangeUser(0.0,30.0);
  h_frame->GetXaxis()->SetNdivisions(505,'N');
  h_frame->GetXaxis()->SetLabelSize(0.04);
  h_frame->GetXaxis()->SetTitle("#sqrt{s_{NN}}  (GeV)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->SetTitleOffset(1.1);
  h_frame->GetXaxis()->CenterTitle();

  h_frame->GetYaxis()->SetRangeUser(0.28,0.42);
  h_frame->GetYaxis()->SetNdivisions(505,'N');
  h_frame->GetYaxis()->SetTitle("#rho_{00} Stat. Error");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->SetTitleOffset(1.1);
  h_frame->GetYaxis()->SetLabelSize(0.04);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->DrawCopy("pE");
  PlotLine(0.0,30.0,1.0/3.0,1.0/3.0,1,2,2);

  TBox *b_phiStatBesI[3]; 
  TBox *b_phiStatBesII[3]; 
  for(int i_energy = 0; i_energy < 3; ++i_energy)
  {
    double x1BesI = beamEnergy[i_energy] + beamShiftPhi - 0.2;
    double y1BesI = 1.0/3.0 - statErrPhi[i_energy];
    double x2BesI = beamEnergy[i_energy] + beamShiftPhi + 0.2;
    double y2BesI = 1.0/3.0 + statErrPhi[i_energy];
    b_phiStatBesI[i_energy] = new TBox(x1BesI,y1BesI,x2BesI,y2BesI);
    b_phiStatBesI[i_energy]->SetFillColor(kMagenta+2);
    b_phiStatBesI[i_energy]->SetFillColorAlpha(kMagenta+2,0.65);
    b_phiStatBesI[i_energy]->SetFillStyle(3004);
    b_phiStatBesI[i_energy]->SetLineStyle(1);
    b_phiStatBesI[i_energy]->SetLineColor(kMagenta+2);
    b_phiStatBesI[i_energy]->SetLineWidth(1);
    if(i_energy != 1) b_phiStatBesI[i_energy]->Draw("l Same");

    double x1BesII = x1BesI+0.4;
    double y1BesII = 1.0/3.0 - statErrPhi[i_energy]*sqrt(numEventBesI[i_energy]/numEventBesII[i_energy]);
    double x2BesII = x2BesI+0.4;
    double y2BesII = 1.0/3.0 + statErrPhi[i_energy]*sqrt(numEventBesI[i_energy]/numEventBesII[i_energy]);
    b_phiStatBesII[i_energy] = new TBox(x1BesII,y1BesII,x2BesII,y2BesII);
    b_phiStatBesII[i_energy]->SetFillColor(kRed);
    b_phiStatBesII[i_energy]->SetFillColorAlpha(kRed,0.65);
    b_phiStatBesII[i_energy]->SetFillStyle(3004);
    b_phiStatBesII[i_energy]->SetLineStyle(1);
    b_phiStatBesII[i_energy]->SetLineColor(kRed);
    b_phiStatBesII[i_energy]->SetLineWidth(1);
    if(i_energy != 1) b_phiStatBesII[i_energy]->Draw("l Same");
  }
  TBox *b_phiStat7;
  {
    double x1BesII = 7.7 - 0.4;
    double y1BesII = 1.0/3.0 - statErrPhi[0]*sqrt(numEventBesI[0]/(0.8*numEvent7));
    double x2BesII = 7.7;
    double y2BesII = 1.0/3.0 + statErrPhi[0]*sqrt(numEventBesI[0]/(0.8*numEvent7));
    b_phiStat7 = new TBox(x1BesII,y1BesII,x2BesII,y2BesII);
    b_phiStat7->SetFillColor(kRed);
    b_phiStat7->SetFillColorAlpha(kRed,0.5);
    b_phiStat7->SetFillStyle(3002);
    b_phiStat7->SetLineStyle(2);
    b_phiStat7->SetLineColor(kRed);
    b_phiStat7->SetLineWidth(1);
    b_phiStat7->Draw("l Same");
  }

  TBox *b_kStarStatBesI[3]; 
  TBox *b_kStarStatBesII[3]; 
  for(int i_energy = 0; i_energy < 3; ++i_energy)
  {
    double x1BesI = beamEnergy[i_energy] + beamShiftKstar - 0.2;
    double y1BesI = 1.0/3.0 - statErrKstar[i_energy];
    double x2BesI = beamEnergy[i_energy] + beamShiftKstar + 0.2;
    double y2BesI = 1.0/3.0 + statErrKstar[i_energy];
    b_kStarStatBesI[i_energy] = new TBox(x1BesI,y1BesI,x2BesI,y2BesI);
    b_kStarStatBesI[i_energy]->SetFillColor(kGray+3);
    b_kStarStatBesI[i_energy]->SetFillColorAlpha(kGray+3,0.65);
    b_kStarStatBesI[i_energy]->SetFillStyle(3005);
    b_kStarStatBesI[i_energy]->SetLineStyle(1);
    b_kStarStatBesI[i_energy]->SetLineColor(kGray+3);
    b_kStarStatBesI[i_energy]->SetLineWidth(1);
    b_kStarStatBesI[i_energy]->Draw("l Same");

    double x1BesII = x1BesI+0.4;
    double y1BesII = 1.0/3.0 - statErrKstar[i_energy]*sqrt(numEventBesI[i_energy]/numEventBesII[i_energy]);
    double x2BesII = x2BesI+0.4;
    double y2BesII = 1.0/3.0 + statErrKstar[i_energy]*sqrt(numEventBesI[i_energy]/numEventBesII[i_energy]);
    b_kStarStatBesII[i_energy] = new TBox(x1BesII,y1BesII,x2BesII,y2BesII);
    b_kStarStatBesII[i_energy]->SetFillColor(kAzure);
    b_kStarStatBesII[i_energy]->SetFillColorAlpha(kAzure,0.65);
    b_kStarStatBesII[i_energy]->SetFillStyle(3005);
    b_kStarStatBesII[i_energy]->SetLineStyle(1);
    b_kStarStatBesII[i_energy]->SetLineColor(kAzure);
    b_kStarStatBesII[i_energy]->SetLineWidth(1);
    b_kStarStatBesII[i_energy]->Draw("l Same");
  }
  TBox *b_kStarStat7;
  {
    double x1BesII = 7.7;
    double y1BesII = 1.0/3.0 - statErrKstar[0]*sqrt(numEventBesI[0]/(0.8*numEvent7));
    double x2BesII = 7.7 + 0.4;
    double y2BesII = 1.0/3.0 + statErrKstar[0]*sqrt(numEventBesI[0]/(0.8*numEvent7));
    b_kStarStat7 = new TBox(x1BesII,y1BesII,x2BesII,y2BesII);
    b_kStarStat7->SetFillColor(kAzure);
    b_kStarStat7->SetFillColorAlpha(kAzure,0.5);
    b_kStarStat7->SetFillStyle(3002);
    b_kStarStat7->SetLineStyle(2);
    b_kStarStat7->SetLineColor(kAzure);
    b_kStarStat7->SetLineWidth(1);
    b_kStarStat7->Draw("l Same");
  }

  TLegend *leg = new TLegend(0.2,0.6,0.6,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(b_phiStatBesI[0],"#phi BES-I 20-60%","F");
  leg->AddEntry(b_phiStatBesII[0],"#phi BES-II 20-60% (Project)","F");
  leg->AddEntry(b_phiStat7,"#phi 7.7 GeV 20-60% (Estimate)","F");
  leg->AddEntry(b_kStarStatBesI[0],"K^{*0} BES-I 20-60%","F");
  leg->AddEntry(b_kStarStatBesII[0],"K^{*0} BES-II 20-60% (Project)","F");
  leg->AddEntry(b_kStarStat7,"K^{*0} 7.7 GeV 20-60% (Estimate)","F");
  leg->Draw("same");

  c_rho00->SaveAs("./c_rho00StatError.png");

  return 1;
}
