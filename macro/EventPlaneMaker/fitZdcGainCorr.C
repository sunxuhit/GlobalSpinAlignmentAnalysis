#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void fitZdcGainCorr(int beamType = 0)
{
  string inputFile = Form("../../data/EventPlaneMaker/%s/file_GainCorr_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  std::string str_mEastWest[2] = {"East","West"};
  std::string str_mVertHori[2] = {"Vert","Hori"};

  TH2F *h_mGainInPut[2][2][8];
  TH1F *h_mGainCorr[2][2][8];
  double meanGain[2][2][8], rmsGain[2][2][8];
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	string histName = Form("h_mZdcGainCorr%s%sSlat%d",str_mEastWest[iEastWest].c_str(),str_mVertHori[iVertHori].c_str(),iSlat);
	h_mGainInPut[iEastWest][iVertHori][iSlat] = (TH2F*)file_InPut->Get(histName.c_str());

	histName = Form("h_mZdcGainCorr%s%sSlat%dADC",str_mEastWest[iEastWest].c_str(),str_mVertHori[iVertHori].c_str(),iSlat);
	h_mGainCorr[iEastWest][iVertHori][iSlat] = (TH1F*)h_mGainInPut[iEastWest][iVertHori][iSlat]->ProjectionY(histName.c_str());
	meanGain[iEastWest][iVertHori][iSlat] = h_mGainCorr[iEastWest][iVertHori][iSlat]->GetMean();
	rmsGain[iEastWest][iVertHori][iSlat] = h_mGainCorr[iEastWest][iVertHori][iSlat]->GetRMS();
	// cout << "iEastWest = " << iEastWest << ", iVertHori = " << iVertHori << ", iSlat = " << iSlat << ", meanGain = " << meanCain[iEastWest][iVertHori][iSlat] << ", rmsGain = " << rmsGain[iEastWest][iVertHori][iSlat] << endl;
      }
    }
  }

  // exponential fit to extract gain correction factor
  double norm[2][2][8], slope[2][2][8];
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	TF1 *f_exp =new TF1("f_exp","[0]*exp(-[1]*x)",0,1000);
	f_exp->SetParameter(0,0.1*h_mGainCorr[iEastWest][iVertHori][iSlat]->GetMaximum());
	f_exp->SetParameter(1,0.01);
	double expStart1st = meanGain[iEastWest][iVertHori][iSlat]+1.0*rmsGain[iEastWest][iVertHori][iSlat];
	double expStop1st  = meanGain[iEastWest][iVertHori][iSlat]+4.0*rmsGain[iEastWest][iVertHori][iSlat];
	f_exp->SetRange(expStart1st,expStop1st);
	h_mGainCorr[iEastWest][iVertHori][iSlat]->Fit(f_exp,"RN"); // first fit

	double expStart2nd = meanGain[iEastWest][iVertHori][iSlat]+3.0*rmsGain[iEastWest][iVertHori][iSlat];
	double expStop2nd  = meanGain[iEastWest][iVertHori][iSlat]+5.0*rmsGain[iEastWest][iVertHori][iSlat];
	f_exp->SetRange(expStart2nd,expStop2nd);
	h_mGainCorr[iEastWest][iVertHori][iSlat]->Fit(f_exp,"RQN"); // second fit
	double chi2 = f_exp->GetChisquare();
	double NDF = f_exp->GetNDF();
	// cout << "chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << endl;

	norm[iEastWest][iVertHori][iSlat] = f_exp->GetParameter(0);
	slope[iEastWest][iVertHori][iSlat] = f_exp->GetParameter(1);
      }
    }
  }

  // QA plots
  TCanvas *c_ZdcGainCorr = new TCanvas("c_ZdcGainCorr","c_ZdcGainCorr",10,10,1600,800);
  c_ZdcGainCorr->Divide(4,2);
  for(int iSlat = 0; iSlat < 8; ++iSlat)
  {
    c_ZdcGainCorr->cd(iSlat+1);
    c_ZdcGainCorr->cd(iSlat+1)->SetLeftMargin(0.15);
    c_ZdcGainCorr->cd(iSlat+1)->SetBottomMargin(0.15);
    c_ZdcGainCorr->cd(iSlat+1)->SetTicks(1,1);
    c_ZdcGainCorr->cd(iSlat+1)->SetGrid(0,0);
    c_ZdcGainCorr->cd(iSlat+1)->SetLogy(1);
  }
  std::string figName = Form("../../figures/EventPlaneMaker/%s/ZdcGainCorr_%s.pdf[",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_ZdcGainCorr->Print(figName.c_str());
  figName = Form("../../figures/EventPlaneMaker/%s/ZdcGainCorr_%s.pdf",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());

  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	c_ZdcGainCorr->cd(iSlat+1);
	h_mGainCorr[iEastWest][iVertHori][iSlat]->DrawCopy("hE");
	TF1 *f_exp =new TF1("f_exp","[0]*exp(-[1]*x)",0,1000);
	f_exp->FixParameter(0,norm[iEastWest][iVertHori][iSlat]);
	f_exp->FixParameter(1,slope[iEastWest][iVertHori][iSlat]);
	double expStart2nd = meanGain[iEastWest][iVertHori][iSlat]+3.0*rmsGain[iEastWest][iVertHori][iSlat];
	double expStop2nd  = meanGain[iEastWest][iVertHori][iSlat]+5.0*rmsGain[iEastWest][iVertHori][iSlat];
	f_exp->SetRange(expStart2nd,expStop2nd);
	f_exp->SetLineColor(2);
	f_exp->SetLineWidth(2);
	f_exp->SetLineStyle(1);
	f_exp->Draw("l same");
      }
      c_ZdcGainCorr->Update();
      c_ZdcGainCorr->Print(figName.c_str());
    }
  }
  figName = Form("../../figures/EventPlaneMaker/%s/ZdcGainCorr_%s.pdf]",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  c_ZdcGainCorr->Print(figName.c_str());

  // print Gain Correcction Factors
  cout << "const double mGainCorr[2][2][8] = {" << endl;
  cout << "    {" << endl;
  cout << "      { ";
  for(int iSlat =0; iSlat < 8; ++iSlat) 
  {
    cout << slope[0][0][0]/slope[0][0][iSlat];
    if(iSlat != 7) cout << ", ";
  }
  cout << "}," << endl; 
  cout << "      { ";
  for(int iSlat =0; iSlat <8; ++iSlat) 
  {
    cout << slope[0][1][0]/slope[0][1][iSlat];
    if(iSlat != 7) cout << ", ";
  }
  cout << "} " << endl;
  cout << "    }," << endl; 

  cout << "    {" << endl; 
  cout << "      {";
  for(int iSlat = 0; iSlat < 8; ++iSlat) 
  {
    cout << slope[1][0][0]/slope[1][0][iSlat];
    if(iSlat != 7) cout << ", ";
  }
  cout << "}," << endl; 
  cout << "      {";
  for(int iSlat =0; iSlat < 8; ++iSlat) 
  {
    cout << slope[1][1][0]/slope[1][1][iSlat];
    if(iSlat != 7) cout << ", ";
  }
  cout << "} " << endl;
  cout << "    }" << endl;
  cout << "}" << endl;

  // save to histograms
  TH1F *h_mGainCorrFactor[2][2][8];
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	string histName = Form("h_mZdcGainFactor%s%sSlat%d",str_mEastWest[iEastWest].c_str(),str_mVertHori[iVertHori].c_str(),iSlat);
	h_mGainCorrFactor[iEastWest][iVertHori][iSlat] = new TH1F(histName.c_str(),histName.c_str(),1,-0.5,0.5);
	h_mGainCorrFactor[iEastWest][iVertHori][iSlat]->SetBinContent(1,slope[iEastWest][iVertHori][0]/slope[iEastWest][iVertHori][iSlat]);
	h_mGainCorrFactor[iEastWest][iVertHori][iSlat]->SetBinError(1,0.1);
      }
    }
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/GainCorrPar/file_ZdcGainCorrFac_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	h_mGainCorrFactor[iEastWest][iVertHori][iSlat]->Write();
      }
    }
  }
  file_OutPut->Close();
  file_InPut->Close();
}
