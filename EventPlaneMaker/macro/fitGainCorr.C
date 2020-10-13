#include "../StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

void fitGainCorr(int energy = 0)
{
  string JobId = "9E5703EB6FAE0E39F93C889E8039552F";
  // string InPutFile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/GainCorr/merged_file/file_%s_GainCorrPar_%s.root",recoEP::mBeamEnergy[energy].c_str(),recoEP::mBeamEnergy[energy].c_str(),JobId.c_str());
  string InPutFile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/GainCorr/test/file_%s_GainCorr_%s.root",recoEP::mBeamEnergy[energy].c_str(),recoEP::mBeamEnergy[energy].c_str(),JobId.c_str());
  // string InPutFile = Form("/star/data01/pwg/sunxuhit/AuAu%s/SpinAlignment/GainCorrParameter/test/file_%s_GainCorr_%s.root",recoEP::mBeamEnergy[energy].c_str(),recoEP::mBeamEnergy[energy].c_str(),JobId.c_str());
  TFile *File_InPut = TFile::Open(InPutFile.c_str());
  if(!File_InPut->IsOpen()) cout << "InPutFile: " << InPutFile.c_str() << "is problematic" << endl;
  cout << "InPutFile sets to: " << InPutFile.c_str() << endl;

  TH2F *h_mGainInPut[2][2][8];
  TH1F *h_mGainCorr[2][2][8];
  float meanGain[2][2][8], rmsGain[2][2][8];
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	string HistName = Form("h_mZdcGainCorr%s%s_%d",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str(),i_slat);
	h_mGainInPut[i_eastwest][i_verthori][i_slat] = (TH2F*)File_InPut->Get(HistName.c_str());
	HistName = Form("h_mZdcGainCorr%s%s_%d_ADC",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str(),i_slat);
	h_mGainCorr[i_eastwest][i_verthori][i_slat] = (TH1F*)h_mGainInPut[i_eastwest][i_verthori][i_slat]->ProjectionY(HistName.c_str());
	meanGain[i_eastwest][i_verthori][i_slat] = h_mGainCorr[i_eastwest][i_verthori][i_slat]->GetMean();
	rmsGain[i_eastwest][i_verthori][i_slat] = h_mGainCorr[i_eastwest][i_verthori][i_slat]->GetRMS();
	// cout << "i_eastwest = " << i_eastwest << ", i_verthori = " << i_verthori << ", i_slat = " << i_slat << ", meanGain = " << meanCain[i_eastwest][i_verthori][i_slat] << ", rmsGain = " << rmsGain[i_eastwest][i_verthori][i_slat] << endl;
      }
    }
  }

  // exponential fit to extract gain correction factor
  float norm[2][2][8], slope[2][2][8];
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	TF1 *f_exp =new TF1("f_exp","[0]*exp(-[1]*x)",0,1000);
	f_exp->SetParameter(0,3e+6);
	f_exp->SetParameter(1,0.02);
	float exp_start_1st = meanGain[i_eastwest][i_verthori][i_slat]+1.0*rmsGain[i_eastwest][i_verthori][i_slat];
	float exp_stop_1st  = meanGain[i_eastwest][i_verthori][i_slat]+4.0*rmsGain[i_eastwest][i_verthori][i_slat];
	f_exp->SetRange(exp_start_1st,exp_stop_1st);
	h_mGainCorr[i_eastwest][i_verthori][i_slat]->Fit(f_exp,"RN"); // first fit

	float exp_start_2nd = meanGain[i_eastwest][i_verthori][i_slat]+3.0*rmsGain[i_eastwest][i_verthori][i_slat];
	float exp_stop_2nd  = meanGain[i_eastwest][i_verthori][i_slat]+5.0*rmsGain[i_eastwest][i_verthori][i_slat];
	f_exp->SetRange(exp_start_2nd,exp_stop_2nd);
	h_mGainCorr[i_eastwest][i_verthori][i_slat]->Fit(f_exp,"RQN"); // second fit
	float chi2 = f_exp->GetChisquare();
	float NDF = f_exp->GetNDF();
	// cout << "chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << endl;

	norm[i_eastwest][i_verthori][i_slat] = f_exp->GetParameter(0);
	slope[i_eastwest][i_verthori][i_slat] = f_exp->GetParameter(1);
      }
    }
  }

  // QA plots
  TCanvas *c_GainCorr[2][2];
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      string CanName = Form("c_mGainCorr%s%s",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str());
      c_GainCorr[i_eastwest][i_verthori] = new TCanvas(CanName.c_str(),CanName.c_str(),10,10,1600,800);
      c_GainCorr[i_eastwest][i_verthori]->Divide(4,2);
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	c_GainCorr[i_eastwest][i_verthori]->cd(i_slat+1);
	c_GainCorr[i_eastwest][i_verthori]->cd(i_slat+1)->SetLeftMargin(0.15);
	c_GainCorr[i_eastwest][i_verthori]->cd(i_slat+1)->SetBottomMargin(0.15);
	c_GainCorr[i_eastwest][i_verthori]->cd(i_slat+1)->SetTicks(1,1);
	c_GainCorr[i_eastwest][i_verthori]->cd(i_slat+1)->SetGrid(0,0);
	c_GainCorr[i_eastwest][i_verthori]->cd(i_slat+1)->SetLogy(1);
	h_mGainCorr[i_eastwest][i_verthori][i_slat]->Draw("hE");
	TF1 *f_exp =new TF1("f_exp","[0]*exp(-[1]*x)",0,1000);
	f_exp->FixParameter(0,norm[i_eastwest][i_verthori][i_slat]);
	f_exp->FixParameter(1,slope[i_eastwest][i_verthori][i_slat]);
	float exp_start_2nd = meanGain[i_eastwest][i_verthori][i_slat]+3.0*rmsGain[i_eastwest][i_verthori][i_slat];
	float exp_stop_2nd  = meanGain[i_eastwest][i_verthori][i_slat]+5.0*rmsGain[i_eastwest][i_verthori][i_slat];
	f_exp->SetRange(exp_start_2nd,exp_stop_2nd);
	f_exp->SetLineColor(2);
	f_exp->SetLineWidth(2);
	f_exp->SetLineStyle(1);
	f_exp->Draw("l same");
      }
      string FigureName = Form("./figures/c_mGainCorr%s%s.eps",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str());
      c_GainCorr[i_eastwest][i_verthori]->SaveAs(FigureName.c_str());
    }
  }

  // print Gain Correcction Factors
  cout << "const float mGainCorr[2][2][8] = {" << endl;
  cout << "    {" << endl;
  cout << "      { ";
  for(int i_slat =0; i_slat < 8; ++i_slat) 
  {
    cout << slope[0][0][0]/slope[0][0][i_slat];
    if(i_slat != 7) cout << ", ";
  }
  cout << "}," << endl; 
  cout << "      { ";
  for(int i_slat =0; i_slat <8; ++i_slat) 
  {
    cout << slope[0][1][0]/slope[0][1][i_slat];
    if(i_slat != 7) cout << ", ";
  }
  cout << "} " << endl;
  cout << "    }," << endl; 

  cout << "    {" << endl; 
  cout << "      {";
  for(int i_slat = 0; i_slat < 8; ++i_slat) 
  {
    cout << slope[1][0][0]/slope[1][0][i_slat];
    if(i_slat != 7) cout << ", ";
  }
  cout << "}," << endl; 
  cout << "      {";
  for(int i_slat =0; i_slat < 8; ++i_slat) 
  {
    cout << slope[1][1][0]/slope[1][1][i_slat];
    if(i_slat != 7) cout << ", ";
  }
  cout << "} " << endl;
  cout << "    }" << endl;
  cout << "}" << endl;

  // save to histograms
  TH1F *h_mGainCorrFactor[2][2][8];
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	string HistName = Form("h_mGainCorrFactor%s%s_%d",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str(),i_slat);
	h_mGainCorrFactor[i_eastwest][i_verthori][i_slat] = new TH1F(HistName.c_str(),HistName.c_str(),1,-0.5,0.5);
	h_mGainCorrFactor[i_eastwest][i_verthori][i_slat]->SetBinContent(1,slope[i_eastwest][i_verthori][0]/slope[i_eastwest][i_verthori][i_slat]);
	h_mGainCorrFactor[i_eastwest][i_verthori][i_slat]->SetBinError(1,0.1);
      }
    }
  }
  string OutPutFile = Form("../StRoot/StEventPlaneUtility/GainCorrPar/file_%s_GainCorrFac.root",recoEP::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	h_mGainCorrFactor[i_eastwest][i_verthori][i_slat]->Write();
      }
    }
  }
  File_OutPut->Close();
}
