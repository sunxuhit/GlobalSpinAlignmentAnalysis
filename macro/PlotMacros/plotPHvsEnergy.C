const int ndatH = 2;
const int nSH   = 5; 
TGraph *gHydro[ndatH];
void LoadHydro();

TGraphErrors *gPolAMPT[2];
void LoadAMPT();

TGraphErrors *gPolBES[2];
TGraphAsymmErrors *gSEPolBES[2];
void LoadBES();

TGraphErrors *gPolRun4[2];
void LoadRun4();

TGraphErrors *gPolRun11[2];
TGraphAsymmErrors *gSEPolRun11[2];

TGraphErrors *gPolRunAll[2];
TGraphAsymmErrors *gSEPolRunAll[2];
void LoadRunAll();

TGraphErrors *gpolQM2019[3];
TGraphAsymmErrors *gSEpolQM2019[3];
void LoadQM2019();

TGraphErrors *gpol7p2GeV;
TGraphAsymmErrors *gSEpol7p2GeV;
void Load7p2GeV();

TGraphErrors *gpolHADES[2];
TGraphAsymmErrors *gSEpolHADES[2];
void LoadHADES();

TGraphErrors *gpolALICE[2];
TGraphAsymmErrors *gSEpolALICE[2];
void LoadALICE();

TGraphErrors *gpol3GeV;
TGraphAsymmErrors *gSEpol3GeV;
void Load3GeV();

TGraphErrors *gpolRuRu200GeV[2];
TGraphAsymmErrors *gSEpolRuRu200GeV[2];
void LoadRuRu200GeV();

TGraphErrors *gpolZrZr200GeV[2];
TGraphAsymmErrors *gSEpolZrZr200GeV[2];
void LoadZrZr200GeV();



void SetAxis( TH1F *ax, float xoff=1.5, float yoff=1.4 );
void SetLegend( TLegend *l );



void plotPHvsEnergy(){

  LoadHydro();
  LoadBES();
  LoadRun4();
  LoadRunAll();
  //Load7p2GeV();
  LoadHADES();
  LoadALICE();
  Load3GeV();
  LoadAMPT();
  //LoadQM2019();
  //LoadRuRu200GeV();
  //LoadZrZr200GeV();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadTopMargin(0.05);

  TH1F *ax;

  TF1 *zero = new TF1("zero","[0]",1,8000);
  zero->SetLineStyle(2);
  zero->SetLineWidth(1);
  zero->SetLineColor(1);

  TLatex lt;
  lt.SetNDC();
  lt.SetTextFont(43);
  lt.SetTextSize(12);

  TLatex lt2;
  lt2.SetNDC();
  lt2.SetTextFont(43);
  lt2.SetTextSize(14);
  lt2.SetTextColor(1);


  TCanvas *c1 = new TCanvas("c1","c1",500,400);
  c1->cd();
  gPad->SetBottomMargin(0.11);
  ax = gPad->DrawFrame( 1.0, -2.5, 6000, 8.5 );  // to show all data
  gPad->SetLogx();
  SetAxis(ax,1.0,1.3);
  ax->SetXTitle("#sqrt{s_{NN}} [GeV] ");
  ax->SetYTitle("P_{H} [%] ");
  zero->Draw("same");

  // Hydro
  gHydro[0]->Draw("L");
  gHydro[1]->Draw("L");

  // AMPT
  gPolAMPT[0]->Draw("L");
  gPolAMPT[1]->Draw("L");


  for( int ich=0; ich<2; ich++ ){
    gSEPolBES[ich]->Draw("2");
    gPolBES[ich]->Draw("PZ");
    gPolRun4[ich]->Draw("PZ");
  }

  //for(int ich=0;ich<3;ich++){
  //		gpolQM2019[ich]         -> Draw("PZ");
  //		gSEpolQM2019[ich]       -> Draw("2");
  //}


  for(int ich=0;ich<2;ich++){
    gPolRunAll[ich]        -> Draw("PZ"); 
    gSEPolRunAll[ich]      -> Draw("2"); 
    gpolALICE[ich]         -> Draw("PZ");
    gSEpolALICE[ich]       -> Draw("2");
    //gpolRuRu200GeV[ich]    -> Draw("PZ");
    //gSEpolRuRu200GeV[ich]  -> Draw("2");
    //gpolZrZr200GeV[ich]    -> Draw("PZ");
    //gSEpolZrZr200GeV[ich]  -> Draw("2");
  }
  //gpol7p2GeV       -> Draw("PZ");
  //gSEpol7p2GeV     -> Draw("2");
  gpolHADES[0]->Draw("PZ");
  gpolHADES[1]->Draw("PZ");
  gSEpolHADES[0]->Draw("2");
  gSEpolHADES[1]->Draw("2");

  gpol3GeV   -> Draw("PZ");
  gSEpol3GeV -> Draw("2");


  lt.DrawLatex( 0.605103,0.9034505, "STAR Au+Au 20%-50%");
  lt2.DrawLatex( 0.700946,0.185341, "#alpha_{#Lambda} = -#alpha_{#bar{#Lambda}} =0.732");

  TLegend *lSTAR = new TLegend(0.5728346,0.726,0.9507874,0.894,NULL,"brNDC");
  SetLegend(lSTAR);
  lSTAR->SetNColumns(3);
  lSTAR->AddEntry( gPolBES[0], "Nature548.62 (2017)  ", "0" );
  lSTAR->AddEntry( gPolBES[0], "#Lambda      ", "P" );
  lSTAR->AddEntry( gPolBES[1], "#bar{#Lambda} ", "P" );
  lSTAR->AddEntry( gPolRun4[0], "PRC76.024915 (2007)  ", "0" );
  lSTAR->AddEntry( gPolRun4[0], "#Lambda      ", "P" );
  lSTAR->AddEntry( gPolRun4[1], "#bar{#Lambda} ", "P" );
  lSTAR->AddEntry( gPolRunAll[0], "PRC98.014910 (2018)  ", "0" );
  lSTAR->AddEntry( gPolRunAll[0], "#Lambda      ", "P" );
  lSTAR->AddEntry( gPolRunAll[1], "#bar{#Lambda} ", "P" );
  lSTAR->AddEntry(gpol3GeV,"PRC104.L061901 (2021)","0");
  lSTAR->AddEntry(gpol3GeV,"#Lambda","P");
  lSTAR->Draw();

  //TLegend *lQM2019 = new TLegend(0.6043307,0.6433333,0.8661417,0.7213333,NULL,"brNDC");
  //SetLegend(lQM2019);
  //lQM2019->SetNColumns(3);
  //lQM2019->SetHeader("STAR prelim.");
  //lQM2019->AddEntry(gpolQM2019[0],"#Lambda+#bar{#Lambda}","P");
  //lQM2019->AddEntry(gpolQM2019[1],"#Lambda","P");
  //lQM2019->AddEntry(gpolQM2019[2],"#bar{#Lambda}","P");
  //lQM2019->Draw();

  //TLegend *lRuZr200GeV = new TLegend(0.5983307,0.5733333,0.8661417,0.643333,NULL,"brNDC");
  //SetLegend(lRuZr200GeV);
  //lRuZr200GeV->SetNColumns(2);
  //lRuZr200GeV->AddEntry(gpolRuRu200GeV[0],"#Lambda      ","P");
  //lRuZr200GeV->AddEntry(gpolRuRu200GeV[1],"#bar{#Lambda}        Ru+Ru 20-50%","P");
  //lRuZr200GeV->AddEntry(gpolZrZr200GeV[0],"#Lambda      ","P");
  //lRuZr200GeV->AddEntry(gpolZrZr200GeV[1],"#bar{#Lambda}        Zr+Zr    20-50%","P");
  //lRuZr200GeV->Draw();

  TLegend *lALICE = new TLegend(0.6003622,0.6413333,0.9606299,0.7266667,NULL,"brNDC");
  SetLegend(lALICE);
  lALICE->SetNColumns(2);
  lALICE->SetHeader("ALICE PRC101.044611 (2020)");
  lALICE->AddEntry(gpolALICE[0],"#Lambda","P");
  lALICE->AddEntry(gpolALICE[1],"#bar{#Lambda}    Pb+Pb 15-50%","P");
  lALICE->Draw();

  TLegend *lHADES = new TLegend(0.6003307,0.5253333,0.8622047,0.630333,NULL,"brNDC");
  SetLegend(lHADES);
  lHADES->SetHeader("HADES PLB835.137506 (2022)");
  lHADES->AddEntry(gpolHADES[0],"#Lambda    Au+Au 20-40%","P");
  lHADES->AddEntry(gpolHADES[1],"#Lambda    Ag+Ag 20-40%","P");
  lHADES->Draw();


  TLegend *lhy = new TLegend(0.1909449,0.136,0.5846457,0.2053333,NULL,"brNDC");
  SetLegend(lhy);
  lhy->SetNColumns(2);
  lhy->SetHeader("UrQMD+vHLLE, #Lambda");
  lhy->AddEntry( gHydro[0], "primary", "L" );
  lhy->AddEntry( gHydro[1], "primary+feed-down", "L" );
  lhy->Draw();

  TLegend *lampt = new TLegend(0.1929134,0.216,0.5807087,0.2693333,NULL,"brNDC");
  SetLegend(lampt);
  lampt->SetNColumns(2);
  lampt->SetHeader("AMPT, #Lambda");
  lampt->AddEntry( gPolAMPT[0], "primary", "L" );
  lampt->AddEntry( gPolAMPT[1], "primary+feed-down", "L" );
  lampt->Draw();

  // future measurement
  TBox *bHighMuB= new TBox(2.0, 0.0, 7.5, 5.0);
  bHighMuB->SetFillColor(kYellow);
  bHighMuB->SetFillStyle(3001);
  bHighMuB->Draw("same");

  // NICA @ JINR: 3 - 11 GeV
  TLatex *lNica;
  lNica = new TLatex(4.5, 3.55,"NICA");
  lNica->SetTextFont(42);
  lNica->SetTextSize(0.023);
  lNica->SetTextColor(1);
  lNica->SetLineWidth(2);
  lNica->Draw();
  TArrow *aNica;
  aNica = new TArrow(3.0, 3.5, 11.0, 3.5, 0.02,"<>");
  aNica->SetFillColor(1);
  aNica->SetFillStyle(1001);
  aNica->SetLineColor(17);
  aNica->SetLineWidth(1);
  aNica->SetLineStyle(1);
  aNica->SetAngle(30);
  aNica->Draw();
  
  //FXT : 3.0 - 7.7
  TLatex *lFxt;
  lFxt = new TLatex(2.75, 2.95,"STAR-FXT");
  lFxt->SetTextFont(42);
  lFxt->SetTextSize(0.023);
  lFxt->SetTextColor(kRed-2);
  lFxt->SetLineWidth(2);
  lFxt->Draw();
  TArrow *aFxt;
  aFxt = new TArrow(3.0, 2.9, 7.7, 2.9, 0.02,"<>");
  aFxt->SetFillColor(kRed-7);
  aFxt->SetFillStyle(1002);
  aFxt->SetLineColor(kRed-7);
  aFxt->SetLineWidth(1);
  aFxt->SetAngle(30);
  aFxt->Draw();
  
  //BES-II: 7.7 - 27.0
  TLatex *lBESII;
  lBESII = new TLatex(8.75, 0.05,"STAR-BESII");
  lBESII->SetTextFont(42);
  lBESII->SetTextSize(0.023);
  lBESII->SetTextColor(kRed-2);
  lBESII->SetLineWidth(2);
  lBESII->Draw();
  TArrow *aBESII;
  aBESII = new TArrow(7.7, 0.35, 27.0, 0.35, 0.02,"<>");
  aBESII->SetFillColor(kRed-7);
  aBESII->SetFillStyle(1002);
  aBESII->SetLineColor(kRed-7);
  aBESII->SetLineWidth(1);
  aBESII->SetAngle(30);
  aBESII->Draw();
  
  // FAIR @ GSI: 2 - 4.9 GeV
  TLatex *lFair;
  lFair = new TLatex(2.40, 2.35,"FAIR");
  lFair->SetTextFont(42);
  lFair->SetTextSize(0.02283105);
  lFair->SetLineWidth(1);
  lFair->Draw();
  TArrow *aFair;
  aFair = new TArrow(2.0, 2.3, 4.9, 2.3, 0.02,"<>");
  aFair->SetFillColor(1);
  aFair->SetFillStyle(1001);
  aFair->SetLineColor(17);
  aFair->SetLineWidth(1);
  aFair->SetAngle(30);
  aFair->Draw();
  
  // CEE @ HIAF: 2 - 4 GeV
  TLatex *lHiaf;
  lHiaf = new TLatex(1.7, 1.75,"CEE@HIAF");
  lHiaf->SetTextFont(42);
  lHiaf->SetTextSize(0.022);
  lHiaf->SetTextColor(kRed-2);
  lHiaf->SetLineWidth(1);
  lHiaf->Draw();
  TArrow *aHiaf;
  aHiaf = new TArrow(2.0, 1.7, 4.0, 1.7, 0.02,"<>");
  aHiaf->SetFillColor(1);
  aHiaf->SetFillStyle(1001);
  aHiaf->SetLineColor(kRed-7);
  aHiaf->SetLineWidth(1);
  aHiaf->SetAngle(30);
  aHiaf->Draw();
  
  // CEE @ HIRFL: 2.11 GeV
  TLatex *lHirfl;
  lHirfl = new TLatex(2.5, 1.15,"CEE@HIRFL");
  lHirfl->SetTextFont(42);
  lHirfl->SetTextSize(0.022);
  lHirfl->SetTextColor(kRed-2);
  lHirfl->SetLineWidth(1);
  lHirfl->Draw();
  TArrow *aHirflL = new TArrow(1.65, 1.0, 2.11, 1.0, 0.02,">");
  aHirflL->SetFillColor(1);
  aHirflL->SetFillStyle(1001);
  aHirflL->SetLineColor(kRed-7);
  aHirflL->SetLineWidth(1);
  aHirflL->SetAngle(30);
  aHirflL->Draw();
  TArrow *aHirflR = new TArrow(2.32, 1.0, 3.0, 1.0, 0.02,"<");
  aHirflR->SetFillColor(1);
  aHirflR->SetFillStyle(1001);
  aHirflR->SetLineColor(kRed-7);
  aHirflR->SetLineWidth(1);
  aHirflR->SetAngle(30);
  aHirflR->Draw();

  c1->SaveAs("../../figures/AnalysisNote/fig_PhFutureMeasurement.eps");
}

void LoadHydro(){

  // from plot digitizer
  //=============================
  ifstream inH; 
  //	inH.open("/Users/niida/Documents/WSU/HydroBesGPol/gpol_hydro.txt");
  //	if( !inH ){ cout <<"Could not open gpol_hydro.txt!!"<< endl; return; }

  // from Yuri
  //=============================
  double PHp[] = { 0.01749583,  0.00988011,  0.00584896,  0.00403138,  0.00170648 }; // primary
  double PHf[] = { 0.01470586,  0.00835941,  0.00493472,  0.00340938,  0.00144373 }; // primary + feed-down

  // Load data1 //
  //=============================
  string tmp;
  for( int i=0; i<ndatH*2; i++ ){ 
    inH >> tmp;
    cout << tmp <<" ";
  }
  cout << endl;

  for( int idat=0; idat<ndatH; idat++ ){
    gHydro[idat] = new TGraph(nSH);
    gHydro[idat]->SetLineColor(kBlue);
    gHydro[idat]->SetLineColor(kRed);
    gHydro[idat]->SetLineStyle(1+idat*6);
    gHydro[idat]->SetLineWidth(2);
  }

  double fx, fval;
  double rSH[5] = { 7.7, 19.6, 39, 62.4, 200 };
  for( int ix=0; ix<nSH; ix++ ){
    for( int idat=0; idat<ndatH; idat++ ){
      //inH >> fx >> fval; 
      //gHydro[idat]->SetPoint( ix, rSH[ix], fval*100. );

      if( idat==0 ) gHydro[idat]->SetPoint( ix, rSH[ix], PHp[ix]*100. );
      else          gHydro[idat]->SetPoint( ix, rSH[ix], PHf[ix]*100. );
    }	}
  //=============================
  inH.close();

}

void LoadBES(){

  Double_t roots[6]   = {7.682,    11.454,   14.5546,  19.564,   26.994,   38.996};
  Double_t xerr[6] = {0};
  Double_t xerr2[6] = {0};
  Double_t ScaleF=0.642/0.732;

  Double_t Plam[6] = {2.03922, 1.34356, 1.32071, 0.950382, 1.04722, 0.506326};
  Double_t dPlam[6] = {0.628254, 0.395796, 0.481546, 0.305166, 0.28208, 0.423624};
  Double_t Palam[6] = {8.6691, 1.80167, 2.27599, 1.51528, 1.24525, 0.938276};
  //Double_t Palam[6] = {10.6691, 1.80167, 2.27599, 1.51528, 1.24525, 0.938276};
  Double_t dPalam[6] = {3.56897, 1.26119, 1.20993, 0.610072, 0.471465, 0.61503};

  Double_t PalamCr[6];
  Double_t dPalamCr[6];

  Double_t PlamSystUpper[8]  = {0.00,     0.00,     0.00,     0.00,     0.00,     0.00,     0.00,    0.00};
  Double_t PlamSystLower[8]  = {0.20,     0.20,     0.30,     0.20,     0.20,     0.20,     0.00,    0.00};
  Double_t PalamSystUpper[8] = {0.00,     0.00,     0.40,     0.00,     0.00,     0.00,     0.00,    0.00};
  Double_t PalamSystLower[8] = {1.00,     0.15,     0.15,     0.15,     0.15,     0.15,     0.00,    0.00};

  Double_t roots2[6];
  Double_t SEPlam[2][6];
  Double_t SEPalam[2][6];
  for( int is=0; is<6; is++ ){ 
    roots2[is] = roots[is]*1.05;
    xerr2 [is] = roots[is]*0.08;
    SEPlam[0][is] =  PlamSystLower[is];  //Plam[is]*PlamSystLower[is];
    SEPlam[1][is] =  PlamSystUpper[is];  //Plam[is]*PlamSystUpper[is];
    SEPalam[0][is] = PalamSystLower[is];//Palam[is]*PalamSystLower[is];
    SEPalam[1][is] = PalamSystUpper[is];//Palam[is]*PalamSystUpper[is];

    PalamCr[is] = Palam[is] * 0.642 / 0.71;
    dPalamCr[is] = dPalam[is] * 0.642 / 0.71;
  }

  for(int i=0;i<6;i++){
    Plam[i]*=ScaleF;
    Palam[i]*=ScaleF;
    dPlam[i]*=ScaleF;
    dPalam[i]*=ScaleF;
    SEPlam[0][i]*=ScaleF;
    SEPlam[1][i]*=ScaleF;
    SEPalam[0][i]*=ScaleF;
    SEPalam[1][i]*=ScaleF;
    cout << "Lambda      x = " << roots[i] << "   PH = " << Plam[i] << " +/- " << dPlam[i] << " + " << SEPlam[1][i] << " - " << SEPlam[0][i] << endl;
    //cout << "LambdaBar   x = " << roots[i] << "   PH = " << Palam[i] << " +/- " << dPalam[i] << " + " << SEPalam[1][i] << " - " << SEPalam[0][i] << endl;
  }

  gPolBES[0] = new TGraphErrors( 6, roots, Plam, xerr, dPlam );
  gPolBES[1] = new TGraphErrors( 6, roots2, Palam, xerr, dPalam );
  gPolBES[0]->SetMarkerStyle(20);
  gPolBES[1]->SetMarkerStyle(24);
  gPolBES[0]->SetMarkerSize(0.8); 
  gPolBES[1]->SetMarkerSize(0.8);

  gSEPolBES[0] = new TGraphAsymmErrors( 6, roots,  Plam,  xerr2, xerr2, SEPlam[0],  SEPlam[1] );
  gSEPolBES[1] = new TGraphAsymmErrors( 6, roots2, Palam, xerr2, xerr2, SEPalam[0], SEPalam[1] );
  gSEPolBES[0]->SetFillStyle(0);
  gSEPolBES[1]->SetFillStyle(0);
  gSEPolBES[0]->SetLineColor(1);
  gSEPolBES[1]->SetLineColor(1);


}

void LoadRun4(){

  Double_t roots[2] = { 62.4, 200 };
  Double_t xerr[2] = {0};

  Double_t Plam [2] = { 0.0133393419, 0.0012491361 };
  Double_t EPlam[2] = { 0.0116727212, 0.0101871395 };
  Double_t Palam [2] = { 0.0171184549, -0.0077602848 };
  Double_t EPalam[2] = { 0.015918979, 0.0115104159 };

  Double_t ScaleF=0.642/0.732;
  Double_t roots2[2];
  for( int is=0; is<2; is++ ){ 
    roots2[is] = roots[is]*1.05;
    Plam[is] *= 100.;
    EPlam[is] *= 100.;
    Palam[is] *= 100.;
    EPalam[is] *= 100.;
    Plam[is] *= ScaleF;
    Palam[is] *= ScaleF;
    EPlam[is] *= ScaleF;
    EPalam[is] *= ScaleF;
    cout << "Lambda      x = " << roots[is] << "    PH = " << Plam[is] << " +/- " << EPlam[is] << endl;
    cout << "LambdaBar   x = " << roots[is] << "    PH = " << Palam[is] << " +/- " << EPalam[is] << endl;
  }
  roots [1] -= 20;
  roots2[1] -= 20;

  gPolRun4[0] = new TGraphErrors( 2, roots, Plam, xerr, EPlam );
  gPolRun4[1] = new TGraphErrors( 2, roots2, Palam, xerr, EPalam );
  gPolRun4[0]->SetMarkerStyle(20);
  gPolRun4[1]->SetMarkerStyle(24);
  gPolRun4[0]->SetMarkerColor(kViolet+1);
  gPolRun4[1]->SetMarkerColor(kViolet+1);
  gPolRun4[0]->SetMarkerSize(0.8);
  gPolRun4[1]->SetMarkerSize(0.8);
}

void LoadAMPT(){

  ifstream in;
  //in.open("polarization_hui.txt");
  if( !in ){
    cout <<"could not open AMPT data file!"<< endl;
    return;
  }

  double PHp[] = {2.98, 2.32, 2.14, 1.76, 1.35, 1.05, 0.736, 0.279};
  double PHf[] = {2.57, 1.91, 1.79, 1.4, 1.13, 0.89, 0.632, 0.239};
  double px[8] = {8.32, 11.9, 12.9, 17.6, 24.6, 36.2, 59.3, 200};
  double fx[8] = {7.76, 12.1, 12.9, 18.4, 24.3, 35.7, 58.1, 193};

  gPolAMPT[0] = new TGraphErrors(8);
  gPolAMPT[1] = new TGraphErrors(8);
  gPolAMPT[0]->SetLineColor(kMagenta-9);
  gPolAMPT[1]->SetLineColor(kMagenta-9);
  gPolAMPT[1]->SetLineStyle(7);
  gPolAMPT[0]->SetLineWidth(5);
  gPolAMPT[1]->SetLineWidth(5);

  //	string tmp;
  //	in >> tmp >> tmp >> tmp >> tmp;
  //
  //	float data[3];
  //	cout << endl <<"AMPT data"<< endl;
  for( int is=0; is<8; is++ ){
    //		in >> data[0] >> data[1] >> data[2];	
    //		cout << data[0] <<" "<< data[1] <<" "<< data[2] << endl;
    //		
    //		gPolAMPT[0]->SetPoint( is, data[0], data[1] );
    //		gPolAMPT[1]->SetPoint( is, data[0], data[2] );
    //
    gPolAMPT[0]->SetPoint( is, px[is], PHp[is] );
    gPolAMPT[1]->SetPoint( is, fx[is], PHf[is] );
  }
  //	cout << endl;
}

void LoadRunAll(){

  for( int ich=0; ich<2; ich++ ){

    gPolRunAll  [ich] = new TGraphErrors(1);
    gSEPolRunAll[ich] = new TGraphAsymmErrors(1);
    gPolRunAll  [ich]->SetMarkerStyle(20+ich*4);
    gPolRunAll  [ich]->SetMarkerSize(0.8);
    gPolRunAll  [ich]->SetMarkerColor(kGreen+2);
    gPolRunAll  [ich]->SetLineColor(kGreen+2);

    gSEPolRunAll[ich]->SetFillStyle(0);
    gSEPolRunAll[ich]->SetFillColor(kGreen+2);
    gSEPolRunAll[ich]->SetLineColor(kGreen+2);

    if( ich==1 ){
      //gPolRunAll  [ich]->SetMarkerColor(kGreen+2);
      //gPolRunAll  [ich]->SetLineColor(kGreen+2);
      //gSEPolRunAll[ich]->SetFillColor(kGreen+2);
      //gSEPolRunAll[ich]->SetLineColor(kGreen+2);
    }
  }

  // 20-50%
  //runAll ich=0 pol[%]=0.277435 +- 0.0402006 + 0.0390459 - 0.0487696
  //runAll ich=1 pol[%]=0.239503 +- 0.0448805 + 0.0606676 - 0.0453513
  float fPH[2];
  float ePH[2];
  fPH[0] = 0.277435;
  fPH[1] = 0.239503;
  ePH[0] = 0.0402006;
  ePH[1] = 0.0448805;

  float seL[2], seH[2];
  seH[0] = 0.0390459;
  seH[1] = 0.0606676;
  seL[0] = 0.0487696;
  seL[1] = 0.0453513;

  double fx = 208;
  Double_t ScaleF=0.642/0.732;
  for( int ich=0; ich<2; ich++ ){
    fPH[ich]*=ScaleF;
    ePH[ich]*=ScaleF;
    seL[ich]*=ScaleF;
    seH[ich]*=ScaleF;
    gPolRunAll  [ich]->SetPoint     ( 0, fx+ich*20, fPH[ich] );
    gPolRunAll  [ich]->SetPointError( 0, 0,         ePH[ich] );
    gSEPolRunAll[ich]->SetPoint     ( 0, fx+ich*20, fPH[ich] );
    gSEPolRunAll[ich]->SetPointError( 0, 15, 20, seL[ich], seH[ich] );
    cout << "ich = " << ich << "  x = " << fx << "    PH = " << fPH[ich] << " +/- " << ePH[ich] << endl;
  }

}

void LoadQM2019(){
  int MarkerS[3] = {22,20,24};
  double size[3] = {0.8,0.8,0.8};
  int color=kRed;

  for( int is=0; is<3; is++){

    gpolQM2019[is] = new TGraphErrors(1);
    gSEpolQM2019[is] = new TGraphAsymmErrors(1);
    gpolQM2019 [is]->SetMarkerStyle(MarkerS[is]);
    gpolQM2019 [is]->SetMarkerSize(size[is]);
    gpolQM2019 [is]->SetFillStyle(1);
    gpolQM2019 [is] ->SetMarkerColor(color);
    gpolQM2019 [is] ->SetLineColor(color);
    gpolQM2019 [is] ->SetFillColor(color);
    gpolQM2019 [is] ->SetLineWidth(1);

    gSEpolQM2019 [is]->SetFillStyle(0);
    gSEpolQM2019 [is]->SetFillColor(color);
    gSEpolQM2019 [is]->SetLineColor(color);
    gSEpolQM2019 [is]->SetLineWidth(1);
  }


  float fPH[3];//0-> 27GeV  1-> 54GeV Lambda 2-> 54GeV antiLambda;
  float ePH[3];
  fPH[0] = 0.00598;
  ePH[0] = 0.00055;
  fPH[1] = 0.0055814;
  ePH[1] = 0.00231977;
  fPH[2] = 0.00713799;
  ePH[2] = 0.00317908;

  double exL[3] = {0.065,0.120938,0.22537486};
  double exH[3] = {0.063,0.120934,0.216585};

  Double_t ScaleF=0.642/0.732;
  for(int ich=0;ich<3;ich++){
    fPH[ich] *= ScaleF;
    ePH[ich] *= ScaleF;
    exL[ich] *= ScaleF;
    exH[ich] *= ScaleF;
  }

  double fx[3] = {27,51.4,57.4};
  //double width[3] = {2.0,2.5,2.5};
  double width[3] = {0, 0, 0};
  for(int ich=0; ich<3; ich++){
    width[ich] = fx[ich]*0.08;
    gpolQM2019[ich]->SetPoint     ( 0, fx[ich],fPH[ich]*100 );
    gpolQM2019[ich]->SetPointError( 0, 0, ePH[ich]*100 );
    gSEpolQM2019[ich]->SetPoint     ( 0, fx[ich],fPH[ich]*100 );
    gSEpolQM2019[ich]->SetPointError( 0, width[ich],width[ich], exL[ich],exH[ich]);
  }


}

void Load7p2GeV(){//7.2 GeV
  int MarkerS = 8;
  float Size = 0.9;
  int Color = kRed;

  gpol7p2GeV   = new TGraphErrors(1);
  gSEpol7p2GeV = new TGraphAsymmErrors(1);
  gpol7p2GeV -> SetMarkerStyle(MarkerS);
  gpol7p2GeV -> SetMarkerSize(Size);
  gpol7p2GeV -> SetFillStyle(1);
  gpol7p2GeV -> SetMarkerColor(Color);
  gpol7p2GeV -> SetLineColor(Color);
  gpol7p2GeV -> SetFillColor(Color);
  gpol7p2GeV -> SetLineWidth(1);
  gSEpol7p2GeV -> SetFillStyle(0);
  gSEpol7p2GeV -> SetLineColor(Color);

  //Lambda
  double fPHLambda = 2.8919;//7.2GeV
  double ePHLambda = 0.336956;


  double fxLambda = 7.2;

  double eBLambda=0.432486;//7.2GeV
  double eTLambda=1.00774;

  double width = 0.5;//7.2GeV

  gpol7p2GeV -> SetPoint     (0, fxLambda, fPHLambda );
  gpol7p2GeV -> SetPointError(0,       0, ePHLambda );
  gSEpol7p2GeV -> SetPoint     (0, fxLambda, fPHLambda );
  gSEpol7p2GeV -> SetPointError(0, width, width, eBLambda, eTLambda );
}


void Load3GeV(){
  gpol3GeV = new TGraphErrors(1);
  gSEpol3GeV = new TGraphAsymmErrors(1);
  gpol3GeV -> SetMarkerStyle(20);
  gpol3GeV -> SetFillStyle(1);
  gpol3GeV -> SetMarkerSize(0.8);
  gpol3GeV -> SetMarkerColor(kAzure+8);
  gpol3GeV ->   SetLineColor(kAzure+8);
  gpol3GeV ->   SetFillColor(kAzure+8);
  gpol3GeV -> SetLineWidth(2);

  gSEpol3GeV -> SetFillStyle(0);
  gSEpol3GeV -> SetFillColor(kAzure+8);
  gSEpol3GeV -> SetLineColor(kAzure+8);
  gSEpol3GeV -> SetLineWidth(1);

  float fx = 3.0;
  float fPH = 4.91;
  float ePH = 0.81;
  float sePH = 0.15;
  float width = fx*0.07;
  gpol3GeV   -> SetPoint     (0,fx,fPH);
  gpol3GeV   -> SetPointError(0,0,ePH);
  gSEpol3GeV -> SetPoint     (0,fx,fPH);
  gSEpol3GeV -> SetPointError(0,width,width*1.1,sePH,sePH);
  cout << "x = " << fx << "    PH = " << fPH << " +/- " << ePH << endl;
}



void LoadHADES(){
  int Mstyle[2]={33,29};
  float Msize[2]={1.3,0.9};
  for(int i=0;i<2;i++){
    gpolHADES[i] = new TGraphErrors(1);
    gSEpolHADES[i] = new TGraphAsymmErrors(1);
    gpolHADES[i] -> SetMarkerStyle(Mstyle[i]);
    gpolHADES[i] -> SetFillStyle(1);
    gpolHADES[i] -> SetMarkerSize(1.3);
    gpolHADES[i] -> SetMarkerColor(kOrange+1);
    gpolHADES[i] ->   SetLineColor(kOrange+1);
    gpolHADES[i] ->   SetFillColor(kOrange+1);
    gpolHADES[i] -> SetLineWidth(2);

    gSEpolHADES[i] -> SetFillStyle(0);
    gSEpolHADES[i] -> SetFillColor(kOrange+1);
    gSEpolHADES[i] -> SetLineColor(kOrange+1);
    gSEpolHADES[i] -> SetLineWidth(1);
  }


  //float fx[2] = {2.4,2.7};
  //float fPH[2] = {4.77,3.19};
  //float ePH[2] = {0.97,0.30};
  //float sePH[2] = {1.57,0.31};
  //float width[2] = {0.13,0.13};
  float fx[2] = {2.35,2.65}; //AuAu, AgAg
  float fPH[2] = {6.8, 6.2};
  float ePH[2] = {1.3, 0.4};
  float sePH[2] = {2.1, 0.6};
  float width[2] = {0.13,0.13};
  for(int i=0;i<2;i++){
    gpolHADES[i]   -> SetPoint     (0,fx[i],fPH[i]);
    gpolHADES[i]   -> SetPointError(0,0,ePH[i]);
    gSEpolHADES[i] -> SetPoint     (0,fx[i],fPH[i]);
    gSEpolHADES[i] -> SetPointError(0,width[i],width[i],sePH[i],sePH[i]);
    cout << "x = " << fx[i] << "    PH = " << fPH[i] << " +/- " << ePH[i] << endl;
  }
}

void LoadALICE(){

  for(int ich=0;ich<2;ich++){

    gpolALICE[ich]   = new TGraphErrors(2);
    gSEpolALICE[ich] = new TGraphAsymmErrors(2);
    gpolALICE  [ich]->SetMarkerStyle(34-ich*6);
    gpolALICE  [ich]->SetMarkerSize(1.0);
    gpolALICE  [ich]->SetMarkerColor(kBlue);
    gpolALICE  [ich]->SetLineColor(kBlue);

    gSEpolALICE[ich]->SetFillStyle(0);
    gSEpolALICE[ich]->SetFillColor(kBlue);
    gSEpolALICE[ich]->SetLineColor(kBlue);
    gSEpolALICE[ich]->SetLineWidth(1);
  }

  float fPHlam[2];
  float fPHalam[2];
  fPHlam[0]=0.08;//2.76TeV
  fPHlam[1]=-0.13;//5.02TeV
  fPHalam[0]=-0.05;//2.76TeV
  fPHalam[1]=0.14;//5.02TeV

  float ePHlam[2];
  float ePHalam[2];
  ePHlam[0]=0.10;
  ePHlam[1]=0.11;
  ePHalam[0]=0.10;
  ePHalam[1]=0.12;

  float sePHlam[2];
  float sePHalam[2];
  sePHlam[0]=0.04;
  sePHlam[1]=0.04;
  sePHalam[0]=0.03;
  sePHalam[1]=0.03;

  float fx[2] = {2760, 5020};
  double ScaleF=0.642/0.732;
  for(int iene=0;iene<2;iene++){
    //scaling
    fPHlam[iene]   *= ScaleF;
    fPHalam[iene]  *= ScaleF;
    ePHlam[iene]   *= ScaleF;
    ePHalam[iene]  *= ScaleF;
    sePHlam[iene]  *= ScaleF;
    sePHalam[iene] *= ScaleF;

    gpolALICE[0]   -> SetPoint     (iene,  fx[iene]+500*iene,  fPHlam[iene]);
    gpolALICE[0]   -> SetPointError(iene,  0,                  ePHlam[iene]);
    gSEpolALICE[0] -> SetPoint     (iene,  fx[iene]+500*iene,  fPHlam[iene]);
    gSEpolALICE[0] -> SetPointError(iene,  50, 70, sePHlam[iene], sePHlam[iene]);

    gpolALICE[1]   -> SetPoint     (iene,  fx[iene]+500*iene,  fPHalam[iene]);
    gpolALICE[1]   -> SetPointError(iene,  0,                  ePHalam[iene]);
    gSEpolALICE[1] -> SetPoint     (iene,  fx[iene]+500*iene,  fPHalam[iene]);
    gSEpolALICE[1] -> SetPointError(iene,  150, 200, sePHalam[iene], sePHalam[iene]);
  }


}

void LoadRuRu200GeV(){

  for(int ich=0;ich<2;ich++){//0 is Lambda , 1 is Lambdabar
    gpolRuRu200GeV[ich]   = new TGraphErrors(2);
    gSEpolRuRu200GeV[ich] = new TGraphAsymmErrors(2);
    gpolRuRu200GeV  [ich]->SetMarkerStyle(21+ich*4);
    gpolRuRu200GeV  [ich]->SetMarkerSize(0.8);
    gpolRuRu200GeV  [ich]->SetMarkerColor(kRed);
    gpolRuRu200GeV  [ich]->SetLineColor(kRed);

    gSEpolRuRu200GeV[ich]->SetFillStyle(0);
    gSEpolRuRu200GeV[ich]->SetFillColor(kRed);
    gSEpolRuRu200GeV[ich]->SetLineColor(kRed);
    gSEpolRuRu200GeV[ich]->SetLineWidth(1);
  }

  float fPHlam=0.210;//Lambda
  float fPHalam=0.221;//Lambdabar
  float ePHlam = 0.069 ;
  float ePHalam = 0.077;
  float sePHlam=0.016;
  float sePHalam=0.014;

  float fx1 = 250, fx2 = 270;

  gpolRuRu200GeV[0]   -> SetPoint     (0,  fx1,  fPHlam);
  gpolRuRu200GeV[0]   -> SetPointError(0,  0,    ePHlam);
  gSEpolRuRu200GeV[0] -> SetPoint     (0,  fx1,  fPHlam);
  gSEpolRuRu200GeV[0] -> SetPointError(0,  15, 20, sePHlam, sePHlam);

  gpolRuRu200GeV[1]   -> SetPoint     (0,  fx2,  fPHalam);
  gpolRuRu200GeV[1]   -> SetPointError(0,  0,    ePHalam);
  gSEpolRuRu200GeV[1] -> SetPoint     (0,  fx2,  fPHalam);
  gSEpolRuRu200GeV[1] -> SetPointError(0,  15, 20, sePHalam, sePHalam);

}


void LoadZrZr200GeV(){

  for(int ich=0;ich<2;ich++){//0 is Lambda , 1 is Lambdabar
    gpolZrZr200GeV[ich]   = new TGraphErrors(2);
    gSEpolZrZr200GeV[ich] = new TGraphAsymmErrors(2);
    gpolZrZr200GeV  [ich]->SetMarkerStyle(33-ich*6);
    gpolZrZr200GeV  [ich]->SetMarkerSize(1.0);
    gpolZrZr200GeV  [ich]->SetMarkerColor(kRed);
    gpolZrZr200GeV  [ich]->SetLineColor(kRed);

    gSEpolZrZr200GeV[ich]->SetFillStyle(0);
    gSEpolZrZr200GeV[ich]->SetFillColor(kRed);
    gSEpolZrZr200GeV[ich]->SetLineColor(kRed);
    gSEpolZrZr200GeV[ich]->SetLineWidth(1);
  }

  float fPHlam=0.279;//Lambda
  float fPHalam=0.176;//Lambdabar
  float ePHlam = 0.062 ;
  float ePHalam = 0.069;
  float sePHlam=0.020;
  float sePHalam=0.010;

  float fx1 = 290, fx2 = 310;

  gpolZrZr200GeV[0]   -> SetPoint     (0,  fx1,  fPHlam);
  gpolZrZr200GeV[0]   -> SetPointError(0,  0,    ePHlam);
  gSEpolZrZr200GeV[0] -> SetPoint     (0,  fx1,  fPHlam);
  gSEpolZrZr200GeV[0] -> SetPointError(0,  15, 20, sePHlam, sePHlam);

  gpolZrZr200GeV[1]   -> SetPoint     (0,  fx2,  fPHalam);
  gpolZrZr200GeV[1]   -> SetPointError(0,  0,    ePHalam);
  gSEpolZrZr200GeV[1] -> SetPoint     (0,  fx2,  fPHalam);
  gSEpolZrZr200GeV[1] -> SetPointError(0,  15, 20, sePHalam, sePHalam);

}

void SetAxis( TH1F *ax, float xff, float yoff ){
  ax->SetTitleFont(43,"X");
  ax->SetTitleSize(18,"X");
  ax->SetTitleOffset(xff,"X");
  ax->SetTitleFont(43,"Y");
  ax->SetTitleSize(18,"Y");
  ax->SetTitleOffset(yoff,"Y");
  ax->SetNdivisions(505,"X");
  ax->SetNdivisions(505,"Y");
  ax->SetLabelFont(43,"X");
  ax->SetLabelFont(43,"Y");
  ax->SetLabelSize(14,"X");
  ax->SetLabelSize(14,"Y");
}

void SetLegend( TLegend *l ){

  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextFont(43);
  l->SetTextSize(11);
}
