#include <string>
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"
#include "../StRoot/StRunQAMaker/StRunQACons.h"

using namespace std;

static const string CutsQA[2] = {"Before","After"};

void findBadRunList(int energy = 0)
{
  string JobId = "low";
  string inputfile = Form("/star/u/sunxuhit/AuAu%s/SpinAlignment/RunQA/merged_file/file_%s_RunQA_%s.root",runQA::mBeamEnergy[energy].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
  TFile *File_InPut = TFile::Open(inputfile.c_str());

  TProfile *p_mQA_RefMult[2]; // 0: before cuts | 1: after cuts
  TProfile *p_mQA_gRefMult[2];
  TProfile *p_mQA_ZdcX[2];
  TProfile *p_mQA_Vz[2];
  TProfile *p_mQA_Vr[2];
  // TProfile *p_mQA_Res1[2];
  // TProfile *p_mQA_Res2[2];
  // TProfile *p_mQA_Res3[2];
  // TProfile *p_mQA_v2[2];
  // TProfile *p_mQA_v3[2];

  TProfile *p_mQA_gDca[2];
  TProfile *p_mQA_nHitsFit[2];
  TProfile *p_mQA_PrimPt[2];
  TProfile *p_mQA_PrimEta[2];
  TProfile *p_mQA_PrimPhi[2];
  TProfile *p_mQA_GlobPt[2];
  TProfile *p_mQA_GlobEta[2];
  TProfile *p_mQA_GlobPhi[2];

  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    std::string ProName = Form("p_mQA_RefMult_%s",CutsQA[i_cut].c_str());
    p_mQA_RefMult[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_RefMult[i_cut]->SetLineColor(i_cut+1);
    p_mQA_RefMult[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_RefMult[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_RefMult[i_cut]->GetYaxis()->SetTitle("<refMult>");

    ProName = Form("p_mQA_gRefMult_%s",CutsQA[i_cut].c_str());
    p_mQA_gRefMult[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_gRefMult[i_cut]->SetLineColor(i_cut+1);
    p_mQA_gRefMult[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_gRefMult[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_gRefMult[i_cut]->GetYaxis()->SetTitle("<gRefMult>");

    ProName = Form("p_mQA_ZdcX_%s",CutsQA[i_cut].c_str());
    p_mQA_ZdcX[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_ZdcX[i_cut]->SetLineColor(i_cut+1);
    p_mQA_ZdcX[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_ZdcX[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_ZdcX[i_cut]->GetYaxis()->SetTitle("<ZdcX>");

    ProName = Form("p_mQA_Vz_%s",CutsQA[i_cut].c_str());
    p_mQA_Vz[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_Vz[i_cut]->SetLineColor(i_cut+1);
    p_mQA_Vz[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_Vz[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_Vz[i_cut]->GetYaxis()->SetTitle("<Vz>");

    ProName = Form("p_mQA_Vr_%s",CutsQA[i_cut].c_str());
    p_mQA_Vr[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_Vr[i_cut]->SetLineColor(i_cut+1);
    p_mQA_Vr[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_Vr[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_Vr[i_cut]->GetYaxis()->SetTitle("<Vr>");

    ProName = Form("p_mQA_gDca_%s",CutsQA[i_cut].c_str());
    p_mQA_gDca[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_gDca[i_cut]->SetLineColor(i_cut+1);
    p_mQA_gDca[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_gDca[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_gDca[i_cut]->GetYaxis()->SetTitle("<gDca>");

    ProName = Form("p_mQA_nHitsFit_%s",CutsQA[i_cut].c_str());
    p_mQA_nHitsFit[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_nHitsFit[i_cut]->SetLineColor(i_cut+1);
    p_mQA_nHitsFit[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_nHitsFit[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_nHitsFit[i_cut]->GetYaxis()->SetTitle("<nHitsFit>");

    ProName = Form("p_mQA_PrimPt_%s",CutsQA[i_cut].c_str());
    p_mQA_PrimPt[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_PrimPt[i_cut]->SetLineColor(i_cut+1);
    p_mQA_PrimPt[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_PrimPt[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_PrimPt[i_cut]->GetYaxis()->SetTitle("<primPt>");

    ProName = Form("p_mQA_PrimEta_%s",CutsQA[i_cut].c_str());
    p_mQA_PrimEta[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_PrimEta[i_cut]->SetLineColor(i_cut+1);
    p_mQA_PrimEta[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_PrimEta[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_PrimEta[i_cut]->GetYaxis()->SetTitle("<primEta>");

    ProName = Form("p_mQA_PrimPhi_%s",CutsQA[i_cut].c_str());
    p_mQA_PrimPhi[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_PrimPhi[i_cut]->SetLineColor(i_cut+1);
    p_mQA_PrimPhi[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_PrimPhi[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_PrimPhi[i_cut]->GetYaxis()->SetTitle("<primPhi>");

    ProName = Form("p_mQA_GlobPt_%s",CutsQA[i_cut].c_str());
    p_mQA_GlobPt[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_GlobPt[i_cut]->SetLineColor(i_cut+1);
    p_mQA_GlobPt[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_GlobPt[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_GlobPt[i_cut]->GetYaxis()->SetTitle("<globPt>");

    ProName = Form("p_mQA_GlobEta_%s",CutsQA[i_cut].c_str());
    p_mQA_GlobEta[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_GlobEta[i_cut]->SetLineColor(i_cut+1);
    p_mQA_GlobEta[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_GlobEta[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_GlobEta[i_cut]->GetYaxis()->SetTitle("<globEta>");

    ProName = Form("p_mQA_GlobPhi_%s",CutsQA[i_cut].c_str());
    p_mQA_GlobPhi[i_cut] = (TProfile*)File_InPut->Get(ProName.c_str());
    p_mQA_GlobPhi[i_cut]->SetLineColor(i_cut+1);
    p_mQA_GlobPhi[i_cut]->SetMarkerColor(i_cut+1);
    p_mQA_GlobPhi[i_cut]->GetXaxis()->SetTitle("runIndex");
    p_mQA_GlobPhi[i_cut]->GetYaxis()->SetTitle("<globPhi>");
  }

  TCanvas *c_EvtLvlRunQA[2];
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    string CanName = Form("c_EvtLvlRunQA_%s",CutsQA[i_cut].c_str());
    c_EvtLvlRunQA[i_cut] = new TCanvas(CanName.c_str(),CanName.c_str(),10,10,1500,600);
    c_EvtLvlRunQA[i_cut]->Divide(5,2);
    for(int i_pad = 0; i_pad < 10; ++i_pad)
    {
      c_EvtLvlRunQA[i_cut]->cd(i_pad+1);
      c_EvtLvlRunQA[i_cut]->cd(i_pad+1)->SetLeftMargin(0.1);
      c_EvtLvlRunQA[i_cut]->cd(i_pad+1)->SetRightMargin(0.1);
      c_EvtLvlRunQA[i_cut]->cd(i_pad+1)->SetBottomMargin(0.1);
      c_EvtLvlRunQA[i_cut]->cd(i_pad+1)->SetGrid(0,0);
      c_EvtLvlRunQA[i_cut]->cd(i_pad+1)->SetTicks(1,1);
    }

    c_EvtLvlRunQA[i_cut]->cd(1);
    p_mQA_RefMult[i_cut]->Draw("pE");

    c_EvtLvlRunQA[i_cut]->cd(2);
    p_mQA_gRefMult[i_cut]->Draw("pE");

    c_EvtLvlRunQA[i_cut]->cd(3);
    p_mQA_ZdcX[i_cut]->Draw("pE");

    c_EvtLvlRunQA[i_cut]->cd(4);
    p_mQA_Vz[i_cut]->Draw("pE");

    c_EvtLvlRunQA[i_cut]->cd(5);
    p_mQA_Vr[i_cut]->Draw("pE");;

    // string FigName = Form("c_EvtLvlRunQA_%s_%s_%s.png",CutsQA[i_cut].c_str(),runQA::mBeamEnergy[energy].c_str(),JobId.c_str());
    // c_EvtLvlRunQA[i_cut]->SaveAs(FigName.c_str());
  }
}
