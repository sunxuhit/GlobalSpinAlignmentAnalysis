#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

#define _RawEp_ 1
#define _WgtEp_ 0
#define _ReCtrEp_ 0
#define _ShiftEp_ 0
#define _ShiftFullEp_ 0

void getEpdEpDist(int beamType = 0)
{
  const int mNumCentrality = 9;

#if _RawEp_
  string inputFileRawEp = Form("../../data/%s/EventPlaneMaker/file_%s_GainCorr.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPutRawEp = TFile::Open(inputFileRawEp.c_str());
  if(!file_InPutRawEp->IsOpen()) cout << "inputFile of Raw EP: " << inputFileRawEp.c_str() << "is problematic" << endl;
  cout << "inputFile of Raw EP sets to: " << inputFileRawEp.c_str() << endl;

  TH2F *h_mEpdEp1RawEast[mNumCentrality]; // 1st raw EP
  TH2F *h_mEpdEp1RawWest[mNumCentrality];
  TH2F *h_mEpdEp1RawFull[mNumCentrality];
  TH2F *h_mEpdEp1RawCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1RawEastCent%d",iCent); // 2nd EP
    h_mEpdEp1RawEast[iCent] = (TH2F*)file_InPutRawEp->Get(histName.c_str());
    histName = Form("h_mEpdEp1RawWestCent%d",iCent);
    h_mEpdEp1RawWest[iCent] = (TH2F*)file_InPutRawEp->Get(histName.c_str());
    histName = Form("h_mEpdEp1RawFullCent%d",iCent);
    h_mEpdEp1RawFull[iCent] = (TH2F*)file_InPutRawEp->Get(histName.c_str());
    histName = Form("h_mEpdEp1RawCorrCent%d",iCent);
    h_mEpdEp1RawCorr[iCent] = (TH2F*)file_InPutRawEp->Get(histName.c_str());
  }

  TCanvas *c_EpdEp1RawDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_EpdEp1RawDistCent%d",iCent);
    c_EpdEp1RawDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1600,400);
    c_EpdEp1RawDist[iCent]->Divide(4,1);
    c_EpdEp1RawDist[iCent]->cd(1); h_mEpdEp1RawEast[iCent]->ProjectionY()->Draw();
    c_EpdEp1RawDist[iCent]->cd(2); h_mEpdEp1RawWest[iCent]->ProjectionY()->Draw();
    c_EpdEp1RawDist[iCent]->cd(3); h_mEpdEp1RawFull[iCent]->ProjectionY()->Draw();
    c_EpdEp1RawDist[iCent]->cd(4); h_mEpdEp1RawCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/EpdRawEpCent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_EpdEp1RawDist[iCent]->SaveAs(figName.c_str());
  }

  string outputFileRawEp = Form("../../data/%s/EventPlaneMaker/file_%s_EpdRawEpDist.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileRawEp.c_str() << endl;
  TFile *file_OutPutRawEp = new TFile(outputFileRawEp.c_str(),"RECREATE");
  file_OutPutRawEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1RawEast[iCent]->Write();
    h_mEpdEp1RawWest[iCent]->Write();
    h_mEpdEp1RawFull[iCent]->Write();
    h_mEpdEp1RawCorr[iCent]->Write();
  }
  file_OutPutRawEp->Close();
#endif

#if _WgtEp_
  TH2F *h_mEpdEp1WgtEast[mNumCentrality]; // 1st weighted EP
  TH2F *h_mEpdEp1WgtWest[mNumCentrality];
  TH2F *h_mEpdEp1WgtFull[mNumCentrality];
  TH2F *h_mEpdEp1WgtCorr[mNumCentrality]; // Psi1East vs Psi1West
#endif

#if _ReCtrEp_
  TH2F *h_mEpdEp1ReCtrEast[mNumCentrality]; // 1st recenter EP
  TH2F *h_mEpdEp1ReCtrWest[mNumCentrality];
  TH2F *h_mEpdEp1ReCtrFull[mNumCentrality];
  TH2F *h_mEpdEp1ReCtrCorr[mNumCentrality]; // Psi1East vs Psi1West
#endif

#if _ShiftEp_
  TH2F *h_mEpdEp1ShiftEast[mNumCentrality]; // 1st shift EP
  TH2F *h_mEpdEp1ShiftWest[mNumCentrality];
  TH2F *h_mEpdEp1ShiftFull[mNumCentrality];
  TH2F *h_mEpdEp1ShiftCorr[mNumCentrality]; // Psi1East vs Psi1West
#endif

#if _ShiftFullEp_
  TH2F *h_mEpdEp1ShiftFullCorr[mNumCentrality]; // 1st shift full EP
#endif
}
