#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getZdcReCtrPar(int beamType = 0)
{
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_ReCenterPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // ReCenter Correction | x axis is runIndex, y axis is Centrality
  TProfile2D *p_mZdcQ1ReCtrVertEast[mNumVzBin];
  TProfile2D *p_mZdcQ1ReCtrHoriEast[mNumVzBin];
  TProfile2D *p_mZdcQ1ReCtrVertWest[mNumVzBin];
  TProfile2D *p_mZdcQ1ReCtrHoriWest[mNumVzBin];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mZdcQ1ReCtrVertEastVz%d",iVz);
    p_mZdcQ1ReCtrVertEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mZdcQ1ReCtrHoriEastVz%d",iVz);
    p_mZdcQ1ReCtrHoriEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mZdcQ1ReCtrVertWestVz%d",iVz);
    p_mZdcQ1ReCtrVertWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mZdcQ1ReCtrHoriWestVz%d",iVz);
    p_mZdcQ1ReCtrHoriWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ReCenterPar/file_ZdcReCenterPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mZdcQ1ReCtrVertEast[iVz]->Write();
    p_mZdcQ1ReCtrHoriEast[iVz]->Write();
    p_mZdcQ1ReCtrVertWest[iVz]->Write();
    p_mZdcQ1ReCtrHoriWest[iVz]->Write();
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mZdcEp1RawEast[mNumCentrality]; // raw EP
  TH2F *h_mZdcEp1RawWest[mNumCentrality];
  TH2F *h_mZdcEp1RawFull[mNumCentrality]; // Qwest-QEast
  TH2F *h_mZdcEp1RawCorr[mNumCentrality]; // Psi1East vs Psi1West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEp1RawEastCent%d",iCent);
    h_mZdcEp1RawEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEp1RawWestCent%d",iCent);
    h_mZdcEp1RawWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mZdcEp1RawFullCent%d",iCent);
    h_mZdcEp1RawFull[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mZdcEp1RawCorrCent%d",iCent);
    h_mZdcEp1RawCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  TCanvas *c_ZdcEp1RawDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_EpdEp1RawDistCent%d",iCent);
    c_ZdcEp1RawDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1600,400);
    c_ZdcEp1RawDist[iCent]->Divide(4,1);
    c_ZdcEp1RawDist[iCent]->cd(1); h_mZdcEp1RawEast[iCent]->ProjectionY()->Draw();
    c_ZdcEp1RawDist[iCent]->cd(2); h_mZdcEp1RawWest[iCent]->ProjectionY()->Draw();
    c_ZdcEp1RawDist[iCent]->cd(3); h_mZdcEp1RawFull[iCent]->ProjectionY()->Draw();
    c_ZdcEp1RawDist[iCent]->cd(4); h_mZdcEp1RawCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/ZdcRawEpCent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_ZdcEp1RawDist[iCent]->SaveAs(figName.c_str());
  }

  string outputFileRawEp = Form("../../data/%s/EventPlaneMaker/file_ZdcRawEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileRawEp.c_str() << endl;
  TFile *file_OutPutRawEp = new TFile(outputFileRawEp.c_str(),"RECREATE");
  file_OutPutRawEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEp1RawEast[iCent]->Write();
    h_mZdcEp1RawWest[iCent]->Write();
    h_mZdcEp1RawFull[iCent]->Write();
    h_mZdcEp1RawCorr[iCent]->Write();
  }
  file_OutPutRawEp->Close();
}
