#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TCanvas.h"

#include "../../Utility/include/StSpinAlignmentCons.h"

void getTpcReCtrPar(int beamType = 0)
{
  const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0
  const int mNumCentrality = 9;

  string inputFile = Form("../../data/%s/EventPlaneMaker/file_ReCenterPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  TFile *file_InPut = TFile::Open(inputFile.c_str());
  if(!file_InPut->IsOpen()) cout << "inputFile: " << inputFile.c_str() << "is problematic" << endl;
  cout << "inputFile sets to: " << inputFile.c_str() << endl;

  // ReCenter Correction | x axis is runIndex, y axis is Centrality
  TProfile2D *p_mTpcQ2ReCtrXEast[mNumVzBin]; // 2nd EP
  TProfile2D *p_mTpcQ2ReCtrYEast[mNumVzBin];
  TProfile2D *p_mTpcQ2ReCtrXWest[mNumVzBin];
  TProfile2D *p_mTpcQ2ReCtrYWest[mNumVzBin];

  TProfile2D *p_mTpcQ3ReCtrXEast[mNumVzBin]; // 3rd EP
  TProfile2D *p_mTpcQ3ReCtrYEast[mNumVzBin];
  TProfile2D *p_mTpcQ3ReCtrXWest[mNumVzBin];
  TProfile2D *p_mTpcQ3ReCtrYWest[mNumVzBin];
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mTpcQ2ReCtrXEastVz%d",iVz); // 2nd EP
    p_mTpcQ2ReCtrXEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCtrYEastVz%d",iVz);
    p_mTpcQ2ReCtrYEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mTpcQ2ReCtrXWestVz%d",iVz);
    p_mTpcQ2ReCtrXWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcQ2ReCtrYWestVz%d",iVz);
    p_mTpcQ2ReCtrYWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mTpcQ3ReCtrXEastVz%d",iVz); // 3rd EP
    p_mTpcQ3ReCtrXEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcQ3ReCtrYEastVz%d",iVz);
    p_mTpcQ3ReCtrYEast[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());

    proName = Form("p_mTpcQ3ReCtrXWestVz%d",iVz);
    p_mTpcQ3ReCtrXWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
    proName = Form("p_mTpcQ3ReCtrYWestVz%d",iVz);
    p_mTpcQ3ReCtrYWest[iVz] = (TProfile2D*)file_InPut->Get(proName.c_str());
  }

  string outputFile = Form("../../Utility/EventPlaneMaker/%s/ReCenterPar/file_TpcReCenterPar_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile: " << outputFile.c_str() << endl;
  TFile *file_OutPut = new TFile(outputFile.c_str(),"RECREATE");
  file_OutPut->cd();
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mTpcQ2ReCtrXEast[iVz]->Write();
    p_mTpcQ2ReCtrYEast[iVz]->Write();
    p_mTpcQ2ReCtrXWest[iVz]->Write();
    p_mTpcQ2ReCtrYWest[iVz]->Write();

    p_mTpcQ3ReCtrXEast[iVz]->Write();
    p_mTpcQ3ReCtrYEast[iVz]->Write();
    p_mTpcQ3ReCtrXWest[iVz]->Write();
    p_mTpcQ3ReCtrYWest[iVz]->Write();
  }
  file_OutPut->Close();

  // Event Plane Distribution | x axis is runIndex, y axis is EP angle
  TH2F *h_mTpcEp2RawEast[mNumCentrality]; // 2nd raw EP
  TH2F *h_mTpcEp2RawWest[mNumCentrality];
  TH2F *h_mTpcEp2RawCorr[mNumCentrality]; // Psi2East vs Psi2West

  TH2F *h_mTpcEp3RawEast[mNumCentrality]; // 3rd raw EP
  TH2F *h_mTpcEp3RawWest[mNumCentrality];
  TH2F *h_mTpcEp3RawCorr[mNumCentrality]; // Psi3East vs Psi3West
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mTpcEp2RawEastCent%d",iCent); // 2nd EP
    h_mTpcEp2RawEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2RawWestCent%d",iCent);
    h_mTpcEp2RawWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp2RawCorrCent%d",iCent);
    h_mTpcEp2RawCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());

    histName = Form("h_mTpcEp3RawEastCent%d",iCent); // 3rd EP
    h_mTpcEp3RawEast[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3RawWestCent%d",iCent);
    h_mTpcEp3RawWest[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
    histName = Form("h_mTpcEp3RawCorrCent%d",iCent);
    h_mTpcEp3RawCorr[iCent] = (TH2F*)file_InPut->Get(histName.c_str());
  }

  TCanvas *c_TpcEpRawDist[mNumCentrality];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string canvName = Form("c_TpcEpRawDistCent%d",iCent);
    c_TpcEpRawDist[iCent] = new TCanvas(canvName.c_str(),canvName.c_str(),10,10,1500,1000);
    c_TpcEpRawDist[iCent]->Divide(3,2);
    c_TpcEpRawDist[iCent]->cd(1); h_mTpcEp2RawEast[iCent]->ProjectionY()->Draw();
    c_TpcEpRawDist[iCent]->cd(2); h_mTpcEp2RawWest[iCent]->ProjectionY()->Draw();
    c_TpcEpRawDist[iCent]->cd(3); h_mTpcEp2RawCorr[iCent]->Draw("colz");

    c_TpcEpRawDist[iCent]->cd(4); h_mTpcEp3RawEast[iCent]->ProjectionY()->Draw();
    c_TpcEpRawDist[iCent]->cd(5); h_mTpcEp3RawWest[iCent]->ProjectionY()->Draw();
    c_TpcEpRawDist[iCent]->cd(6); h_mTpcEp3RawCorr[iCent]->Draw("colz");

    std::string figName = Form("../../figures/%s/EventPlaneMaker/TpcRawEpCent%d_%s.pdf",globCons::str_mBeamType[beamType].c_str(),iCent,globCons::str_mBeamType[beamType].c_str());
    c_TpcEpRawDist[iCent]->SaveAs(figName.c_str());
  }

  string outputFileRawEp = Form("../../data/%s/EventPlaneMaker/file_TpcRawEpDist_%s.root",globCons::str_mBeamType[beamType].c_str(),globCons::str_mBeamType[beamType].c_str());
  cout << "outputFile of Raw EP: " << outputFileRawEp.c_str() << endl;
  TFile *file_OutPutRawEp = new TFile(outputFileRawEp.c_str(),"RECREATE");
  file_OutPutRawEp->cd();
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mTpcEp2RawEast[iCent]->Write();
    h_mTpcEp2RawWest[iCent]->Write();
    h_mTpcEp2RawCorr[iCent]->Write();

    h_mTpcEp3RawEast[iCent]->Write();
    h_mTpcEp3RawWest[iCent]->Write();
    h_mTpcEp3RawCorr[iCent]->Write();
  }
  file_OutPutRawEp->Close();
}
