#include <TH2F.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>

#include "StRoot/StEventPlaneMaker/StEventPlaneHistoManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

ClassImp(StEventPlaneHistoManager)

//-------------------------------------------------------------------------------------------

StEventPlaneHistoManager::StEventPlaneHistoManager()
{
}

//-------------------------------------------------------------------------------------------

StEventPlaneHistoManager::~StEventPlaneHistoManager()
{
  /* */
}

//-------------------------------------------------------------------------------------------
// ZDC EP
void StEventPlaneHistoManager::initZdcGainCorr()
{
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	string HistName = Form("h_mZdcGainCorr%s%s_%d",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str(),i_slat);
	h_mZdcGainCorr[i_eastwest][i_verthori][i_slat] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,5000,-4.5,4995.5);
      }
    }
  }
}

void StEventPlaneHistoManager::fillZdcGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, float zdcsmd)
{
  h_mZdcGainCorr[i_eastwest][i_verthori][i_slat]->Fill((float)runIndex,zdcsmd);
}

void StEventPlaneHistoManager::writeZdcGainCorr()
{
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	h_mZdcGainCorr[i_eastwest][i_verthori][i_slat]->Write();
      }
    }
  }
}

// raw EP
void StEventPlaneHistoManager::initZdcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    string HistName = Form("h_mZdcRawEast_%d",i_cent);
    h_mZdcRawEast[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    HistName = Form("h_mZdcRawWest_%d",i_cent);
    h_mZdcRawWest[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
    HistName = Form("h_mZdcRawFull_%d",i_cent);
    h_mZdcRawFull[i_cent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5);
  }
}

void StEventPlaneHistoManager::fillZdcRawEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex)
{
  float PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcRawEast[Cent9]->Fill(PsiEast,runIndex);
  float PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcRawWest[Cent9]->Fill(PsiWest,runIndex);
  float PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcRawFull[Cent9]->Fill(PsiFull,runIndex);
}

void StEventPlaneHistoManager::writeZdcRawEP()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    h_mZdcRawEast[i_cent]->Write();
    h_mZdcRawWest[i_cent]->Write();
    h_mZdcRawFull[i_cent]->Write();
  }
}

//-------------------------------------------------------------------------------------------
