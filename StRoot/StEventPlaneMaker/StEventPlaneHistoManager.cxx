#include <TH2F.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneHistoManager.h"

ClassImp(StEventPlaneHistoManager)

//-------------------------------------------------------------------------------------------

StEventPlaneHistoManager::StEventPlaneHistoManager(int beamType) : mType(beamType)
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
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	string HistName = Form("h_mZdcGainCorr%s%sSlat%d",str_mEastWest[iEastWest].c_str(),str_mVertHori[iVertHori].c_str(),iSlat);
	h_mZdcGainCorr[iEastWest][iVertHori][iSlat] = new TH2F(HistName.c_str(),HistName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,5000,-4.5,4995.5);
      }
    }
  }
}

void StEventPlaneHistoManager::fillZdcGainCorr(int iEastWest, int iVertHori, int iSlat, int runIndex, double zdcsmd)
{
  h_mZdcGainCorr[iEastWest][iVertHori][iSlat]->Fill((double)runIndex,zdcsmd);
}

void StEventPlaneHistoManager::writeZdcGainCorr()
{
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	h_mZdcGainCorr[iEastWest][iVertHori][iSlat]->Write();
      }
    }
  }
}

// raw EP
void StEventPlaneHistoManager::initZdcRawEP()
{
  for(int iCent = 0; iCent < 9; ++iCent)
  {
    string HistName = Form("h_mZdcRawEastCent%d",iCent);
    h_mZdcRawEast[iCent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5);
    HistName = Form("h_mZdcRawWestCent%d",iCent);
    h_mZdcRawWest[iCent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5);
    HistName = Form("h_mZdcRawFullCent%d",iCent);
    h_mZdcRawFull[iCent] = new TH2F(HistName.c_str(),HistName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5);
  }
}

void StEventPlaneHistoManager::fillZdcRawEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcRawEast[Cent9]->Fill(PsiEast,runIndex);
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcRawWest[Cent9]->Fill(PsiWest,runIndex);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcRawFull[Cent9]->Fill(PsiFull,runIndex);
}

void StEventPlaneHistoManager::writeZdcRawEP()
{
  for(int iCent = 0; iCent < 9; ++iCent)
  {
    h_mZdcRawEast[iCent]->Write();
    h_mZdcRawWest[iCent]->Write();
    h_mZdcRawFull[iCent]->Write();
  }
}

//-------------------------------------------------------------------------------------------
