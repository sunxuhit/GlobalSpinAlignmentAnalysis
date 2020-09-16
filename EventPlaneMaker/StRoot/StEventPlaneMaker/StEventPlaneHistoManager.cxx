#include "StRoot/StEventPlaneMaker/StEventPlaneHistoManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

#include <TH2F.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>

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
	h_mZdcGainCorr[i_eastwest][i_verthori][i_slat] = new TH2F(HistName.c_str(),HistName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,1000,-4.5,995.5);
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
//-------------------------------------------------------------------------------------------
