#include "StRoot/StVecMesonMaker/StVecMesonHistoManager.h"

#include <TH2F.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>

ClassImp(StVecMesonHistoManager)

//-------------------------------------------------------------------------------------------

StVecMesonHistoManager::StVecMesonHistoManager()
{
}

//-------------------------------------------------------------------------------------------

StVecMesonHistoManager::~StVecMesonHistoManager()
{
  /* */
}

//-------------------------------------------------------------------------------------------
void StVecMesonHistoManager::Init_EventQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    std::string HistName = Form("h_mRefMult_%s",mCutsQA[i_cut].c_str());
    h_mRefMult[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5);

    HistName = Form("h_mCentrality9_%s",mCutsQA[i_cut].c_str());
    h_mCentrality9[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),12,-1.5,10.5);

    HistName = Form("h_mRefMultTofMatch_%s",mCutsQA[i_cut].c_str());
    h_mRefMultTofMatch[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5,2000,-0.5,1999.5);

    HistName = Form("h_mRefMultTofHits_%s",mCutsQA[i_cut].c_str());
    h_mRefMultTofHits[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5,2000,-0.5,1999.5);

    HistName = Form("h_mVertexXY_%s",mCutsQA[i_cut].c_str());
    h_mVertexXY[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05,201,-10.05,10.05);

    HistName = Form("h_mVertexZ_%s",mCutsQA[i_cut].c_str());
    h_mVertexZ[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),201,-100.5,100.5);

    HistName = Form("h_mVzVzVpd_%s",mCutsQA[i_cut].c_str());
    h_mVzVzVpd[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),201,-100.5,100.5,201,-100.5,100.5);

    HistName = Form("h_mDiffVzVzVpd_%s",mCutsQA[i_cut].c_str());
    h_mDiffVzVzVpd[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05);
  }
}

void StVecMesonHistoManager::Fill_EventQA_RefMult(float refMult, float cent9, float tofHits, float tofMatch, int cutSelection)
{
  h_mRefMult[cutSelection]->Fill(refMult);
  h_mCentrality9[cutSelection]->Fill(cent9);
  h_mRefMultTofMatch[cutSelection]->Fill(refMult,tofMatch);
  h_mRefMultTofHits[cutSelection]->Fill(refMult,tofHits);
}

void StVecMesonHistoManager::Fill_EventQA_Vertex(float vx, float vy, float vz, float vzVpd, int cutSelection)
{
  h_mVertexXY[cutSelection]->Fill(vx,vy);
  h_mVertexZ[cutSelection]->Fill(vz);
  h_mVzVzVpd[cutSelection]->Fill(vz,vzVpd);
  h_mDiffVzVzVpd[cutSelection]->Fill(vz-vzVpd);
}

void StVecMesonHistoManager::Write_EventQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    h_mRefMult[i_cut]->Write();
    h_mCentrality9[i_cut]->Write();
    h_mRefMultTofMatch[i_cut]->Write();
    h_mRefMultTofHits[i_cut]->Write();

    h_mVertexXY[i_cut]->Write();
    h_mVertexZ[i_cut]->Write();
    h_mVzVzVpd[i_cut]->Write();
    h_mDiffVzVzVpd[i_cut]->Write();
  }
}

//-------------------------------------------------------------------------------------------
/*
void StVecMesonHistoManager::InitEP()
{
  h_mEastReCenter = new TH1F("h_mEastReCenter","h_mEastReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mWestReCenter = new TH1F("h_mWestReCenter","h_mWestReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mRanAReCenter = new TH1F("h_mRanAReCenter","h_mRanAReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mRanBReCenter = new TH1F("h_mRanBReCenter","h_mRanBReCenter",360,-TMath::Pi(),TMath::Pi());
  h_mFullReCenter = new TH1F("h_mFullReCenter","h_mFullReCenter",360,-TMath::Pi(),TMath::Pi());

  h_mEastShift = new TH1F("h_mEastShift","h_mEastShift",360,-TMath::Pi(),TMath::Pi());
  h_mWestShift = new TH1F("h_mWestShift","h_mWestShift",360,-TMath::Pi(),TMath::Pi());
  h_mRanAShift = new TH1F("h_mRanAShift","h_mRanAShift",360,-TMath::Pi(),TMath::Pi());
  h_mRanBShift = new TH1F("h_mRanBShift","h_mRanBShift",360,-TMath::Pi(),TMath::Pi());
  h_mFullShift = new TH1F("h_mFullShift","h_mFullShift",360,-TMath::Pi(),TMath::Pi());
}

void StVecMesonHistoManager::FillEP_Eta(Float_t Psi2_East, Float_t Psi2_West)
{
  h_mEastRaw->Fill(Psi2_East);
  h_mWestRaw->Fill(Psi2_West);
}

void StVecMesonHistoManager::FillEP_Full(Float_t Psi2_Full)
{
  h_mFullRaw->Fill(Psi2_Full);
}

void StVecMesonHistoManager::FillEP_Sub(Float_t Psi2East_ReCenter, Float_t Psi2East_Shift, Float_t Psi2West_ReCenter, Float_t Psi2West_Shift)
{
  h_mEastReCenter->Fill(Psi2East_ReCenter);
  h_mEastShift->Fill(Psi2East_Shift);

  h_mWestReCenter->Fill(Psi2West_ReCenter);
  h_mWestShift->Fill(Psi2West_Shift);
}

void StVecMesonHistoManager::FillEP_Ran(Float_t Psi2RanA_ReCenter, Float_t Psi2RanA_Shift, Float_t Psi2RanB_ReCenter, Float_t Psi2RanB_Shift, Float_t Psi2Full_ReCenter, Float_t Psi2Full_Shift)
{
  h_mRanAReCenter->Fill(Psi2RanA_ReCenter);
  h_mRanAShift->Fill(Psi2RanA_Shift);

  h_mRanBReCenter->Fill(Psi2RanB_ReCenter);
  h_mRanBShift->Fill(Psi2RanB_Shift);

  h_mFullReCenter->Fill(Psi2Full_ReCenter);
  h_mFullShift->Fill(Psi2Full_Shift);
}

void StVecMesonHistoManager::WriteEP()
{
  h_mEastReCenter->Write();
  h_mWestReCenter->Write();
  h_mRanAReCenter->Write();
  h_mRanBReCenter->Write();
  h_mFullReCenter->Write();

  h_mEastShift->Write();
  h_mWestShift->Write();
  h_mRanAShift->Write();
  h_mRanBShift->Write();
  h_mFullShift->Write();
}
*/
