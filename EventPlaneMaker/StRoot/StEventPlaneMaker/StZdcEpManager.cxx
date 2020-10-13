// #include "StPicoDstMaker/StPicoDstMaker.h"
// #include "StPicoEvent/StPicoDst.h"
// #include "StPicoEvent/StPicoEvent.h"
// #include "StPicoEvent/StPicoTrack.h"
// #include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StMessMgr.h"

#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TF1.h"

Double_t Resolution_Full(Double_t *x_val, Double_t *par)
{
  Double_t y;
  Double_t chi = x_val[0];
  Double_t arg = chi*chi/4.0;
  Double_t norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

ClassImp(StZdcEpManager)

//---------------------------------------------------------------------------------

StZdcEpManager::StZdcEpManager(int energy)
{
  mEnergy = energy;
  clearZdcEp();
}

StZdcEpManager::~StZdcEpManager()
{
  /* */
}

void StZdcEpManager::clearZdcEp()
{
  mCent9 = -1;
  mRunIndex = -1;
  mVz_sign = -1;
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	mZdcSmd[i_eastwest][i_verthori][i_slat] = 0.0;
	mGainCorrFactor[i_eastwest][i_verthori][i_slat] = -999.9;
      }
    }
  }
  mCenterEastVertical   = -999.9;
  mCenterEastHorizontal = -999.9;
  mCenterWestVertical   = -999.9;
  mCenterWestHorizontal = -999.9;
  for(int i_cent = 0; i_cent < 9; ++i_cent) mResolution[i_cent] = 0.0;
}

void StZdcEpManager::initZdcEp(int Cent9, int RunIndex)
{
  mCent9 = Cent9;
  mRunIndex = RunIndex;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::setZdcSmd(int eastwest, int verthori, int slat, const float zdcsmd) 
{
  mZdcSmd[eastwest][verthori][slat] = (zdcsmd > 0.) ? zdcsmd : 0.;
}

float StZdcEpManager::getZdcSmd(int eastwest, int verthori, int slat)
{
  return mZdcSmd[eastwest][verthori][slat];
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readGainCorr()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/GainCorrPar/file_%s_ZdcGainCorrFac.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_GainCorrPar = TFile::Open(InPutFile.c_str());
  for(int i_eastwest = 0; i_eastwest < 2; ++i_eastwest)
  {
    for(int i_verthori = 0; i_verthori < 2; ++i_verthori)
    {
      for(int i_slat = 0; i_slat < 8; ++i_slat)
      {
	string HistName = Form("h_mGainCorrFactor%s%s_%d",recoEP::mEastWest[i_eastwest].c_str(),recoEP::mVertHori[i_verthori].c_str(),i_slat);
	TH1F *h_GainCorrFac = (TH1F*)mFile_GainCorrPar->Get(HistName.c_str());
	mGainCorrFactor[i_eastwest][i_verthori][i_slat] = h_GainCorrFac->GetBinContent(1);
	// cout << "i_eastwest = " << i_eastwest << ", i_verthori = " << i_verthori << ", i_slat = " << i_slat << ", mGainCorrFactor = " << mGainCorrFactor[i_eastwest][i_verthori][i_slat] << endl;
      }
    }
  }
}

void StZdcEpManager::setZdcSmdGainCorr(int eastwest, int verthori, int slat, const float zdcsmd)
{
  mZdcSmd[eastwest][verthori][slat] = (zdcsmd > 0.) ? zdcsmd/mGainCorrFactor[eastwest][verthori][slat] : 0.;
  // cout << "input zdc = " << zdcsmd << ", mGainCorrFactor = " << mGainCorrFactor[eastwest][verthori][slat] << ", GainCorred = " << mZdcSmd[eastwest][verthori][slat] << endl;
}

float StZdcEpManager::getZdcSmdGainCorr(int eastwest, int verthori, int slat)
{
  return mZdcSmd[eastwest][verthori][slat];
}

float StZdcEpManager::getPosition(int eastwest, int verthori, int slat, int mode)
{
  //get position of each slat

  float zdcsmd_vert[7] = {0.5,2,3.5,5,6.5,8,9.5};
  float zdcsmd_hori[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

  if(mode > 1) // with beam center corrected
  {
    setZdcSmdCenter();
    if(mCenterEastVertical < -100.0 || mCenterEastHorizontal < -100.0 || mCenterWestVertical < -100.0 || mCenterWestHorizontal < -100.0) 
    {
      cout << "Forgot Re-Center!!!!" << endl;
      return 0;
    }
    if(eastwest == 0 && verthori == 0) return zdcsmd_vert[slat]-mCenterEastVertical;
    if(eastwest == 1 && verthori == 0) return mCenterWestVertical-zdcsmd_vert[slat];
    if(eastwest == 0 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.)-mCenterEastHorizontal;
    if(eastwest == 1 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.)-mCenterWestHorizontal;
  }
  else // raw beam center returned
  {
    if(eastwest == 0 && verthori == 0) return zdcsmd_vert[slat];
    if(eastwest == 1 && verthori == 0) return -zdcsmd_vert[slat];
    if(eastwest == 0 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.);
    if(eastwest == 1 && verthori == 1) return zdcsmd_hori[slat]/sqrt(2.);
  }

  return 0;
}

TVector2 StZdcEpManager::getQEast(int mode) 
{
  TVector2 qVector(0.0,0.0);
  float qXsum = 0.; float qYsum = 0.;
  float qXwgt = 0.; float qYwgt = 0.;

  for(int i_vert = 0; i_vert < 7; ++i_vert) // vertical
  {
    qXsum += getPosition(0,0,i_vert,mode)*getZdcSmdGainCorr(0,0,i_vert);
    qXwgt += getZdcSmdGainCorr(0,0,i_vert);
  }
  for(int i_hori = 0; i_hori < 8; ++i_hori) // horizontal
  {
    qYsum += getPosition(0,1,i_hori,mode)*getZdcSmdGainCorr(0,1,i_hori);
    qYwgt += getZdcSmdGainCorr(0,1,i_hori);
  }

  if(qXwgt > 0.0 && qYwgt > 0.0) qVector.Set(qXsum/qXwgt,qYsum/qYwgt);
  if(mode > 2)  qVector = ApplyZdcSmdShiftCorrEast(qVector);

  return qVector;
}

TVector2 StZdcEpManager::getQWest(int mode)
{
  TVector2 qVector(0.0,0.0);
  float qXsum = 0.; float qYsum = 0.;
  float qXwgt = 0.; float qYwgt = 0.;

  for(int i_vert = 0; i_vert < 7; ++i_vert) // vertical
  {
    qXsum += getPosition(1,0,i_vert,mode)*getZdcSmdGainCorr(1,0,i_vert);
    qXwgt += getZdcSmdGainCorr(1,0,i_vert);
  }
  for(int i_hori = 0; i_hori < 8; ++i_hori) // horizontal
  {
    qYsum += getPosition(1,1,i_hori,mode)*getZdcSmdGainCorr(1,1,i_hori);
    qYwgt += getZdcSmdGainCorr(1,1,i_hori);
  }

  if(qXwgt > 0.0 && qYwgt > 0.0) qVector.Set(qXsum/qXwgt,qYsum/qYwgt);
  if(mode > 2) qVector = ApplyZdcSmdShiftCorrWest(qVector);

  return qVector;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readReCenterCorr()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/ReCenterParameter/file_%s_ZdcReCenterParameter.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_ReCenterPar = TFile::Open(InPutFile.c_str());

  string ProName;

  ProName = "p_mQEastVertical";
  p_mQEastVertical = (TProfile2D*)mFile_ReCenterPar->Get(ProName.c_str());
  ProName = "p_mQEastHorizontal";
  p_mQEastHorizontal = (TProfile2D*)mFile_ReCenterPar->Get(ProName.c_str());

  ProName = "p_mQWestVertical";
  p_mQWestVertical = (TProfile2D*)mFile_ReCenterPar->Get(ProName.c_str());
  ProName = "p_mQWestHorizontal";
  p_mQWestHorizontal = (TProfile2D*)mFile_ReCenterPar->Get(ProName.c_str());
}

void StZdcEpManager::setZdcSmdCenter()
{
  int binEastVertical = p_mQEastVertical->FindBin((double)mRunIndex,(double)mCent9);
  mCenterEastVertical = p_mQEastVertical->GetBinContent(binEastVertical);

  int binEastHorizontal = p_mQEastHorizontal->FindBin((double)mRunIndex,(double)mCent9);
  mCenterEastHorizontal = p_mQEastHorizontal->GetBinContent(binEastHorizontal);

  int binWestVertical = p_mQWestVertical->FindBin((double)mRunIndex,(double)mCent9);
  mCenterWestVertical = -1.0*p_mQWestVertical->GetBinContent(binWestVertical);

  int binWestHorizontal = p_mQWestHorizontal->FindBin((double)mRunIndex,(double)mCent9);
  mCenterWestHorizontal = p_mQWestHorizontal->GetBinContent(binWestHorizontal);
  // cout << "mCenterEastVertical = " << mCenterEastVertical << ", mCenterEastHorizontal = " << mCenterEastHorizontal << ", mCenterWestVertical = " << mCenterWestVertical << ", mCenterWestHorizontal = " << mCenterWestHorizontal << endl;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readShiftCorr()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/ShiftParameter/file_%s_ZdcShiftParameter.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_ShiftPar = TFile::Open(InPutFile.c_str());

  for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
  {
    string ProName;

    ProName = Form("p_mQEastCos_%d",i_shift);
    p_mQEastCos[i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
    ProName = Form("p_mQEastSin_%d",i_shift);
    p_mQEastSin[i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());

    ProName = Form("p_mQWestCos_%d",i_shift);
    p_mQWestCos[i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
    ProName = Form("p_mQWestSin_%d",i_shift);
    p_mQWestSin[i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
  }
}

TVector2 StZdcEpManager::ApplyZdcSmdShiftCorrEast(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  float Psi_ReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  float delta_Psi = 0.0;
  float Psi_Shift;

  for(Int_t i_shift = 0; i_shift < 20; ++i_shift) // Shift Order loop
  {
    int bin_Cos = p_mQEastCos[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Cos = p_mQEastCos[i_shift]->GetBinContent(bin_Cos);

    int bin_Sin = p_mQEastSin[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Sin = p_mQEastSin[i_shift]->GetBinContent(bin_Sin);

    delta_Psi += (2.0/((float)i_shift+1.0))*(-1.0*mean_Sin*TMath::Cos(((float)i_shift+1.0)*Psi_ReCenter)+mean_Cos*TMath::Sin(((float)i_shift+1.0)*Psi_ReCenter));
  }

  float Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  qVecShift.Set(TMath::Cos(Psi_Shift),TMath::Sin(Psi_Shift));

  return qVecShift;
}

TVector2 StZdcEpManager::ApplyZdcSmdShiftCorrWest(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  float Psi_ReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  float delta_Psi = 0.0;
  float Psi_Shift;

  for(Int_t i_shift = 0; i_shift < 20; ++i_shift) // Shift Order loop
  {
    int bin_Cos = p_mQWestCos[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Cos = p_mQWestCos[i_shift]->GetBinContent(bin_Cos);

    int bin_Sin = p_mQWestSin[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Sin = p_mQWestSin[i_shift]->GetBinContent(bin_Sin);

    delta_Psi += (2.0/((float)i_shift+1.0))*(-1.0*mean_Sin*TMath::Cos(((float)i_shift+1.0)*Psi_ReCenter)+mean_Cos*TMath::Sin(((float)i_shift+1.0)*Psi_ReCenter));
  }

  float Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  qVecShift.Set(TMath::Cos(Psi_Shift),TMath::Sin(Psi_Shift));

  return qVecShift;
}

float StZdcEpManager::AngleShift(float Psi_raw)
{
  float Psi_Corr = Psi_raw;
  if(Psi_raw > 1.0*TMath::Pi())
  {
    Psi_Corr = Psi_raw - TMath::TwoPi();
  }
  if(Psi_raw < -1.0*TMath::Pi())
  {
    Psi_Corr = Psi_raw + TMath::TwoPi();
  }

  return Psi_Corr;
}

//---------------------------------------------------------------------------------


void StZdcEpManager::readShiftCorrFull()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/ShiftParameterFull/file_%s_ZdcShiftParameterFull.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_ShiftPar = TFile::Open(InPutFile.c_str());

  for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
  {
    string ProName;

    ProName = Form("p_mQFullCos_%d",i_shift);
    p_mQFullCos[i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
    ProName = Form("p_mQFullSin_%d",i_shift);
    p_mQFullSin[i_shift] = (TProfile2D*)mFile_ShiftPar->Get(ProName.c_str());
  }
}

TVector2 StZdcEpManager::ApplyZdcSmdShiftCorrFull(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  float Psi_ReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  float delta_Psi = 0.0;
  float Psi_Shift;

  for(Int_t i_shift = 0; i_shift < 20; ++i_shift) // Shift Order loop
  {
    int bin_Cos = p_mQFullCos[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Cos = p_mQFullCos[i_shift]->GetBinContent(bin_Cos);

    int bin_Sin = p_mQFullSin[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    float mean_Sin = p_mQFullSin[i_shift]->GetBinContent(bin_Sin);

    delta_Psi += (2.0/((float)i_shift+1.0))*(-1.0*mean_Sin*TMath::Cos(((float)i_shift+1.0)*Psi_ReCenter)+mean_Cos*TMath::Sin(((float)i_shift+1.0)*Psi_ReCenter));
  }

  float Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  qVecShift.Set(TMath::Cos(Psi_Shift),TMath::Sin(Psi_Shift));

  return qVecShift;
}

TVector2 StZdcEpManager::getQFull(TVector2 QEast, TVector2 QWest)
{
  TVector2 qVector = QWest-QEast;
  TVector2 qVecShift = ApplyZdcSmdShiftCorrFull(qVector);

  return qVecShift;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readResolution()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/Resolution/file_%s_ZdcResolution.root",recoEP::mBeamEnergy[mEnergy].c_str());
  mFile_Resolution = TFile::Open(InPutFile.c_str());
  p_mResolution = (TProfile*)mFile_Resolution->Get("p_mResolution");
}

void StZdcEpManager::calResolution()
{
  TF1 *f_Res_Full = new TF1("f_Res_Full",Resolution_Full,0,10,0);
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    float val_res_full, err_res_full;
    float val_res_raw = p_mResolution->GetBinContent(i_cent+1);
    float err_res_raw = p_mResolution->GetBinError(i_cent+1);
    //    cout << "val_res_raw = " << val_res_raw << ", err_res_raw = " << err_res_raw << endl;
    if(val_res_raw <= 0)
    {
      val_res_full = -999.9;
      err_res_full = 1.0;
    }
    else
    {
      float val_res_sub = TMath::Sqrt(val_res_raw);
      float err_res_sub = err_res_raw/(2*val_res_sub);

      // calculate full event plane resolution
      float chi_sub = f_Res_Full->GetX(val_res_sub);
      float chi_full = chi_sub*TMath::Sqrt(2.0);
      val_res_full = f_Res_Full->Eval(chi_full);
      // error propagation
      float err_chi_sub = err_res_sub/f_Res_Full->Derivative(chi_sub);
      err_res_full = f_Res_Full->Derivative(chi_full)*err_chi_sub*TMath::Sqrt(2.0);
    }
    mResolution[i_cent] = val_res_full;
  }
}

float StZdcEpManager::getResolution(int Cent9)
{
  return mResolution[Cent9];
}
//---------------------------------------------------------------------------------
