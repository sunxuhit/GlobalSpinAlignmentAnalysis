#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TF1.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "Utility/include/StSpinAlignmentFunctions.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"
// #include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

ClassImp(StZdcEpManager)

//---------------------------------------------------------------------------------

StZdcEpManager::StZdcEpManager(int beamType) : mType(beamType)
{
  // mType = energy;
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
  mVzBin = -1;
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	mZdcSmd[iEastWest][iVertHori][iSlat] = 0.0;
	mGainCorrFactor[iEastWest][iVertHori][iSlat] = -999.9;
      }
    }
  }
  mCenterEastVertical   = -999.9;
  mCenterEastHorizontal = -999.9;
  mCenterWestVertical   = -999.9;
  mCenterWestHorizontal = -999.9;
  for(int i_cent = 0; i_cent < 9; ++i_cent) mResolution[i_cent] = 0.0;
}

void StZdcEpManager::initZdcEp(int Cent9, int RunIndex, int vzBin)
{
  mCent9    = Cent9;
  mRunIndex = RunIndex;
  mVzBin    = vzBin;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::setZdcSmd(int eastwest, int verthori, int slat, const double zdcsmd) 
{
  mZdcSmd[eastwest][verthori][slat] = (zdcsmd > 0.) ? zdcsmd : 0.;
}

double StZdcEpManager::getZdcSmd(int eastwest, int verthori, int slat)
{
  return mZdcSmd[eastwest][verthori][slat];
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readGainCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/GainCorrPar/file_%s_ZdcGainCorrFac.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mGainCorrPar = TFile::Open(inputFile.c_str());
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	std::string histName = Form("h_mZdcGainCorrFactor%s%sSlat%d",str_mEastWest[iEastWest].c_str(),str_mVertHori[iVertHori].c_str(),iSlat);
	TH1F *h_GainCorrFac = (TH1F*)file_mGainCorrPar->Get(histName.c_str());
	mGainCorrFactor[iEastWest][iVertHori][iSlat] = h_GainCorrFac->GetBinContent(1);
	// cout << "iEastWest = " << iEastWest << ", iVertHori = " << iVertHori << ", iSlat = " << iSlat << ", mGainCorrFactor = " << mGainCorrFactor[iEastWest][iVertHori][iSlat] << endl;
      }
    }
  }
  file_mGainCorrPar->Close();
}

void StZdcEpManager::setZdcSmdGainCorr(int eastwest, int verthori, int slat, const double zdcsmd)
{
  mZdcSmd[eastwest][verthori][slat] = (zdcsmd > 0. && mGainCorrFactor[eastwest][verthori][slat] > 0.) ? zdcsmd/mGainCorrFactor[eastwest][verthori][slat] : 0.;
  // cout << "input zdc = " << zdcsmd << ", mGainCorrFactor = " << mGainCorrFactor[eastwest][verthori][slat] << ", GainCorred = " << mZdcSmd[eastwest][verthori][slat] << endl;
}

double StZdcEpManager::getZdcSmdGainCorr(int eastwest, int verthori, int slat)
{
  return mZdcSmd[eastwest][verthori][slat];
}

double StZdcEpManager::getPosition(int eastwest, int verthori, int slat, int mode)
{ //get position of each slat
  double zdcsmdVert[7] = {0.5,2,3.5,5,6.5,8,9.5};
  double zdcsmdHori[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

  if(mode > 1) // with beam center corrected
  {
    setZdcSmdCenter();
    if(mCenterEastVertical < -100.0 || mCenterEastHorizontal < -100.0 || mCenterWestVertical < -100.0 || mCenterWestHorizontal < -100.0) 
    {
      std::cout << "Forgot Re-Center!!!!" << std::endl;
      return 0;
    }
    if(eastwest == 0 && verthori == 0) return zdcsmdVert[slat]-mCenterEastVertical;
    if(eastwest == 1 && verthori == 0) return mCenterWestVertical-zdcsmdVert[slat];
    if(eastwest == 0 && verthori == 1) return zdcsmdHori[slat]/sqrt(2.)-mCenterEastHorizontal;
    if(eastwest == 1 && verthori == 1) return zdcsmdHori[slat]/sqrt(2.)-mCenterWestHorizontal;
  }
  else // raw beam center returned
  {
    if(eastwest == 0 && verthori == 0) return zdcsmdVert[slat];
    if(eastwest == 1 && verthori == 0) return -zdcsmdVert[slat];
    if(eastwest == 0 && verthori == 1) return zdcsmdHori[slat]/sqrt(2.);
    if(eastwest == 1 && verthori == 1) return zdcsmdHori[slat]/sqrt(2.);
  }

  return 0;
}

TVector2 StZdcEpManager::getQEast(int mode) 
{
  TVector2 qVector(0.0,0.0);
  double qXsum = 0.; double qYsum = 0.;
  double qXwgt = 0.; double qYwgt = 0.;

  for(int iVert = 0; iVert < 7; ++iVert) // vertical
  {
    qXsum += getPosition(0,0,iVert,mode)*getZdcSmdGainCorr(0,0,iVert);
    qXwgt += getZdcSmdGainCorr(0,0,iVert);
  }
  for(int iHori = 0; iHori < 8; ++iHori) // horizontal
  {
    qYsum += getPosition(0,1,iHori,mode)*getZdcSmdGainCorr(0,1,iHori);
    qYwgt += getZdcSmdGainCorr(0,1,iHori);
  }

  if(qXwgt > 0.0 && qYwgt > 0.0) qVector.Set(qXsum/qXwgt,qYsum/qYwgt);
  // if(mode > 2)  qVector = ApplyZdcSmdShiftCorrEast(qVector);

  return qVector;
}

TVector2 StZdcEpManager::getQWest(int mode)
{
  TVector2 qVector(0.0,0.0);
  double qXsum = 0.; double qYsum = 0.;
  double qXwgt = 0.; double qYwgt = 0.;

  for(int iVert = 0; iVert < 7; ++iVert) // vertical
  {
    qXsum += getPosition(1,0,iVert,mode)*getZdcSmdGainCorr(1,0,iVert);
    qXwgt += getZdcSmdGainCorr(1,0,iVert);
  }
  for(int iHori = 0; iHori < 8; ++iHori) // horizontal
  {
    qYsum += getPosition(1,1,iHori,mode)*getZdcSmdGainCorr(1,1,iHori);
    qYwgt += getZdcSmdGainCorr(1,1,iHori);
  }

  if(qXwgt > 0.0 && qYwgt > 0.0) qVector.Set(qXsum/qXwgt,qYsum/qYwgt);
  // if(mode > 2) qVector = ApplyZdcSmdShiftCorrWest(qVector);

  return qVector;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readReCenterCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_ZdcReCenterPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mReCenterPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string ProName = Form("p_mZdcQEastVertVz%d",iVz);
    p_mZdcQEastVertical[iVz] = (TProfile2D*)file_mReCenterPar->Get(ProName.c_str());
    ProName = Form("p_mZdcQEastHoriVz%d",iVz);
    p_mZdcQEastHorizontal[iVz] = (TProfile2D*)file_mReCenterPar->Get(ProName.c_str());

    ProName = Form("p_mZdcQWestVertVz%d",iVz);
    p_mZdcQWestVertical[iVz] = (TProfile2D*)file_mReCenterPar->Get(ProName.c_str());
    ProName = Form("p_mZdcQWestHoriVz%d",iVz);
    p_mZdcQWestHorizontal[iVz] = (TProfile2D*)file_mReCenterPar->Get(ProName.c_str());
  }
}

void StZdcEpManager::setZdcSmdCenter()
{
  const int binEastVertical = p_mZdcQEastVertical[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterEastVertical       = p_mZdcQEastVertical[mVzBin]->GetBinContent(binEastVertical);

  const int binEastHorizontal = p_mZdcQEastHorizontal[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterEastHorizontal       = p_mZdcQEastHorizontal[mVzBin]->GetBinContent(binEastHorizontal);

  const int binWestVertical = p_mZdcQWestVertical[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterWestVertical       = -1.0*p_mZdcQWestVertical[mVzBin]->GetBinContent(binWestVertical);

  const int binWestHorizontal = p_mZdcQWestHorizontal[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterWestHorizontal       = p_mZdcQWestHorizontal[mVzBin]->GetBinContent(binWestHorizontal);
  // cout << "mCenterEastVertical = " << mCenterEastVertical << ", mCenterEastHorizontal = " << mCenterEastHorizontal << ", mCenterWestVertical = " << mCenterWestVertical << ", mCenterWestHorizontal = " << mCenterWestHorizontal << endl;
}

#if 0
//---------------------------------------------------------------------------------

void StZdcEpManager::readShiftCorr()
{
  string InPutFile = Form("StRoot/StEventPlaneUtility/ShiftParameter/file_%s_ZdcShiftParameter.root",recoEP::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(InPutFile.c_str());

  for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
  {
    string ProName;

    ProName = Form("p_mQEastCos_%d",i_shift);
    p_mQEastCos[i_shift] = (TProfile2D*)file_mShiftPar->Get(ProName.c_str());
    ProName = Form("p_mQEastSin_%d",i_shift);
    p_mQEastSin[i_shift] = (TProfile2D*)file_mShiftPar->Get(ProName.c_str());

    ProName = Form("p_mQWestCos_%d",i_shift);
    p_mQWestCos[i_shift] = (TProfile2D*)file_mShiftPar->Get(ProName.c_str());
    ProName = Form("p_mQWestSin_%d",i_shift);
    p_mQWestSin[i_shift] = (TProfile2D*)file_mShiftPar->Get(ProName.c_str());
  }
}

TVector2 StZdcEpManager::ApplyZdcSmdShiftCorrEast(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  double Psi_ReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  double delta_Psi = 0.0;
  double Psi_Shift;

  for(Int_t i_shift = 0; i_shift < 20; ++i_shift) // Shift Order loop
  {
    int bin_Cos = p_mQEastCos[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    double mean_Cos = p_mQEastCos[i_shift]->GetBinContent(bin_Cos);

    int bin_Sin = p_mQEastSin[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    double mean_Sin = p_mQEastSin[i_shift]->GetBinContent(bin_Sin);

    delta_Psi += (2.0/((double)i_shift+1.0))*(-1.0*mean_Sin*TMath::Cos(((double)i_shift+1.0)*Psi_ReCenter)+mean_Cos*TMath::Sin(((double)i_shift+1.0)*Psi_ReCenter));
  }

  double Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  qVecShift.Set(TMath::Cos(Psi_Shift),TMath::Sin(Psi_Shift));

  return qVecShift;
}

TVector2 StZdcEpManager::ApplyZdcSmdShiftCorrWest(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  double Psi_ReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  double delta_Psi = 0.0;
  double Psi_Shift;

  for(Int_t i_shift = 0; i_shift < 20; ++i_shift) // Shift Order loop
  {
    int bin_Cos = p_mQWestCos[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    double mean_Cos = p_mQWestCos[i_shift]->GetBinContent(bin_Cos);

    int bin_Sin = p_mQWestSin[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    double mean_Sin = p_mQWestSin[i_shift]->GetBinContent(bin_Sin);

    delta_Psi += (2.0/((double)i_shift+1.0))*(-1.0*mean_Sin*TMath::Cos(((double)i_shift+1.0)*Psi_ReCenter)+mean_Cos*TMath::Sin(((double)i_shift+1.0)*Psi_ReCenter));
  }

  double Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  qVecShift.Set(TMath::Cos(Psi_Shift),TMath::Sin(Psi_Shift));

  return qVecShift;
}

double StZdcEpManager::AngleShift(double Psi_raw)
{
  double Psi_Corr = Psi_raw;
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
  string InPutFile = Form("StRoot/StEventPlaneUtility/ShiftParameterFull/file_%s_ZdcShiftParameterFull.root",recoEP::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(InPutFile.c_str());

  for(int i_shift = 0; i_shift < 20; ++i_shift) // Shift Order
  {
    string ProName;

    ProName = Form("p_mQFullCos_%d",i_shift);
    p_mQFullCos[i_shift] = (TProfile2D*)file_mShiftPar->Get(ProName.c_str());
    ProName = Form("p_mQFullSin_%d",i_shift);
    p_mQFullSin[i_shift] = (TProfile2D*)file_mShiftPar->Get(ProName.c_str());
  }
}

TVector2 StZdcEpManager::ApplyZdcSmdShiftCorrFull(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  double Psi_ReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  double delta_Psi = 0.0;
  double Psi_Shift;

  for(Int_t i_shift = 0; i_shift < 20; ++i_shift) // Shift Order loop
  {
    int bin_Cos = p_mQFullCos[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    double mean_Cos = p_mQFullCos[i_shift]->GetBinContent(bin_Cos);

    int bin_Sin = p_mQFullSin[i_shift]->FindBin((double)mRunIndex,(double)mCent9);
    double mean_Sin = p_mQFullSin[i_shift]->GetBinContent(bin_Sin);

    delta_Psi += (2.0/((double)i_shift+1.0))*(-1.0*mean_Sin*TMath::Cos(((double)i_shift+1.0)*Psi_ReCenter)+mean_Cos*TMath::Sin(((double)i_shift+1.0)*Psi_ReCenter));
  }

  double Psi_Shift_raw = Psi_ReCenter + delta_Psi;
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
  string InPutFile = Form("StRoot/StEventPlaneUtility/Resolution/file_%s_ZdcResolution.root",recoEP::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(InPutFile.c_str());
  p_mResolution = (TProfile*)file_mResolution->Get("p_mResolution");
}

void StZdcEpManager::calResolution()
{
  TF1 *f_Res_Full = new TF1("f_Res_Full",Resolution_FullZdc,0,10,0);
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    double val_res_full, err_res_full;
    double val_res_raw = p_mResolution->GetBinContent(i_cent+1);
    double err_res_raw = p_mResolution->GetBinError(i_cent+1);
    //    cout << "val_res_raw = " << val_res_raw << ", err_res_raw = " << err_res_raw << endl;
    if(val_res_raw <= 0)
    {
      val_res_full = -999.9;
      err_res_full = 1.0;
    }
    else
    {
      double val_res_sub = TMath::Sqrt(val_res_raw);
      double err_res_sub = err_res_raw/(2*val_res_sub);

      // calculate full event plane resolution
      double chi_sub = f_Res_Full->GetX(val_res_sub);
      double chi_full = chi_sub*TMath::Sqrt(2.0);
      val_res_full = f_Res_Full->Eval(chi_full);
      // error propagation
      double err_chi_sub = err_res_sub/f_Res_Full->Derivative(chi_sub);
      err_res_full = f_Res_Full->Derivative(chi_full)*err_chi_sub*TMath::Sqrt(2.0);
    }
    mResolution[i_cent] = val_res_full;
  }
}

double StZdcEpManager::getResolution(int Cent9)
{
  return mResolution[Cent9];
}
#endif
//---------------------------------------------------------------------------------
