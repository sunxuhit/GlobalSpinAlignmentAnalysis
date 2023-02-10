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
  clearZdcEpManager();
}

StZdcEpManager::~StZdcEpManager()
{
  /* */
}

void StZdcEpManager::clearZdcEpManager()
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
	mGainFactor[iEastWest][iVertHori][iSlat] = -999.9;
      }
    }
  }
  mCenterVertEast = -999.9;
  mCenterHoriEast = -999.9;
  mCenterVertWest = -999.9;
  mCenterHoriWest = -999.9;
}

void StZdcEpManager::initZdcEpManager(int cent9, int runIndex, int vzBin)
{
  mCent9    = cent9;
  mRunIndex = runIndex;
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
// Gain Correction
void StZdcEpManager::initZdcGainCorr()
{
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	std::string histName = Form("h_mZdcGainCorr%s%sSlat%d",str_mEastWest[iEastWest].c_str(),str_mVertHori[iVertHori].c_str(),iSlat);
	h_mZdcGainCorr[iEastWest][iVertHori][iSlat] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,5000,-4.5,4995.5);
      }
    }
  }
}

void StZdcEpManager::fillZdcGainCorr(int iEastWest, int iVertHori, int iSlat, double zdcsmd)
{
  h_mZdcGainCorr[iEastWest][iVertHori][iSlat]->Fill((double)mRunIndex,zdcsmd);
}

void StZdcEpManager::writeZdcGainCorr()
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

void StZdcEpManager::readZdcGainCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/GainCorrPar/file_%s_ZdcGainCorrFac.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mGainCorrPar = TFile::Open(inputFile.c_str());
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	std::string histName = Form("h_mZdcGainFactor%s%sSlat%d",str_mEastWest[iEastWest].c_str(),str_mVertHori[iVertHori].c_str(),iSlat);
	TH1F *h_GainFac = (TH1F*)file_mGainCorrPar->Get(histName.c_str());
	mGainFactor[iEastWest][iVertHori][iSlat] = h_GainFac->GetBinContent(1);
	// cout << "iEastWest = " << iEastWest << ", iVertHori = " << iVertHori << ", iSlat = " << iSlat << ", mGainFactor = " << mGainFactor[iEastWest][iVertHori][iSlat] << endl;
      }
    }
  }
  file_mGainCorrPar->Close();
}

void StZdcEpManager::setZdcSmdGainCorr(int eastwest, int verthori, int slat, const double zdcsmd)
{
  mZdcSmd[eastwest][verthori][slat] = (zdcsmd > 0. && mGainFactor[eastwest][verthori][slat] > 0.) ? zdcsmd/mGainFactor[eastwest][verthori][slat] : 0.;
  // cout << "input zdc = " << zdcsmd << ", mGainFactor = " << mGainFactor[eastwest][verthori][slat] << ", GainCorred = " << mZdcSmd[eastwest][verthori][slat] << endl;
}

double StZdcEpManager::getZdcSmdGainCorr(int eastwest, int verthori, int slat)
{
  return mZdcSmd[eastwest][verthori][slat];
}
//---------------------------------------------------------------------------------
// ReCenter Correction
void StZdcEpManager::initZdcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mZdcQReCenterVertEastVz%d",iVz);
    p_mZdcQReCenterVertEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    proName = Form("p_mZdcQReCenterHoriEastVz%d",iVz);
    p_mZdcQReCenterHoriEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);

    proName = Form("p_mZdcQReCenterVertWestVz%d",iVz);
    p_mZdcQReCenterVertWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    proName = Form("p_mZdcQReCenterHoriWestVz%d",iVz);
    p_mZdcQReCenterHoriWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
  }
}

void StZdcEpManager::fillZdcReCenterEast(TVector2 QVector)
{
  const double Qx = QVector.X();
  const double Qy = QVector.Y();

  p_mZdcQReCenterVertEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,Qx);
  p_mZdcQReCenterHoriEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,Qy);
}

void StZdcEpManager::fillZdcReCenterWest(TVector2 QVector)
{
  const double Qx = QVector.X();
  const double Qy = QVector.Y();

  p_mZdcQReCenterVertWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,Qx);
  p_mZdcQReCenterHoriWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,Qy);
}

void StZdcEpManager::writeZdcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mZdcQReCenterVertEast[iVz]->Write();
    p_mZdcQReCenterHoriEast[iVz]->Write();
    p_mZdcQReCenterVertWest[iVz]->Write();
    p_mZdcQReCenterHoriWest[iVz]->Write();
  }
}

void StZdcEpManager::readZdcReCenterCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_ZdcReCenterPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mReCenterPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mZdcQReCenterVertEastVz%d",iVz);
    p_mZdcQReCenterVertEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mZdcQReCenterHoriEastVz%d",iVz);
    p_mZdcQReCenterHoriEast[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mZdcQReCenterVertWestVz%d",iVz);
    p_mZdcQReCenterVertWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mZdcQReCenterHoriWestVz%d",iVz);
    p_mZdcQReCenterHoriWest[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
  }
}

void StZdcEpManager::setZdcSmdCenter()
{
  const int binVertEast = p_mZdcQReCenterVertEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterVertEast       = p_mZdcQReCenterVertEast[mVzBin]->GetBinContent(binVertEast);

  const int binHoriEast = p_mZdcQReCenterHoriEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterHoriEast       = p_mZdcQReCenterHoriEast[mVzBin]->GetBinContent(binHoriEast);

  const int binVertWest = p_mZdcQReCenterVertWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterVertWest       = -1.0*p_mZdcQReCenterVertWest[mVzBin]->GetBinContent(binVertWest);

  const int binHoriWest = p_mZdcQReCenterHoriWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterHoriWest       = p_mZdcQReCenterHoriWest[mVzBin]->GetBinContent(binHoriWest);
  // cout << "mCenterVertEast = " << mCenterVertEast << ", mCenterHoriEast = " << mCenterHoriEast << ", mCenterVertWest = " << mCenterVertWest << ", mCenterHoriWest = " << mCenterHoriWest << endl;
}

//---------------------------------------------------------------------------------
// Shift Correction East/West EP
void StZdcEpManager::initZdcShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQShiftCos%dEastVz%d",iShift,iVz);
      p_mZdcQShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
      proName = Form("p_mZdcQShiftSin%dEastVz%d",iShift,iVz);
      p_mZdcQShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);

      proName = Form("p_mZdcQShiftCos%dWestVz%d",iShift,iVz);
      p_mZdcQShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
      proName = Form("p_mZdcQShiftSin%dWestVz%d",iShift,iVz);
      p_mZdcQShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    }
  }
}

void StZdcEpManager::fillZdcShiftEast(TVector2 QVector)
{
  const double Psi = TMath::ATan2(QVector.Y(),QVector.X());
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double PsiCos = TMath::Cos(((double)iShift+1.0)*Psi);
    const double PsiSin = TMath::Sin(((double)iShift+1.0)*Psi);
    p_mZdcQShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiCos);
    p_mZdcQShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiSin);
  }
}

void StZdcEpManager::fillZdcShiftWest(TVector2 QVector)
{
  const double Psi = TMath::ATan2(QVector.Y(),QVector.X());
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double PsiCos = TMath::Cos(((double)iShift+1.0)*Psi);
    const double PsiSin = TMath::Sin(((double)iShift+1.0)*Psi);
    p_mZdcQShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiCos);
    p_mZdcQShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiSin);
  }
}

void StZdcEpManager::writeZdcShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mZdcQShiftCosEast[iVz][iShift]->Write();
      p_mZdcQShiftSinEast[iVz][iShift]->Write();
      p_mZdcQShiftCosWest[iVz][iShift]->Write();
      p_mZdcQShiftSinWest[iVz][iShift]->Write();
    }
  }
}

void StZdcEpManager::readZdcShiftCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftParameter/file_%s_ZdcShiftPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(inputFile.c_str());

  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQShiftCos%dEastVz%d",iShift,iVz);
      p_mZdcQShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mZdcQShiftSin%dEastVz%d",iShift,iVz);
      p_mZdcQShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mZdcQShiftCos%dWestVz%d",iShift,iVz);
      p_mZdcQShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mZdcQShiftSin%dWestVz%d",iShift,iVz);
      p_mZdcQShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

TVector2 StZdcEpManager::applyZdcSmdShiftCorrEast(TVector2 QVector)
{
  const double PsiReCenter = TMath::ATan2(QVector.Y(),QVector.X());

  double deltaPsi = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mZdcQShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mZdcQShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mZdcQShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mZdcQShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(((double)iShift+1.0)*PsiReCenter));
  }

  const double PsiShiftRaw = PsiReCenter + deltaPsi;
  const double PsiShift = angleShift(PsiShiftRaw);

  TVector2 QVecShift(0.0,0.0);
  QVecShift.Set(TMath::Cos(PsiShift),TMath::Sin(PsiShift));

  return QVecShift;
}

TVector2 StZdcEpManager::applyZdcSmdShiftCorrWest(TVector2 QVector)
{
  const double PsiReCenter = TMath::ATan2(QVector.Y(),QVector.X());

  double deltaPsi = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mZdcQShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mZdcQShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mZdcQShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mZdcQShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(((double)iShift+1.0)*PsiReCenter));
  }

  const double PsiShiftRaw = PsiReCenter + deltaPsi;
  const double PsiShift = angleShift(PsiShiftRaw);

  TVector2 QVecShift(0.0,0.0);
  QVecShift.Set(TMath::Cos(PsiShift),TMath::Sin(PsiShift));

  return QVecShift;
}

double StZdcEpManager::angleShift(double PsiRaw)
{
  double PsiCorr = PsiRaw;
  if(PsiRaw > 1.0*TMath::Pi())
  {
    PsiCorr = PsiRaw - TMath::TwoPi();
  }
  if(PsiRaw < -1.0*TMath::Pi())
  {
    PsiCorr = PsiRaw + TMath::TwoPi();
  }

  return PsiCorr;
}
//---------------------------------------------------------------------------------
// Shift Correction Full EP
void StZdcEpManager::initZdcShiftFull()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQShiftCos%dFullVz%d",iShift,iVz);
      p_mZdcQShiftCosFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
      proName = Form("p_mZdcQShiftSin%dFullVz%d",iShift,iVz);
      p_mZdcQShiftSinFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    }
  }
}

void StZdcEpManager::fillZdcShiftFull(TVector2 QVector)
{
  const double Psi = TMath::ATan2(QVector.Y(),QVector.X());
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double PsiCos = TMath::Cos(((double)iShift+1.0)*Psi);
    const double PsiSin = TMath::Sin(((double)iShift+1.0)*Psi);
    p_mZdcQShiftCosFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiCos);
    p_mZdcQShiftSinFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiSin);
  }
}

void StZdcEpManager::writeZdcShiftFull()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mZdcQShiftCosFull[iVz][iShift]->Write();
      p_mZdcQShiftSinFull[iVz][iShift]->Write();
    }
  }
}

void StZdcEpManager::readZdcShiftCorrFull()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftParameter/file_%s_ZdcShiftParFull.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(inputFile.c_str());

  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < 20; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQShiftCos%dFullVz%d",iShift,iVz);
      p_mZdcQShiftCosFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mZdcQShiftSin%dFullVz%d",iShift,iVz);
      p_mZdcQShiftSinFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

TVector2 StZdcEpManager::applyZdcSmdShiftCorrFull(TVector2 QVector)
{
  const double PsiReCenter = TMath::ATan2(QVector.Y(),QVector.X());

  double deltaPsi = 0.0;
  for(Int_t iShift = 0; iShift < 20; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mZdcQShiftCosFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mZdcQShiftCosFull[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mZdcQShiftSinFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mZdcQShiftSinFull[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(((double)iShift+1.0)*PsiReCenter));
  }

  const double PsiShiftRaw = PsiReCenter + deltaPsi;
  const double PsiShift = angleShift(PsiShiftRaw);

  TVector2 QVecShift(0.0,0.0);
  QVecShift.Set(TMath::Cos(PsiShift),TMath::Sin(PsiShift));

  return QVecShift;
}
//---------------------------------------------------------------------------------
// QVector
double StZdcEpManager::getPosition(int eastwest, int verthori, int slat, int mode)
{ //get position of each slat
  double zdcsmdVert[7] = {0.5,2,3.5,5,6.5,8,9.5};
  double zdcsmdHori[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

  if(mode > 1) // with beam center corrected
  {
    setZdcSmdCenter();
    if(mCenterVertEast < -100.0 || mCenterHoriEast < -100.0 || mCenterVertWest < -100.0 || mCenterHoriWest < -100.0) 
    {
      std::cout << "Forgot Re-Center!!!!" << std::endl;
      return 0;
    }
    if(eastwest == 0 && verthori == 0) return zdcsmdVert[slat]-mCenterVertEast;
    if(eastwest == 1 && verthori == 0) return mCenterVertWest-zdcsmdVert[slat];
    if(eastwest == 0 && verthori == 1) return zdcsmdHori[slat]/sqrt(2.)-mCenterHoriEast;
    if(eastwest == 1 && verthori == 1) return zdcsmdHori[slat]/sqrt(2.)-mCenterHoriWest;
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
  TVector2 QVector(0.0,0.0);
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

  if(qXwgt > 0.0 && qYwgt > 0.0) QVector.Set(qXsum/qXwgt,qYsum/qYwgt);
  if(mode > 2)  QVector = applyZdcSmdShiftCorrEast(QVector);

  return QVector;
}

TVector2 StZdcEpManager::getQWest(int mode)
{
  TVector2 QVector(0.0,0.0);
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

  if(qXwgt > 0.0 && qYwgt > 0.0) QVector.Set(qXsum/qXwgt,qYsum/qYwgt);
  if(mode > 2) QVector = applyZdcSmdShiftCorrWest(QVector);

  return QVector;
}

TVector2 StZdcEpManager::getQFull(TVector2 QEast, TVector2 QWest, int mode)
{
  // TVector2 QVector = QWest-QEast;
  TVector2 QVecShift = QWest - QEast;
  if(mode > 2) QVecShift = applyZdcSmdShiftCorrFull(QVecShift);

  return QVecShift;
}
//---------------------------------------------------------------------------------
// Full Event Plane Resolution
void StZdcEpManager::readZdcResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_%s_ZdcEpResolution.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mZdcEpResolution = (TProfile*)file_mResolution->Get("p_mZdcEpResolution");

  // for(int iCent = 0; iCent < 9; ++iCent) 
  // {
  //   mZdcResFullVal[iCent] = 0.0;
  //   mZdcResFullErr[iCent] = 0.0;
  // }

  calZdcResolution();
}

void StZdcEpManager::calZdcResolution()
{
  TF1 *f_ResFull = new TF1("f_ResFull",ResolutionFull,0,10,0);
  for(Int_t iCent = 0; iCent < 9; iCent++)
  {
    double valResFull, errResFull;
    double valResRaw = p_mZdcEpResolution->GetBinContent(iCent+1);
    double errResRaw = p_mZdcEpResolution->GetBinError(iCent+1);
    //    cout << "valResRaw = " << valResRaw << ", errResRaw = " << errResRaw << endl;
    if(valResRaw <= 0)
    {
      valResFull = -999.9;
      errResFull = 1.0;
    }
    else
    {
      double valResSub = TMath::Sqrt(valResRaw);
      double errResSub = errResRaw/(2*valResSub);

      // calculate full event plane resolution
      double chiSub = f_ResFull->GetX(valResSub);
      double chiFull = chiSub*TMath::Sqrt(2.0);
      valResFull = f_ResFull->Eval(chiFull);
      // error propagation
      double errChiSub = errResSub/f_ResFull->Derivative(chiSub);
      errResFull = f_ResFull->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
    }
    mZdcResFullVal[iCent] = valResFull;
    mZdcResFullErr[iCent] = errResFull;
  }
}

double StZdcEpManager::getZdcResFullVal(int cent9)
{
  return mZdcResFullVal[cent9];
}

double StZdcEpManager::getZdcResFullErr(int cent9)
{
  return mZdcResFullErr[cent9];
}
//---------------------------------------------------------------------------------
// raw EP
void StZdcEpManager::initZdcRawEP()
{
  for(int iCent = 0; iCent < 9; ++iCent)
  {
    std::string histName = Form("h_mZdcRawEpEastCent%d",iCent);
    h_mZdcRawEpEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5);
    histName = Form("h_mZdcRawEpWestCent%d",iCent);
    h_mZdcRawEpWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5);
    histName = Form("h_mZdcRawEpFullCent%d",iCent);
    h_mZdcRawEpFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5);
  }
}

void StZdcEpManager::fillZdcRawEP(TVector2 QEast, TVector2 QWest, TVector2 QFull)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcRawEpEast[mCent9]->Fill(PsiEast,mRunIndex);
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcRawEpWest[mCent9]->Fill(PsiWest,mRunIndex);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcRawEpFull[mCent9]->Fill(PsiFull,mRunIndex);
}

void StZdcEpManager::writeZdcRawEP()
{
  for(int iCent = 0; iCent < 9; ++iCent)
  {
    h_mZdcRawEpEast[iCent]->Write();
    h_mZdcRawEpWest[iCent]->Write();
    h_mZdcRawEpFull[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
