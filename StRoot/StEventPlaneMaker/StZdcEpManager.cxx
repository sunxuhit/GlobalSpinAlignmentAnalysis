#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TF1.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StEventPlaneMaker/StZdcEpManager.h"

double funcZdcEpResFull(double *x_val, double *par)
{
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  double resFull = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return resFull;
}

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
void StZdcEpManager::initZdcGain()
{
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	std::string histName = Form("h_mZdcGainCorr%s%sSlat%d",str_mEastWest[iEastWest].c_str(),str_mVertHori[iVertHori].c_str(),iSlat);
	h_mZdcGainCorr[iEastWest][iVertHori][iSlat] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,1500,-4.5,1495.5);
      }
    }
  }
}

void StZdcEpManager::fillZdcGain(int iEastWest, int iVertHori, int iSlat, double zdcsmd)
{
  h_mZdcGainCorr[iEastWest][iVertHori][iSlat]->Fill((double)mRunIndex,zdcsmd);
}

void StZdcEpManager::writeZdcGain()
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

void StZdcEpManager::readZdcGain()
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
    p_mZdcQReCenterVertEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mZdcQReCenterHoriEastVz%d",iVz);
    p_mZdcQReCenterHoriEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mZdcQReCenterVertWestVz%d",iVz);
    p_mZdcQReCenterVertWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mZdcQReCenterHoriWestVz%d",iVz);
    p_mZdcQReCenterHoriWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
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

void StZdcEpManager::readZdcReCenter()
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
      p_mZdcQShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mZdcQShiftSin%dEastVz%d",iShift,iVz);
      p_mZdcQShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mZdcQShiftCos%dWestVz%d",iShift,iVz);
      p_mZdcQShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mZdcQShiftSin%dWestVz%d",iShift,iVz);
      p_mZdcQShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
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

void StZdcEpManager::readZdcShift()
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
  const double PsiShift = transPsi1(PsiShiftRaw);

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
  const double PsiShift = transPsi1(PsiShiftRaw);

  TVector2 QVecShift(0.0,0.0);
  QVecShift.Set(TMath::Cos(PsiShift),TMath::Sin(PsiShift));

  return QVecShift;
}

double StZdcEpManager::transPsi1(double Psi)
{
  double PsiCorr = Psi;
  if(Psi >  1.0*TMath::Pi()) PsiCorr = Psi - TMath::TwoPi();
  if(Psi < -1.0*TMath::Pi()) PsiCorr = Psi + TMath::TwoPi();

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
      p_mZdcQShiftCosFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mZdcQShiftSin%dFullVz%d",iShift,iVz);
      p_mZdcQShiftSinFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
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

void StZdcEpManager::readZdcShiftFull()
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
  const double PsiShift = transPsi1(PsiShiftRaw);

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

TVector2 StZdcEpManager::getQ1VecEast(int mode) 
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

TVector2 StZdcEpManager::getQ1VecWest(int mode)
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

TVector2 StZdcEpManager::getQ1VecFull(TVector2 QEast, TVector2 QWest, int mode)
{
  // TVector2 QVector = QWest-QEast;
  TVector2 QVecShift = QWest - QEast;
  if(mode > 2) QVecShift = applyZdcSmdShiftCorrFull(QVecShift);

  return QVecShift;
}
//---------------------------------------------------------------------------------
// Full Event Plane Resolution
void StZdcEpManager::initZdcResolution()
{
  p_mZdcSubEpRes = new TProfile("p_mZdcSubEpRes","p_mZdcSubEpRes",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
}

void StZdcEpManager::fillZdcResolution(TVector2 QEast, TVector2 QWest)
{
  double PsiEast = TMath::ATan2(-1.0*QEast.Y(),-1.0*QEast.X()); // flip the sign for Q1VecEast
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X());
  double res1Sub = TMath::Cos(PsiWest-PsiEast);
  p_mZdcSubEpRes->Fill((double)mCent9,res1Sub);
}

void StZdcEpManager::writeZdcResolution()
{
  p_mZdcSubEpRes->Write();
}

void StZdcEpManager::readZdcResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_%s_ZdcEpResolution.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mZdcSubEpRes = (TProfile*)file_mResolution->Get("p_mZdcSubEpRes");

  for(int iCent = 0; iCent < 9; ++iCent) 
  {
    mZdcFullEpResVal[iCent] = 0.0;
    mZdcFullEpResErr[iCent] = 0.0;
  }

  TF1 *f_ZdcEpResFull = new TF1("f_ZdcEpResFull",funcZdcEpResFull,0,10,0);
  for(Int_t iCent = 0; iCent < mNumCentrality; iCent++)
  {
    double valResSub  = -999.9;
    double errResSub  = 1.0;
    double valResFull = -999.9;
    double errResFull = 1.0;
    double valResRaw = p_mZdcSubEpRes->GetBinContent(iCent+1);
    double errResRaw = p_mZdcSubEpRes->GetBinError(iCent+1);
    //    cout << "valResRaw = " << valResRaw << ", errResRaw = " << errResRaw << endl;
    if(valResRaw > 0)
    {
      double valResSub = TMath::Sqrt(valResRaw);
      double errResSub = errResRaw/(2*valResSub);

      // calculate full event plane resolution
      double chiSub = f_ZdcEpResFull->GetX(valResSub);
      double chiFull = chiSub*TMath::Sqrt(2.0);
      valResFull = f_ZdcEpResFull->Eval(chiFull);
      // error propagation
      double errChiSub = errResSub/f_ZdcEpResFull->Derivative(chiSub);
      errResFull = f_ZdcEpResFull->Derivative(chiFull)*errChiSub*TMath::Sqrt(2.0);
    }
    mZdcSubEpResVal[iCent]  = valResSub;
    mZdcSubEpResErr[iCent]  = errResSub;
    mZdcFullEpResVal[iCent] = valResFull;
    mZdcFullEpResErr[iCent] = errResFull;
  }

  file_mResolution->Close();
}

double StZdcEpManager::getZdcSubEpResVal(int cent9)
{
  return mZdcSubEpResVal[cent9];
}

double StZdcEpManager::getZdcSubEpResErr(int cent9)
{
  return mZdcSubEpResErr[cent9];
}

double StZdcEpManager::getZdcFullEpResVal(int cent9)
{
  return mZdcFullEpResVal[cent9];
}

double StZdcEpManager::getZdcFullEpResErr(int cent9)
{
  return mZdcFullEpResErr[cent9];
}
//---------------------------------------------------------------------------------
// raw EP
void StZdcEpManager::initZdcSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEpRawEastCent%d",iCent);
    h_mZdcEpRawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEpRawWestCent%d",iCent);
    h_mZdcEpRawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEpRawFullCent%d",iCent);
    h_mZdcEpRawFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEpRawCorrCent%d",iCent);
    h_mZdcEpRawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StZdcEpManager::fillZdcSubEpRaw(TVector2 QEast, TVector2 QWest, TVector2 QFull)
{
  double PsiEast = TMath::ATan2(-1.0*QEast.Y(),-1.0*QEast.X()); h_mZdcEpRawEast[mCent9]->Fill(mRunIndex,PsiEast); // flip the sign for Q1VecEast
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcEpRawWest[mCent9]->Fill(mRunIndex,PsiWest);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcEpRawFull[mCent9]->Fill(mRunIndex,PsiFull);
  h_mZdcEpRawCorr[mCent9]->Fill(PsiEast,PsiWest);
}

void StZdcEpManager::writeZdcSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEpRawEast[iCent]->Write();
    h_mZdcEpRawWest[iCent]->Write();
    h_mZdcEpRawFull[iCent]->Write();
    h_mZdcEpRawCorr[iCent]->Write();
  }
}

// recenter EP
void StZdcEpManager::initZdcSubEpReCenter()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEpReCenterEastCent%d",iCent);
    h_mZdcEpReCenterEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEpReCenterWestCent%d",iCent);
    h_mZdcEpReCenterWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEpReCenterFullCent%d",iCent);
    h_mZdcEpReCenterFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEpReCenterCorrCent%d",iCent);
    h_mZdcEpReCenterCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StZdcEpManager::fillZdcSubEpReCenter(TVector2 QEast, TVector2 QWest, TVector2 QFull)
{
  double PsiEast = TMath::ATan2(-1.0*QEast.Y(),-1.0*QEast.X()); h_mZdcEpReCenterEast[mCent9]->Fill(mRunIndex,PsiEast); // flip the sign for Q1VecEast
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcEpReCenterWest[mCent9]->Fill(mRunIndex,PsiWest);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcEpReCenterFull[mCent9]->Fill(mRunIndex,PsiFull);
  h_mZdcEpReCenterCorr[mCent9]->Fill(PsiEast,PsiWest);
}

void StZdcEpManager::writeZdcSubEpReCenter()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEpReCenterEast[iCent]->Write();
    h_mZdcEpReCenterWest[iCent]->Write();
    h_mZdcEpReCenterFull[iCent]->Write();
    h_mZdcEpReCenterCorr[iCent]->Write();
  }
}

// shift EP
void StZdcEpManager::initZdcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEpShiftEastCent%d",iCent);
    h_mZdcEpShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEpShiftWestCent%d",iCent);
    h_mZdcEpShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEpShiftFullCent%d",iCent);
    h_mZdcEpShiftFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEpShiftCorrCent%d",iCent);
    h_mZdcEpShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),360,-1.0*TMath::Pi(),TMath::Pi(),360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StZdcEpManager::fillZdcSubEpShift(TVector2 QEast, TVector2 QWest, TVector2 QFull)
{
  double PsiEast = TMath::ATan2(-1.0*QEast.Y(),-1.0*QEast.X()); h_mZdcEpShiftEast[mCent9]->Fill(mRunIndex,PsiEast); // flip the sign for Q1VecEast
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcEpShiftWest[mCent9]->Fill(mRunIndex,PsiWest);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcEpShiftFull[mCent9]->Fill(mRunIndex,PsiFull);
  h_mZdcEpShiftCorr[mCent9]->Fill(PsiEast,PsiWest);
}

void StZdcEpManager::writeZdcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEpShiftEast[iCent]->Write();
    h_mZdcEpShiftWest[iCent]->Write();
    h_mZdcEpShiftFull[iCent]->Write();
    h_mZdcEpShiftCorr[iCent]->Write();
  }
}

// shift Full EP
void StZdcEpManager::initZdcFullEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEpShiftFullCorrCent%d",iCent);
    h_mZdcEpShiftFullCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StZdcEpManager::fillZdcFullEpShift(TVector2 QFull)
{
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); 
  h_mZdcEpShiftFullCorr[mCent9]->Fill(mRunIndex,PsiFull);
}

void StZdcEpManager::writeZdcFullEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEpShiftFullCorr[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
