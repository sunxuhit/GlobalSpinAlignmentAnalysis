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
  for(int iEastWest = 0; iEastWest < 2; ++iEastWest)
  {
    for(int iVertHori = 0; iVertHori < 2; ++iVertHori)
    {
      for(int iSlat = 0; iSlat < 8; ++iSlat)
      {
	mGainFactor[iEastWest][iVertHori][iSlat] = -999.9;
      }
    }
  }
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
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/GainCorrPar/file_ZdcGainCorrFac_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
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
void StZdcEpManager::initZdcReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mZdcQ1ReCtrVertEastVz%d",iVz);
    p_mZdcQ1ReCtrVertEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mZdcQ1ReCtrHoriEastVz%d",iVz);
    p_mZdcQ1ReCtrHoriEast[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

    proName = Form("p_mZdcQ1ReCtrVertWestVz%d",iVz);
    p_mZdcQ1ReCtrVertWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
    proName = Form("p_mZdcQ1ReCtrHoriWestVz%d",iVz);
    p_mZdcQ1ReCtrHoriWest[iVz] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  }
}

void StZdcEpManager::fillZdcReCtrEast(TVector2 QVector)
{
  const double Qx = QVector.X();
  const double Qy = QVector.Y();

  p_mZdcQ1ReCtrVertEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,Qx);
  p_mZdcQ1ReCtrHoriEast[mVzBin]->Fill((double)mRunIndex,(double)mCent9,Qy);
}

void StZdcEpManager::fillZdcReCtrWest(TVector2 QVector)
{
  const double Qx = QVector.X();
  const double Qy = QVector.Y();

  p_mZdcQ1ReCtrVertWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,Qx);
  p_mZdcQ1ReCtrHoriWest[mVzBin]->Fill((double)mRunIndex,(double)mCent9,Qy);
}

void StZdcEpManager::writeZdcReCtr()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mZdcQ1ReCtrVertEast[iVz]->Write();
    p_mZdcQ1ReCtrHoriEast[iVz]->Write();
    p_mZdcQ1ReCtrVertWest[iVz]->Write();
    p_mZdcQ1ReCtrHoriWest[iVz]->Write();
  }
}

void StZdcEpManager::readZdcReCtr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_ZdcReCenterPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mReCtrPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mZdcQ1ReCtrVertEastVz%d",iVz);
    p_mZdcQ1ReCtrVertEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mZdcQ1ReCtrHoriEastVz%d",iVz);
    p_mZdcQ1ReCtrHoriEast[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());

    proName = Form("p_mZdcQ1ReCtrVertWestVz%d",iVz);
    p_mZdcQ1ReCtrVertWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
    proName = Form("p_mZdcQ1ReCtrHoriWestVz%d",iVz);
    p_mZdcQ1ReCtrHoriWest[iVz] = (TProfile2D*)file_mReCtrPar->Get(proName.c_str());
  }
}

void StZdcEpManager::setZdcSmdCenter()
{
  const int binVertEast = p_mZdcQ1ReCtrVertEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterVertEast       = p_mZdcQ1ReCtrVertEast[mVzBin]->GetBinContent(binVertEast);

  const int binHoriEast = p_mZdcQ1ReCtrHoriEast[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterHoriEast       = p_mZdcQ1ReCtrHoriEast[mVzBin]->GetBinContent(binHoriEast);

  const int binVertWest = p_mZdcQ1ReCtrVertWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterVertWest       = -1.0*p_mZdcQ1ReCtrVertWest[mVzBin]->GetBinContent(binVertWest);

  const int binHoriWest = p_mZdcQ1ReCtrHoriWest[mVzBin]->FindBin((double)mRunIndex,(double)mCent9);
  mCenterHoriWest       = p_mZdcQ1ReCtrHoriWest[mVzBin]->GetBinContent(binHoriWest);
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
      std::string proName = Form("p_mZdcQ1ShiftCos%dEastVz%d",iShift,iVz);
      p_mZdcQ1ShiftCosEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mZdcQ1ShiftSin%dEastVz%d",iShift,iVz);
      p_mZdcQ1ShiftSinEast[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);

      proName = Form("p_mZdcQ1ShiftCos%dWestVz%d",iShift,iVz);
      p_mZdcQ1ShiftCosWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mZdcQ1ShiftSin%dWestVz%d",iShift,iVz);
      p_mZdcQ1ShiftSinWest[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
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
    p_mZdcQ1ShiftCosEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiCos);
    p_mZdcQ1ShiftSinEast[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiSin);
  }
}

void StZdcEpManager::fillZdcShiftWest(TVector2 QVector)
{
  const double Psi = TMath::ATan2(QVector.Y(),QVector.X());
  for(int iShift = 0; iShift < mNumShiftCorr; ++iShift)
  {
    const double PsiCos = TMath::Cos(((double)iShift+1.0)*Psi);
    const double PsiSin = TMath::Sin(((double)iShift+1.0)*Psi);
    p_mZdcQ1ShiftCosWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiCos);
    p_mZdcQ1ShiftSinWest[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiSin);
  }
}

void StZdcEpManager::writeZdcShift()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mZdcQ1ShiftCosEast[iVz][iShift]->Write();
      p_mZdcQ1ShiftSinEast[iVz][iShift]->Write();
      p_mZdcQ1ShiftCosWest[iVz][iShift]->Write();
      p_mZdcQ1ShiftSinWest[iVz][iShift]->Write();
    }
  }
}

void StZdcEpManager::readZdcShift()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftPar/file_ZdcShiftPar_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(inputFile.c_str());

  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQ1ShiftCos%dEastVz%d",iShift,iVz);
      p_mZdcQ1ShiftCosEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mZdcQ1ShiftSin%dEastVz%d",iShift,iVz);
      p_mZdcQ1ShiftSinEast[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mZdcQ1ShiftCos%dWestVz%d",iShift,iVz);
      p_mZdcQ1ShiftCosWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mZdcQ1ShiftSin%dWestVz%d",iShift,iVz);
      p_mZdcQ1ShiftSinWest[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

TVector2 StZdcEpManager::applyZdcSmdShiftCorrEast(TVector2 QVector)
{
  const double PsiReCenter = TMath::ATan2(QVector.Y(),QVector.X());

  double deltaPsi = 0.0;
  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mZdcQ1ShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mZdcQ1ShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mZdcQ1ShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mZdcQ1ShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

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
    const int binCos     = p_mZdcQ1ShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mZdcQ1ShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mZdcQ1ShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mZdcQ1ShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

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
      std::string proName = Form("p_mZdcQ1ShiftCos%dFullVz%d",iShift,iVz);
      p_mZdcQ1ShiftCosFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
      proName = Form("p_mZdcQ1ShiftSin%dFullVz%d",iShift,iVz);
      p_mZdcQ1ShiftSinFull[iVz][iShift] = new TProfile2D(proName.c_str(),proName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,mNumCentrality,-0.5,(double)mNumCentrality-0.5);
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
    p_mZdcQ1ShiftCosFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiCos);
    p_mZdcQ1ShiftSinFull[mVzBin][iShift]->Fill((double)mRunIndex,(double)mCent9,PsiSin);
  }
}

void StZdcEpManager::writeZdcShiftFull()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      p_mZdcQ1ShiftCosFull[iVz][iShift]->Write();
      p_mZdcQ1ShiftSinFull[iVz][iShift]->Write();
    }
  }
}

void StZdcEpManager::readZdcShiftFull()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftPar/file_ZdcShiftParFull_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(inputFile.c_str());

  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < 20; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQ1ShiftCos%dFullVz%d",iShift,iVz);
      p_mZdcQ1ShiftCosFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mZdcQ1ShiftSin%dFullVz%d",iShift,iVz);
      p_mZdcQ1ShiftSinFull[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

TVector2 StZdcEpManager::applyZdcSmdShiftCorrFull(TVector2 QVector)
{
  const double PsiReCenter = TMath::ATan2(QVector.Y(),QVector.X());

  double deltaPsi = 0.0;
  for(Int_t iShift = 0; iShift < 20; ++iShift) // Shift Order loop
  {
    const int binCos     = p_mZdcQ1ShiftCosFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanCos = p_mZdcQ1ShiftCosFull[mVzBin][iShift]->GetBinContent(binCos);

    const int binSin     = p_mZdcQ1ShiftSinFull[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    const double meanSin = p_mZdcQ1ShiftSinFull[mVzBin][iShift]->GetBinContent(binSin);

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
  p_mZdcSubEp1Res = new TProfile("p_mZdcSubEp1Res","p_mZdcSubEp1Res",mNumCentrality,-0.5,(double)mNumCentrality-0.5);
}

void StZdcEpManager::fillZdcResolution(TVector2 QEast, TVector2 QWest)
{
  double PsiEast = TMath::ATan2(-1.0*QEast.Y(),-1.0*QEast.X()); // flip the sign for Q1VecEast
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X());
  double res1Sub = TMath::Cos(PsiWest-PsiEast);
  p_mZdcSubEp1Res->Fill((double)mCent9,res1Sub);
}

void StZdcEpManager::writeZdcResolution()
{
  p_mZdcSubEp1Res->Write();
}

void StZdcEpManager::readZdcResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_ZdcEpResolution_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mZdcSubEp1Res = (TProfile*)file_mResolution->Get("p_mZdcSubEp1Res");

  for(int iCent = 0; iCent < 9; ++iCent) 
  {
    mZdcFullEp1ResVal[iCent] = 0.0;
    mZdcFullEp1ResErr[iCent] = 0.0;
  }

  TF1 *f_ZdcEpResFull = new TF1("f_ZdcEpResFull",funcZdcEpResFull,0,10,0);
  for(Int_t iCent = 0; iCent < mNumCentrality; iCent++)
  {
    double valResSub  = -999.9;
    double errResSub  = 1.0;
    double valResFull = -999.9;
    double errResFull = 1.0;
    double valResRaw = p_mZdcSubEp1Res->GetBinContent(iCent+1);
    double errResRaw = p_mZdcSubEp1Res->GetBinError(iCent+1);
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
    mZdcSubEp1ResVal[iCent]  = valResSub;
    mZdcSubEp1ResErr[iCent]  = errResSub;
    mZdcFullEp1ResVal[iCent] = valResFull;
    mZdcFullEp1ResErr[iCent] = errResFull;
  }

  file_mResolution->Close();
}

double StZdcEpManager::getZdcSubEpResVal(int cent9)
{
  return mZdcSubEp1ResVal[cent9];
}

double StZdcEpManager::getZdcSubEpResErr(int cent9)
{
  return mZdcSubEp1ResErr[cent9];
}

double StZdcEpManager::getZdcFullEpResVal(int cent9)
{
  return mZdcFullEp1ResVal[cent9];
}

double StZdcEpManager::getZdcFullEpResErr(int cent9)
{
  return mZdcFullEp1ResErr[cent9];
}
//---------------------------------------------------------------------------------
// Charged Hadron Directed Flow
void StZdcEpManager::initZdcFullEpFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    std::string proName = Form("p_mZdcFullEpDFlowCent%d",i_cent);
    p_mZdcFullEpDFlow[i_cent] = new TProfile(proName.c_str(),proName.c_str(),40,-6.0,6.0);
  }
}

void StZdcEpManager::fillZdcFullEpDFlow(double eta, double pt, double v1, double reweight)
{
  if(pt > 0.15 && pt < 2.0)
  { // pT cut for comparison
    p_mZdcFullEpDFlow[mCent9]->Fill(eta, v1, reweight);
  }
}

void StZdcEpManager::writeZdcFullEpFlow()
{
  for(int i_cent = 0; i_cent < 9; ++i_cent)
  {
    p_mZdcFullEpDFlow[i_cent]->Write();
  }
}
//---------------------------------------------------------------------------------
// raw EP
void StZdcEpManager::initZdcSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEp1RawEastCent%d",iCent);
    h_mZdcEp1RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEp1RawWestCent%d",iCent);
    h_mZdcEp1RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEp1RawFullCent%d",iCent);
    h_mZdcEp1RawFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEp1RawCorrCent%d",iCent);
    h_mZdcEp1RawCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),90,-1.0*TMath::Pi(),TMath::Pi(),90,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StZdcEpManager::fillZdcSubEpRaw(TVector2 QEast, TVector2 QWest, TVector2 QFull)
{
  double PsiEast = TMath::ATan2(-1.0*QEast.Y(),-1.0*QEast.X()); h_mZdcEp1RawEast[mCent9]->Fill(mRunIndex,PsiEast); // flip the sign for Q1VecEast
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcEp1RawWest[mCent9]->Fill(mRunIndex,PsiWest);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcEp1RawFull[mCent9]->Fill(mRunIndex,PsiFull);
  h_mZdcEp1RawCorr[mCent9]->Fill(PsiEast,PsiWest);
}

void StZdcEpManager::writeZdcSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEp1RawEast[iCent]->Write();
    h_mZdcEp1RawWest[iCent]->Write();
    h_mZdcEp1RawFull[iCent]->Write();
    h_mZdcEp1RawCorr[iCent]->Write();
  }
}

// recenter EP
void StZdcEpManager::initZdcSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEp1ReCtrEastCent%d",iCent);
    h_mZdcEp1ReCtrEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEp1ReCtrWestCent%d",iCent);
    h_mZdcEp1ReCtrWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEp1ReCtrFullCent%d",iCent);
    h_mZdcEp1ReCtrFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEp1ReCtrCorrCent%d",iCent);
    h_mZdcEp1ReCtrCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),90,-1.0*TMath::Pi(),TMath::Pi(),90,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StZdcEpManager::fillZdcSubEpReCtr(TVector2 QEast, TVector2 QWest, TVector2 QFull)
{
  double PsiEast = TMath::ATan2(-1.0*QEast.Y(),-1.0*QEast.X()); h_mZdcEp1ReCtrEast[mCent9]->Fill(mRunIndex,PsiEast); // flip the sign for Q1VecEast
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcEp1ReCtrWest[mCent9]->Fill(mRunIndex,PsiWest);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcEp1ReCtrFull[mCent9]->Fill(mRunIndex,PsiFull);
  h_mZdcEp1ReCtrCorr[mCent9]->Fill(PsiEast,PsiWest);
}

void StZdcEpManager::writeZdcSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEp1ReCtrEast[iCent]->Write();
    h_mZdcEp1ReCtrWest[iCent]->Write();
    h_mZdcEp1ReCtrFull[iCent]->Write();
    h_mZdcEp1ReCtrCorr[iCent]->Write();
  }
}

// shift EP
void StZdcEpManager::initZdcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEp1ShiftEastCent%d",iCent);
    h_mZdcEp1ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEp1ShiftWestCent%d",iCent);
    h_mZdcEp1ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEp1ShiftFullCent%d",iCent);
    h_mZdcEp1ShiftFull[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
    histName = Form("h_mZdcEp1ShiftCorrCent%d",iCent);
    h_mZdcEp1ShiftCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),90,-1.0*TMath::Pi(),TMath::Pi(),90,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StZdcEpManager::fillZdcSubEpShift(TVector2 QEast, TVector2 QWest, TVector2 QFull)
{
  double PsiEast = TMath::ATan2(-1.0*QEast.Y(),-1.0*QEast.X()); h_mZdcEp1ShiftEast[mCent9]->Fill(mRunIndex,PsiEast); // flip the sign for Q1VecEast
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcEp1ShiftWest[mCent9]->Fill(mRunIndex,PsiWest);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcEp1ShiftFull[mCent9]->Fill(mRunIndex,PsiFull);
  h_mZdcEp1ShiftCorr[mCent9]->Fill(PsiEast,PsiWest);
}

void StZdcEpManager::writeZdcSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEp1ShiftEast[iCent]->Write();
    h_mZdcEp1ShiftWest[iCent]->Write();
    h_mZdcEp1ShiftFull[iCent]->Write();
    h_mZdcEp1ShiftCorr[iCent]->Write();
  }
}

// shift Full EP
void StZdcEpManager::initZdcFullEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mZdcEp1ShiftFullCorrCent%d",iCent);
    h_mZdcEp1ShiftFullCorr[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,360,-1.0*TMath::Pi(),TMath::Pi());
  }
}

void StZdcEpManager::fillZdcFullEpShift(TVector2 QFull)
{
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); 
  h_mZdcEp1ShiftFullCorr[mCent9]->Fill(mRunIndex,PsiFull);
}

void StZdcEpManager::writeZdcFullEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mZdcEp1ShiftFullCorr[iCent]->Write();
  }
}
//---------------------------------------------------------------------------------
