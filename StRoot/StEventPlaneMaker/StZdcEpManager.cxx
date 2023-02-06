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
	mGainFactor[iEastWest][iVertHori][iSlat] = -999.9;
      }
    }
  }
  mCenterEastVertical   = -999.9;
  mCenterEastHorizontal = -999.9;
  mCenterWestVertical   = -999.9;
  mCenterWestHorizontal = -999.9;
  for(int iCent = 0; iCent < 9; ++iCent) mResolution[iCent] = 0.0;
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

void StZdcEpManager::fillZdcGainCorr(int iEastWest, int iVertHori, int iSlat, int runIndex, double zdcsmd)
{
  h_mZdcGainCorr[iEastWest][iVertHori][iSlat]->Fill((double)runIndex,zdcsmd);
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
  if(mode > 2)  qVector = applyZdcSmdShiftCorrEast(qVector);

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
  if(mode > 2) qVector = applyZdcSmdShiftCorrWest(qVector);

  return qVector;
}

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

void StZdcEpManager::fillZdcRawEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex)
{
  double PsiEast = TMath::ATan2(QEast.Y(),QEast.X()); h_mZdcRawEpEast[Cent9]->Fill(PsiEast,runIndex);
  double PsiWest = TMath::ATan2(QWest.Y(),QWest.X()); h_mZdcRawEpWest[Cent9]->Fill(PsiWest,runIndex);
  double PsiFull = TMath::ATan2(QFull.Y(),QFull.X()); h_mZdcRawEpFull[Cent9]->Fill(PsiFull,runIndex);
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
// ReCenter Correction
void StZdcEpManager::initZdcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string ProName = Form("p_mZdcQEastVertVz%d",iVz);
    p_mZdcQEastVertical[iVz] = new TProfile2D(ProName.c_str(),ProName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    ProName = Form("p_mZdcQEastHoriVz%d",iVz);
    p_mZdcQEastHorizontal[iVz] = new TProfile2D(ProName.c_str(),ProName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);

    ProName = Form("p_mZdcQWestVertVz%d",iVz);
    p_mZdcQWestVertical[iVz] = new TProfile2D(ProName.c_str(),ProName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    ProName = Form("p_mZdcQWestHoriVz%d",iVz);
    p_mZdcQWestHorizontal[iVz] = new TProfile2D(ProName.c_str(),ProName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
  }
}

void StZdcEpManager::fillZdcReCenterEast(TVector2 qVector, int Cent9, int RunIndex, int vzBin)
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  p_mZdcQEastVertical[vzBin]->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mZdcQEastHorizontal[vzBin]->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StZdcEpManager::fillZdcReCenterWest(TVector2 qVector, int Cent9, int RunIndex, int vzBin)
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  p_mZdcQWestVertical[vzBin]->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mZdcQWestHorizontal[vzBin]->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StZdcEpManager::writeZdcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mZdcQEastVertical[iVz]->Write();
    p_mZdcQEastHorizontal[iVz]->Write();
    p_mZdcQWestVertical[iVz]->Write();
    p_mZdcQWestHorizontal[iVz]->Write();
  }
}
void StZdcEpManager::readReCenterCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_ZdcReCenterPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mReCenterPar = TFile::Open(inputFile.c_str());
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string proName = Form("p_mZdcQEastVertVz%d",iVz);
    p_mZdcQEastVertical[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mZdcQEastHoriVz%d",iVz);
    p_mZdcQEastHorizontal[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());

    proName = Form("p_mZdcQWestVertVz%d",iVz);
    p_mZdcQWestVertical[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
    proName = Form("p_mZdcQWestHoriVz%d",iVz);
    p_mZdcQWestHorizontal[iVz] = (TProfile2D*)file_mReCenterPar->Get(proName.c_str());
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

//---------------------------------------------------------------------------------
void StZdcEpManager::readShiftCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftParameter/file_%s_ZdcShiftPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(inputFile.c_str());

  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQEastVz%dCos%d",iVz,iShift);
      p_mZdcQEastCos[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mZdcQEastVz%dSin%d",iVz,iShift);
      p_mZdcQEastSin[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());

      proName = Form("p_mZdcQWestVz%dCos%d",iVz,iShift);
      p_mZdcQWestCos[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mZdcQWestVz%dSin%d",iVz,iShift);
      p_mZdcQWestSin[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

TVector2 StZdcEpManager::applyZdcSmdShiftCorrEast(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  double PsiReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  double deltaPsi = 0.0;
  double PsiShift;

  for(Int_t iShift = 0; iShift < 20; ++iShift) // Shift Order loop
  {
    int binCos     = p_mZdcQEastCos[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanCos = p_mZdcQEastCos[mVzBin][iShift]->GetBinContent(binCos);

    int binSin     = p_mZdcQEastSin[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanSin = p_mZdcQEastSin[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(((double)iShift+1.0)*PsiReCenter));
  }

  double PsiShiftRaw = PsiReCenter + deltaPsi;
  PsiShift = angleShift(PsiShiftRaw);

  qVecShift.Set(TMath::Cos(PsiShift),TMath::Sin(PsiShift));

  return qVecShift;
}

TVector2 StZdcEpManager::applyZdcSmdShiftCorrWest(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  double PsiReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  double deltaPsi = 0.0;
  double PsiShift;

  for(Int_t iShift = 0; iShift < 20; ++iShift) // Shift Order loop
  {
    int binCos     = p_mZdcQWestCos[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanCos = p_mZdcQWestCos[mVzBin][iShift]->GetBinContent(binCos);

    int binSin     = p_mZdcQWestSin[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanSin = p_mZdcQWestSin[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(((double)iShift+1.0)*PsiReCenter));
  }

  double PsiShiftRaw = PsiReCenter + deltaPsi;
  PsiShift = angleShift(PsiShiftRaw);

  qVecShift.Set(TMath::Cos(PsiShift),TMath::Sin(PsiShift));

  return qVecShift;
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


void StZdcEpManager::readShiftCorrFull()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftParameter/file_%s_ZdcShiftParFull.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(inputFile.c_str());

  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    for(int iShift = 0; iShift < 20; ++iShift) // Shift Order
    {
      std::string proName = Form("p_mZdcQFullVz%dCos%d",iVz,iShift);
      p_mZdcQFullCos[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
      proName = Form("p_mZdcQFullVz%dSin%d",iVz,iShift);
      p_mZdcQFullSin[iVz][iShift] = (TProfile2D*)file_mShiftPar->Get(proName.c_str());
    }
  }
}

TVector2 StZdcEpManager::applyZdcSmdShiftCorrFull(TVector2 qVector)
{
  TVector2 qVecShift(0.0,0.0);
  double PsiReCenter = TMath::ATan2(qVector.Y(),qVector.X());
  double deltaPsi = 0.0;
  double PsiShift;

  for(Int_t iShift = 0; iShift < 20; ++iShift) // Shift Order loop
  {
    int binCos     = p_mZdcQFullCos[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanCos = p_mZdcQFullCos[mVzBin][iShift]->GetBinContent(binCos);

    int binSin     = p_mZdcQFullSin[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanSin = p_mZdcQFullSin[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(((double)iShift+1.0)*PsiReCenter));
  }

  double PsiShiftRaw = PsiReCenter + deltaPsi;
  PsiShift = angleShift(PsiShiftRaw);

  qVecShift.Set(TMath::Cos(PsiShift),TMath::Sin(PsiShift));

  return qVecShift;
}

TVector2 StZdcEpManager::getQFull(TVector2 QEast, TVector2 QWest)
{
  TVector2 qVector = QWest-QEast;
  TVector2 qVecShift = applyZdcSmdShiftCorrFull(qVector);

  return qVecShift;
}

//---------------------------------------------------------------------------------

void StZdcEpManager::readResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_%s_ZdcEpResolution.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  p_mZdcEpResolution = (TProfile*)file_mResolution->Get("p_mZdcEpResolution");
}

void StZdcEpManager::calResolution()
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
    mResolution[iCent] = valResFull;
  }
}

double StZdcEpManager::getResolution(int Cent9)
{
  return mResolution[Cent9];
}
//---------------------------------------------------------------------------------
