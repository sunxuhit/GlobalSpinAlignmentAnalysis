#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"

#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "Utility/include/StSpinAlignmentFunctions.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"
// #include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

ClassImp(StTpcEpManager)

//---------------------------------------------------------------------------------
StTpcEpManager::StTpcEpManager(int beamType) : mType(beamType)
{
  // mEnergy = energy;
  clearTpcEp();
}

StTpcEpManager::~StTpcEpManager()
{
  /* */
}
//---------------------------------------------------------------------------------
void StTpcEpManager::clearTpcEp()
{
  mCent9 = -1;
  mRunIndex = -1;
  mVzBin = -1;

  mQCouEastRaw      = 0; // TPC EP East
  mQCouEastReCenter = 0;
  v_mQ2EastRaw.Set(0.0,0.0);
  v_mQ2EastReCenter.Set(0.0,0.0);

  mQCouWestRaw      = 0; // TPC EP West
  mQCouWestReCenter = 0;
  v_mQ2WestRaw.Set(0.0,0.0);
  v_mQ2WestReCenter.Set(0.0,0.0);
}

void StTpcEpManager::initTpcEp(int cent9, int runIndex, int vzBin)
{
  mCent9    = cent9;
  mRunIndex = runIndex;
  mVzBin    = vzBin;
}
//---------------------------------------------------------------------------------


void StTpcEpManager::readReCenterCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ReCenterPar/file_%s_TpcReCenterPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());

  file_mReCenterPar = TFile::Open(inputFile.c_str());
}
//---------------------------------------------------------------------------------

void StTpcEpManager::readShiftCorr()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ShiftParameter/file_%s_TpcShiftPar.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mShiftPar = TFile::Open(InPutFile_Shift.Data());
}

//---------------------------------------------------------------------------------

TVector2 StTpcEpManager::calq2Vector(StPicoTrack *picoTrack)
{
  const double phi = picoTrack->pMom().Phi(); // -pi to pi
  TVector2 q2Vector(0.0,0.0);

  const double q2x = TMath::Cos(2.0*phi);
  const double q2y = TMath::Sin(2.0*phi);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

double StTpcEpManager::getWeight(StPicoTrack *picoTrack)
{
  const double pt = picoTrack->pMom().Perp();
  double wgt;
  if(pt <= anaUtils::mPrimPtEpWeight[mType])
  {
    wgt = pt;
  }
  if(pt > anaUtils::mPrimPtEpWeight[mType])
  {
    wgt = anaUtils::mPrimPtEpWeight[mType];
  }

  return wgt;
}
//---------------------------------------------------------------------------------

void StTpcEpManager::addTrackEastRaw(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2EastRaw += wgt*calq2Vector(picoTrack);
  mQCouEastRaw++;
}

void StTpcEpManager::addTrackEast(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2EastReCenter += wgt*(calq2Vector(picoTrack) - getReCenterParEast());
  mQCouEastReCenter++;
}

TVector2 StTpcEpManager::getReCenterParEast()
{
  double mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_East",mOrder.Data(),mVStr[mVzBin].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_East",mOrder.Data(),mVStr[mVzBin].Data());

  TProfile2D *p_x = (TProfile2D*)file_mReCenterPar->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)file_mReCenterPar->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)mRunIndex,(Double_t)mCent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)mRunIndex,(Double_t)mCent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}


//---------------------------------------------------------------------------------

void StTpcEpManager::addTrackWestRaw(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2WestRaw += wgt*calq2Vector(picoTrack);
  mQCouWestRaw++;
}

void StTpcEpManager::addTrackWest(StPicoTrack *picoTrack)
{
  const double wgt = getWeight(picoTrack);
  v_mQ2WestReCenter += wgt*(calq2Vector(picoTrack) - getReCenterParWest());
  mQCouWestReCenter++;
}

TVector2 StTpcEpManager::getReCenterParWest()
{
  double mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_West",mOrder.Data(),mVStr[mVzBin].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_West",mOrder.Data(),mVStr[mVzBin].Data());

  TProfile2D *p_x = (TProfile2D*)file_mReCenterPar->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)file_mReCenterPar->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)mRunIndex,(Double_t)mCent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)mRunIndex,(Double_t)mCent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}
//---------------------------------------------------------------------------------

void StTpcEpManager::print(TVector2 vector)
{
  cout << "qx = " << vector.X() << endl;
  cout << "qy = " << vector.Y() << endl;
  cout << endl;
}

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// 2nd
#if 0
TVector2 StTpcEpManager::calPsi2_East_EP(int k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  double Qx = v_mQ2EastReCenter.X();
  double Qy = v_mQ2EastReCenter.Y();
  double Psi = TMath::ATan2(Qy,Qx)/2.0;
  double Psi_Sin = TMath::Sin(globCons::mShiftOrder[k]*Psi);
  double Psi_Cos = TMath::Cos(globCons::mShiftOrder[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTpcEpManager::calPsi2_West_EP(int k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  double Qx = v_mQ2WestReCenter.X();
  double Qy = v_mQ2WestReCenter.Y();
  double Psi = TMath::ATan2(Qy,Qx)/2.0;
  double Psi_Sin = TMath::Sin(globCons::mShiftOrder[k]*Psi);
  double Psi_Cos = TMath::Cos(globCons::mShiftOrder[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}
#endif

//---------------------------------------------------------------------------------

double StTpcEpManager::angleShift(double PsiRaw)
{
  double PsiCorr = PsiRaw;
  if(PsiRaw > 0.5*TMath::Pi())
  {
    PsiCorr = PsiRaw - TMath::Pi();
  }
  if(PsiRaw < -0.5*TMath::Pi())
  {
    PsiCorr = PsiRaw + TMath::Pi();
  }

  return PsiCorr;
}

//---------------------------------------------------------------------------------
double StTpcEpManager::calShiftAngle2East()
{
  double PsiReCenter = TMath::ATan2(v_mQ2EastReCenter.Y(),v_mQ2EastReCenter.X())/2.0; // -pi/2 to pi/2
  double deltaPsi = 0.0;
  double PsiShift;

  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    int binCos     = p_mTpcQShiftCosEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanCos = p_mTpcQShiftCosEast[mVzBin][iShift]->GetBinContent(binCos);

    int binSin     = p_mTpcQShiftSinEast[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanSin = p_mTpcQShiftSinEast[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(((double)iShift+1.0)*PsiReCenter)); // TODO: update for 2nd EP
  }

  double PsiShiftRaw = PsiReCenter + deltaPsi;
  PsiShift = angleShift(PsiShiftRaw);

  return PsiShift;
}

double StTpcEpManager::calShiftAngle2West()
{
  double PsiReCenter = TMath::ATan2(v_mQ2WestReCenter.Y(),v_mQ2WestReCenter.X())/2.0;
  double deltaPsi = 0.0;
  double PsiShift;

  for(Int_t iShift = 0; iShift < mNumShiftCorr; ++iShift) // Shift Order loop
  {
    int binCos     = p_mZdcQShiftCosWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanCos = p_mZdcQShiftCosWest[mVzBin][iShift]->GetBinContent(binCos);

    int binSin     = p_mZdcQShiftSinWest[mVzBin][iShift]->FindBin((double)mRunIndex,(double)mCent9);
    double meanSin = p_mZdcQShiftSinWest[mVzBin][iShift]->GetBinContent(binSin);

    deltaPsi += (2.0/((double)iShift+1.0))*(-1.0*meanSin*TMath::Cos(((double)iShift+1.0)*PsiReCenter)+meanCos*TMath::Sin(((double)iShift+1.0)*PsiReCenter)); // TODO: update for 2nd EP
  }

  double PsiShiftRaw = PsiReCenter + deltaPsi;
  PsiShift = angleShift(PsiShiftRaw);

  return PsiShift;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::readResolution()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_%s_TpcEpResolution.root",globCons::str_mBeamType[mType].c_str(),g  lobCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.Data());
  p_mTpcEpResolution = (TProfile*)file_mResolution->Get("p_mRes2_Sub");

  for(int iCent = 0; iCent < 9; ++iCent)
  {
    double valResSub, errResSub;
    double valResRaw = p_mTpcEpResolution->GetBinContent(iCent+1);
    double errResRaw = p_mTpcEpResolution->GetBinError(iCent+1);
    if(valResRaw <= 0)
    {
      valResSub = -999.9;
      errResSub = 1.0;
    }
    else
    {
      valResSub = TMath::Sqrt(valResRaw);
      errResSub = errResRaw/(2*valResSub);
    }

    mResolutionVal[iCent] = valResSub;
    mResolutionErr[iCent] = errResSub;
  }
}

double StTpcEpManager::getResolutionVal(int cent9)
{
  return mResolutionVal[cent9];
}

double StTpcEpManager::getResolutionErr(int cent9)
{
  return mResolutionErr[cent9];
}

//---------------------------------------------------------------------------------
TVector2 StTpcEpManager::getQVector(int epMode) // east/west
{
  if(epMode == 0) return v_mQ2EastReCenter;
  if(epMode == 1) return v_mQ2WestReCenter;
}

TVector2 StTpcEpManager::getQVectorRaw(int epMode) // east/west
{
  if(epMode == 0) return v_mQ2EastRaw;
  if(epMode == 1) return v_mQ2WestRaw;
}

int StTpcEpManager::getNumTrack(int epMode) // east/west
{
  if(epMode == 0) return mQCouEastReCenter;
  if(epMode == 1) return mQCouWestReCenter;
}
//---------------------------------------------------------------------------------
