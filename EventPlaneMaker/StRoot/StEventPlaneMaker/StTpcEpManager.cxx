#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"

#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StEventPlaneMaker/StTpcEpManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

Double_t Resolution_Full(Double_t *x_val, Double_t *par)
{
  Double_t y;
  Double_t chi = x_val[0];
  Double_t arg = chi*chi/4.0;
  Double_t norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

ClassImp(StTpcEpManager)

TString StTpcEpManager::mVStr[2] = {"pos","neg"};
TString StTpcEpManager::mOrder = "2nd";
//---------------------------------------------------------------------------------

StTpcEpManager::StTpcEpManager(Int_t energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StTpcEpManager::~StTpcEpManager()
{
  /* */
}

//---------------------------------------------------------------------------------

void StTpcEpManager::InitReCenterCorrection()
{
  TString InPutFile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ReCenterParameter/file_%s_ReCenterPar.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamEnergy[mEnergy].c_str());

  mInPutFile = TFile::Open(InPutFile.Data());

  mQ2Vector_East_EP.Set(0.0,0.0);
  mQCounter_East = 0;

  mQ2Vector_West_EP.Set(0.0,0.0);
  mQCounter_West = 0;

  mQ2Vector_Full_EP.Set(0.0,0.0);
  mQCounter_Full = 0;
  mQCounter_Full_East = 0;
  mQCounter_Full_West = 0;

  mQ2Vector_A_EP.Set(0.0,0.0);
  mQCounter_A = 0;

  mQ2Vector_B_EP.Set(0.0,0.0);
  mQCounter_B = 0;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::InitShiftCorrection()
{
  TString InPutFile_Shift = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/ShiftParameter/file_%s_ShiftPar.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamEnergy[mEnergy].c_str());
  mInPutFile_Shift = TFile::Open(InPutFile_Shift.Data());
}

//---------------------------------------------------------------------------------

bool StTpcEpManager::passTrackEtaEast(StPicoTrack *track) // neg
{
  Float_t eta = track->pMom().pseudoRapidity();

  // eta cut
  // eta_gap between two sub event plane is 2*mEta_Gap
  if(!(eta > -1.0*vmsa::mEtaMax && eta < -1.0*vmsa::mEta_Gap))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StTpcEpManager::passTrackEtaWest(StPicoTrack *track) // pos 
{
  Float_t eta = track->pMom().pseudoRapidity();

  // eta cut
  // eta_gap between two sub event plane is 2*mEta_Gap
  if(!(eta > vmsa::mEta_Gap && eta < vmsa::mEtaMax))
  {
    return kFALSE;
  }

  return kTRUE;
}
//---------------------------------------------------------------------------------

bool StTpcEpManager::passTrackFull(StPicoTrack *track) // Full Event Plane pt Cut
{
  // pt cut 0.2 - 2.0 GeV/c
  Float_t pt = track->pMom().perp();
  if(!(pt > vmsa::mPrimPtMin[mEnergy] && pt < vmsa::mPrimPtMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------
TVector2 StTpcEpManager::calq2Vector(StPicoTrack *track)
{
  const Float_t phi = track->pMom().phi();
  TVector2 q2Vector(0.0,0.0);

  const Float_t q2x = TMath::Cos(2.0*phi);
  const Float_t q2y = TMath::Sin(2.0*phi);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

Float_t StTpcEpManager::getWeight(StPicoTrack *track)
{
  Float_t pt = track->pMom().perp();
  Float_t w;
  if(pt <= vmsa::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > vmsa::mPrimPtWeight)
  {
    w = vmsa::mPrimPtWeight;
  }

  return w;
}
//---------------------------------------------------------------------------------

TVector2 StTpcEpManager::getReCenterPar_East( Int_t Cent9, Int_t RunIndex, Int_t vz_sign)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_East",mOrder.Data(),mVStr[vz_sign].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_East",mOrder.Data(),mVStr[vz_sign].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

//---------------------------------------------------------------------------------

TVector2 StTpcEpManager::getReCenterPar_West(Int_t Cent9, Int_t RunIndex, Int_t vz_sign)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_West",mOrder.Data(),mVStr[vz_sign].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_West",mOrder.Data(),mVStr[vz_sign].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

//---------------------------------------------------------------------------------

TVector2 StTpcEpManager::getReCenterPar_Full(Int_t Cent9, Int_t RunIndex, Int_t vz_sign)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_Full",mOrder.Data(),mVStr[vz_sign].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_Full",mOrder.Data(),mVStr[vz_sign].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::addTrack_EastRaw(StPicoTrack *track, Int_t Cent9, Int_t RunIndex)
{
  Float_t w = getWeight(track);
  mQ2Vector_EastRaw_EP += w*calq2Vector(track);
  mQCounter_RawEast++;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::addTrack_WestRaw(StPicoTrack *track, Int_t Cent9, Int_t RunIndex)
{
  Float_t w = getWeight(track);
  mQ2Vector_WestRaw_EP += w*calq2Vector(track);
  mQCounter_RawWest++;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::addTrack_FullRaw(StPicoTrack *track, Int_t Cent9, Int_t RunIndex)
{
  Float_t w = getWeight(track);
  mQ2Vector_FullRaw_EP += w*calq2Vector(track);
  mQCounter_RawFull++;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::addTrack_East(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i)
{
  Float_t w = getWeight(track);
  mQ2Vector_East_EP += w*(calq2Vector(track) - getReCenterPar_East(Cent9,RunIndex,i));

  mQCounter_East++;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::addTrack_West(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i)
{
  Float_t w = getWeight(track);
  mQ2Vector_West_EP += w*(calq2Vector(track) - getReCenterPar_West(Cent9,RunIndex,i));

  mQCounter_West++;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::addTrack_Full(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i)
{
  Float_t w = getWeight(track);
  mQ2Vector_Full_EP += w*(calq2Vector(track) - getReCenterPar_Full(Cent9,RunIndex,i));

  mQCounter_Full++;

  Float_t eta = track->pMom().pseudoRapidity();
  if(eta >= 0.0)
  {
    mQCounter_Full_West++;
  }
  if(eta < 0.0)
  {
    mQCounter_Full_East++;
  }
}

void StTpcEpManager::addTrack_A(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i)
{
  Float_t w = getWeight(track);
  mQ2Vector_A_EP += w*(calq2Vector(track) - getReCenterPar_Full(Cent9,RunIndex,i));

  mQCounter_A++;
}

void StTpcEpManager::addTrack_B(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i)
{
  Float_t w = getWeight(track);
  mQ2Vector_B_EP += w*(calq2Vector(track) - getReCenterPar_Full(Cent9,RunIndex,i));

  mQCounter_B++;
}

void StTpcEpManager::Randomization()
{
  TRandom3 Ran;
  TVector2 Q2Switch_EP;
  Int_t CSwitch;
  Ran.SetSeed();
  Float_t ran = Ran.Rndm(); // random number between [0,1]
  if(ran < 0.5)
  {
    // switch Event Plane Q Vector
    Q2Switch_EP = mQ2Vector_A_EP;
    mQ2Vector_A_EP = mQ2Vector_B_EP;
    mQ2Vector_B_EP = Q2Switch_EP;

    // switch Counter
    CSwitch = mQCounter_A;
    mQCounter_A = mQCounter_B;
    mQCounter_B = CSwitch;
  }
}

//---------------------------------------------------------------------------------

void StTpcEpManager::print(TVector2 vector)
{
  cout << "qx = " << vector.X() << endl;
  cout << "qy = " << vector.Y() << endl;
  cout << endl;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::clear()
{
  mQ2Vector_EastRaw_EP.Set(0.0,0.0);
  mQCounter_RawEast = 0;
  mQ2Vector_East_EP.Set(0.0,0.0);
  mQCounter_East = 0;

  mQ2Vector_WestRaw_EP.Set(0.0,0.0);
  mQCounter_RawWest = 0;
  mQ2Vector_West_EP.Set(0.0,0.0);
  mQCounter_West = 0;
  
  mQ2Vector_FullRaw_EP.Set(0.0,0.0);
  mQCounter_RawFull = 0;
  mQ2Vector_Full_EP.Set(0.0,0.0);
  mQCounter_Full = 0;
  mQCounter_Full_East = 0;
  mQCounter_Full_West = 0;

  mQ2Vector_A_EP.Set(0.0,0.0);
  mQCounter_A = 0;

  mQ2Vector_B_EP.Set(0.0,0.0);
  mQCounter_B = 0;
}

//---------------------------------------------------------------------------------

bool StTpcEpManager::passTrackEtaNumRawCut()
{
  if(!(mQCounter_RawEast > vmsa::mTrackMin && mQCounter_RawWest > vmsa::mTrackMin))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StTpcEpManager::passTrackFullNumRawCut()
{
  if(!(mQCounter_RawFull > vmsa::mTrackMin_Full))
  {
    return kFALSE;
  }
  
  return kTRUE;
}

bool StTpcEpManager::passTrackEtaNumCut()
{
  if(!(mQCounter_East > vmsa::mTrackMin && mQCounter_West > vmsa::mTrackMin))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StTpcEpManager::passTrackFullNumCut()
{
  if(!(mQCounter_Full > vmsa::mTrackMin_Full && mQCounter_Full_East > 0 && mQCounter_Full_West > 0))
  {
    return kFALSE;
  }
  
  return kTRUE;
}

//---------------------------------------------------------------------------------
// 2nd
TVector2 StTpcEpManager::calPsi2_East_EP(Int_t k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_East_EP.X();
  Float_t Qy = mQ2Vector_East_EP.Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(vmsa::mShiftOrder[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(vmsa::mShiftOrder[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTpcEpManager::calPsi2_West_EP(Int_t k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_West_EP.X();
  Float_t Qy = mQ2Vector_West_EP.Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(vmsa::mShiftOrder[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(vmsa::mShiftOrder[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

TVector2 StTpcEpManager::calPsi2_Full_EP(Int_t k) // k = ShiftOrder
{
  TVector2 PsiVector(0.0,0.0); 

  Float_t Qx = mQ2Vector_Full_EP.X();
  Float_t Qy = mQ2Vector_Full_EP.Y();
  Float_t Psi = TMath::ATan2(Qy,Qx)/2.0;
  Float_t Psi_Sin = TMath::Sin(vmsa::mShiftOrder[k]*Psi);
  Float_t Psi_Cos = TMath::Cos(vmsa::mShiftOrder[k]*Psi);

  PsiVector.Set(Psi_Cos,Psi_Sin);

  return PsiVector;
}

//---------------------------------------------------------------------------------

Float_t StTpcEpManager::AngleShift(Float_t Psi_raw)
{
  Float_t Psi_Corr = Psi_raw;
  if(Psi_raw > 0.5*TMath::Pi())
  {
    Psi_Corr = Psi_raw - TMath::Pi();
  }
  if(Psi_raw < -0.5*TMath::Pi())
  {
    Psi_Corr = Psi_raw + TMath::Pi();
  }

  return Psi_Corr;
}

//---------------------------------------------------------------------------------
Float_t StTpcEpManager::calShiftAngle2East_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_East_EP.Y(),mQ2Vector_East_EP.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_East",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_East",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

Float_t StTpcEpManager::calShiftAngle2West_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_West_EP.Y(),mQ2Vector_West_EP.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_West",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_West",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

Float_t StTpcEpManager::calShiftAngle2A_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_A_EP.Y(),mQ2Vector_A_EP.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

Float_t StTpcEpManager::calShiftAngle2B_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_B_EP.Y(),mQ2Vector_B_EP.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

Float_t StTpcEpManager::calShiftAngle2Full_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_Full_EP.Y(),mQ2Vector_Full_EP.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

Float_t StTpcEpManager::calShiftAngle2Full_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track)
{
  TVector2 QVector_sub = mQ2Vector_Full_EP;
  if(passTrackFull(track))
  {
    Float_t w = getWeight(track);
    QVector_sub = mQ2Vector_Full_EP - w*(calq2Vector(track) - getReCenterPar_Full(Cent9,runIndex,vz_sign));
  }
  Float_t Psi_ReCenter = TMath::ATan2(QVector_sub.Y(),QVector_sub.X())/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  Float_t Psi_Shift;

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(vmsa::mShiftOrder[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(vmsa::mShiftOrder[k]*Psi_ReCenter));
  }

  Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
  Psi_Shift = AngleShift(Psi_Shift_raw);

  return Psi_Shift;
}

//---------------------------------------------------------------------------------

void StTpcEpManager::InitResolutionCorr()
{
  TString InPutFile_Res = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Resolution/file_%s_Resolution.root",vmsa::mBeamEnergy[mEnergy].c_str(),vmsa::mBeamEnergy[mEnergy].c_str());
  mInPutFile_Res = TFile::Open(InPutFile_Res.Data());
}

Float_t StTpcEpManager::getResolution2_EP(Int_t Cent9)
{
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get("p_mRes2_Sub");
  Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(Cent9));
  if(Res_raw <= 0)
  {
    return -999.9;
  }
  else
  {
    Float_t Res = TMath::Sqrt(Res_raw);
    return Res;
  }
}

Float_t StTpcEpManager::getResolution2_Full_EP(Int_t Cent9)
{
  TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get("p_mRes2_Ran");
  Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(Cent9));
  if(Res_raw <= 0)
  {
    return -999.9;
  }
  else
  {
    Float_t Res_sub = TMath::Sqrt(Res_raw);
    TF1 *f_res = new TF1("f_res",Resolution_Full,0,10,0);
    Float_t chi_sub = f_res->GetX(Res_sub);
    Float_t chi_full = chi_sub*TMath::Sqrt(2.0);
    Float_t Res_full = f_res->Eval(chi_full);
    return Res_full;
  }
}
//---------------------------------------------------------------------------------
TVector2 StTpcEpManager::getQVector(Int_t l) // east/west
{
  if(l == 0) return mQ2Vector_East_EP;
  if(l == 1) return mQ2Vector_West_EP;
  if(l == 2) return mQ2Vector_Full_EP;
  if(l == 3) return mQ2Vector_A_EP;
  if(l == 4) return mQ2Vector_B_EP;
}

TVector2 StTpcEpManager::getQVectorRaw(Int_t l) // east/west
{
  if(l == 0) return mQ2Vector_EastRaw_EP;
  if(l == 1) return mQ2Vector_WestRaw_EP;
  if(l == 2) return mQ2Vector_FullRaw_EP;
}

Int_t StTpcEpManager::getNumTrack(Int_t l) // east/west
{
  if(l == 0) return mQCounter_East;
  if(l == 1) return mQCounter_West;
  if(l == 2) return mQCounter_Full;
  if(l == 3) return mQCounter_Full_East;
  if(l == 4) return mQCounter_Full_West;
}
//---------------------------------------------------------------------------------
