#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TF1.h"

#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StEventPlaneMaker/StMixEpManager.h"

double funcMixEp1Res1(double *x_val, double *par)
{
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  double res1Full = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return res1Full;
}

double funcMixEp1Res2(double *x_val, double *par)
{
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi())/(2.0*TMath::Sqrt(2.0));
  double besselOneHalf = TMath::Sqrt(2.0*arg/TMath::Pi())*TMath::SinH(arg)/arg;
  double besselThreeHalf = TMath::Sqrt(2.0*arg/TMath::Pi())*(TMath::CosH(arg)/arg-TMath::SinH(arg)/(arg*arg));

  double res12Sub = norm*chi*TMath::Exp(-1.0*arg)*(besselOneHalf+besselThreeHalf);

  return res12Sub;
}

ClassImp(StMixEpManager)

//---------------------------------------------------------------------------------
StMixEpManager::StMixEpManager(int beamType) : mType(beamType)
{
  clearMixEpManager();
}

StMixEpManager::~StMixEpManager()
{
  /* */
}
//---------------------------------------------------------------------------------
void StMixEpManager::clearMixEpManager()
{
  mCent9    = -1;
  mRunIndex = -1;
  mVzBin    = -1;
}

void StMixEpManager::initMixEpManager(int cent9, int runIndex, int vzBin)
{
  mCent9    = cent9;
  mRunIndex = runIndex;
  mVzBin    = vzBin;
}
//---------------------------------------------------------------------------------
// 3 Sub Event Plane Resolution
void StMixEpManager::initMixEpRes()
{
  for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
  {
    std::string proName = Form("p_mMixSubEp1ResGrp%d",iGrp); // 1st EP
    p_mMixSubEp1Res[iGrp] = new TProfile(proName.c_str(),proName.c_str(),mNumCentrality,-0.5,(double)mNumCentrality-0.5);
  }
}

void StMixEpManager::fillMixEpRes(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest)
{
  p_mMixSubEp1Res[0]->Fill((double)mCent9,TMath::Cos(Psi1EpdGrp0-Psi1TpcEast)); // 0: EpdEpGrp0 vs. TpcEpEast
  p_mMixSubEp1Res[1]->Fill((double)mCent9,TMath::Cos(Psi1EpdGrp0-Psi1TpcWest)); // 1: EpdEpGrp0 vs. TpcEpWest
  p_mMixSubEp1Res[2]->Fill((double)mCent9,TMath::Cos(Psi1EpdGrp1-Psi1TpcEast)); // 2: EpdEpGrp1 vs. TpcEpEast
  p_mMixSubEp1Res[3]->Fill((double)mCent9,TMath::Cos(Psi1EpdGrp1-Psi1TpcWest)); // 3: EpdEpGrp1 vs. TpcEpWest
  p_mMixSubEp1Res[4]->Fill((double)mCent9,TMath::Cos(Psi1EpdGrp0-Psi1EpdGrp1)); // 4: EpdEpGrp0 vs. EpdEpGrp1
  p_mMixSubEp1Res[5]->Fill((double)mCent9,TMath::Cos(Psi1TpcEast-Psi1TpcWest)); // 5: TpcEpEast vs. TpcEpWest
}

void StMixEpManager::writeMixEpRes()
{
  for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
  {
    p_mMixSubEp1Res[iGrp]->Write();
  }
}

void StMixEpManager::readMixEpRes()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/Resolution/file_MixEpResolution_%s.root",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  file_mResolution = TFile::Open(inputFile.c_str());
  for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
  {
    std::string proName = Form("p_mMixSubEp1ResGrp%d",iGrp); // 1st EP
    p_mMixSubEp1Res[iGrp] = (TProfile*)file_mResolution->Get(proName.c_str())->Clone();
  }

  double valMixEp[mNumCentrality][mNumEpGroup];
  double errMixEp[mNumCentrality][mNumEpGroup];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      valMixEp[iCent][iGrp] = p_mMixSubEp1Res[iGrp]->GetBinContent(iCent+1);
      errMixEp[iCent][iGrp] = p_mMixSubEp1Res[iGrp]->GetBinError(iCent+1);
    }
  }

  // double subRes1Grp0 = TMath::Cos(Psi1EpdGrp0-Psi1EpdGrp1)*TMath::Cos(Psi1EpdGrp0-Psi1TpcWest)/TMath::Cos(Psi1EpdGrp1-Psi1TpcWest);
  // double subRes1Grp1 = TMath::Cos(Psi1EpdGrp0-Psi1EpdGrp1)*TMath::Cos(Psi1EpdGrp0-Psi1TpcEast)/TMath::Cos(Psi1EpdGrp1-Psi1TpcEast);
  // double subRes1Grp2 = TMath::Cos(Psi1EpdGrp0-Psi1TpcEast)*TMath::Cos(Psi1EpdGrp0-Psi1TpcWest)/TMath::Cos(Psi1TpcEast-Psi1TpcWest);
  // double subRes1Grp3 = TMath::Cos(Psi1EpdGrp1-Psi1EpdGrp0)*TMath::Cos(Psi1EpdGrp1-Psi1TpcWest)/TMath::Cos(Psi1EpdGrp0-Psi1TpcWest);
  // double subRes1Grp4 = TMath::Cos(Psi1EpdGrp1-Psi1EpdGrp0)*TMath::Cos(Psi1EpdGrp1-Psi1TpcEast)/TMath::Cos(Psi1EpdGrp0-Psi1TpcEast);
  // double subRes1Grp5 = TMath::Cos(Psi1EpdGrp1-Psi1TpcEast)*TMath::Cos(Psi1EpdGrp1-Psi1TpcWest)/TMath::Cos(Psi1TpcEast-Psi1TpcWest);

  // 0: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpWest (default) | 1: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpEast | 2: EpdEpGrp0 vs. TpcEpEast && TpcEpWest
  // 3: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpWest (mainSys) | 4: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpEast | 5: EpdEpGrp1 vs. TpcEpEast && TpcEpWest
  double valRes1Temp[mNumCentrality][mNumEpGroup];
  double errRes1Temp[mNumCentrality][mNumEpGroup];
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      valRes1Temp[iCent][0] = valMixEp[iCent][4]*valMixEp[iCent][1]/valMixEp[iCent][3];
      valRes1Temp[iCent][1] = valMixEp[iCent][4]*valMixEp[iCent][0]/valMixEp[iCent][2];
      valRes1Temp[iCent][2] = valMixEp[iCent][0]*valMixEp[iCent][1]/valMixEp[iCent][5];
      valRes1Temp[iCent][3] = valMixEp[iCent][4]*valMixEp[iCent][3]/valMixEp[iCent][1];
      valRes1Temp[iCent][4] = valMixEp[iCent][4]*valMixEp[iCent][2]/valMixEp[iCent][0];
      valRes1Temp[iCent][5] = valMixEp[iCent][2]*valMixEp[iCent][3]/valMixEp[iCent][5];
      errRes1Temp[iCent][0] = propMixEpResErr(valMixEp[iCent][4],errMixEp[iCent][4],valMixEp[iCent][1],errMixEp[iCent][1],valMixEp[iCent][3],errMixEp[iCent][3]);
      errRes1Temp[iCent][1] = propMixEpResErr(valMixEp[iCent][4],errMixEp[iCent][4],valMixEp[iCent][0],errMixEp[iCent][0],valMixEp[iCent][2],errMixEp[iCent][2]);
      errRes1Temp[iCent][2] = propMixEpResErr(valMixEp[iCent][0],errMixEp[iCent][0],valMixEp[iCent][1],errMixEp[iCent][1],valMixEp[iCent][5],errMixEp[iCent][5]);
      errRes1Temp[iCent][3] = propMixEpResErr(valMixEp[iCent][4],errMixEp[iCent][4],valMixEp[iCent][3],errMixEp[iCent][3],valMixEp[iCent][1],errMixEp[iCent][1]);
      errRes1Temp[iCent][4] = propMixEpResErr(valMixEp[iCent][4],errMixEp[iCent][4],valMixEp[iCent][2],errMixEp[iCent][2],valMixEp[iCent][0],errMixEp[iCent][0]);
      errRes1Temp[iCent][5] = propMixEpResErr(valMixEp[iCent][2],errMixEp[iCent][2],valMixEp[iCent][3],errMixEp[iCent][3],valMixEp[iCent][5],errMixEp[iCent][5]);
    }
  }

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      mMixSubEp1Res1Val[iCent][iGrp] = 0.0;
      mMixSubEp1Res1Err[iCent][iGrp] = 0.0;
      mMixSubEp1Res2Val[iCent][iGrp] = 0.0;
      mMixSubEp1Res2Err[iCent][iGrp] = 0.0;
    }
  }

  TF1 *f_MixEp1Res1 = new TF1("f_MixEp1Res1",funcMixEp1Res1,0,10,0);
  TF1 *f_MixEp1Res2 = new TF1("f_MixEp1Res2",funcMixEp1Res2,0,10,0);
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      double valRes1Sub = -999.9;
      double errRes1Sub = 1.0;
      double valRes1Raw = valRes1Temp[iCent][iGrp];
      double errRes1Raw = errRes1Temp[iCent][iGrp];

      double valRes12Sub = -999.9;
      double errRes12Sub = 1.0;

      if(valRes1Raw > 0 && valRes1Raw < 1.0)
      {
	valRes1Sub = TMath::Sqrt(valRes1Raw);
	errRes1Sub = errRes1Raw/(2.0*valRes1Sub);

	double chi1Sub    = f_MixEp1Res1->GetX(valRes1Sub); // calculate 2nd EP Res
	valRes12Sub       = f_MixEp1Res2->Eval(chi1Sub);
	double errChi1Sub = errRes1Sub/f_MixEp1Res1->Derivative(chi1Sub); // error propagation
	errRes12Sub       = f_MixEp1Res2->Derivative(chi1Sub)*errChi1Sub;
      }

      mMixSubEp1Res1Val[iCent][iGrp] = valRes1Sub;
      mMixSubEp1Res1Err[iCent][iGrp] = errRes1Sub;
      mMixSubEp1Res2Val[iCent][iGrp] = valRes12Sub;
      mMixSubEp1Res2Err[iCent][iGrp] = errRes12Sub;
    }
  }
  file_mResolution->Close();
}

double StMixEpManager::propMixEpResErr(double valA, double sigA, double valB, double sigB, double valC, double sigC)
{ // return the error of valA*valB/valC
  double errA = sigA/valA;
  double errB = sigB/valB;
  double errC = sigC/valC;
  double valAxB = valA*valB;
  double sigAxB = valAxB*TMath::Sqrt(errA*errA+errB*errB); // valA*valB*sqrt((sigA/valA)^2+(sigB/valB)^2)
  double errAB = sigAxB/valAxB;

  double sigABdivC = valAxB/valC*TMath::Sqrt(errAB*errAB+errC*errC); // valAB/valC*sqrt((sigAB/valAB)^2+(sigC/valC)^2);

  return sigABdivC;
}

double StMixEpManager::getMixSubEp1Res1Val(int cent9, int grpId)
{
  return mMixSubEp1Res1Val[cent9][grpId];
}

double StMixEpManager::getMixSubEp1Res1Err(int cent9, int grpId)
{
  return mMixSubEp1Res1Err[cent9][grpId];
}

double StMixEpManager::getMixSubEp1Res2Val(int cent9, int grpId)
{
  return mMixSubEp1Res2Val[cent9][grpId];
}

double StMixEpManager::getMixSubEp1Res2Err(int cent9, int grpId)
{
  return mMixSubEp1Res2Err[cent9][grpId];
}
//---------------------------------------------------------------------------------
// deutron efficiency
void StMixEpManager::readDeuEfficiency()
{
  std::string inputFile = Form("Utility/EventPlaneMaker/%s/ChargedFlow/efficiency.root",globCons::str_mBeamType[mType].c_str());
  file_mDeuEfficiency = TFile::Open(inputFile.c_str());
  tracking_d = (TH2D*)file_mDeuEfficiency->Get("tracking_d");
  tpc_d = (TH1D*)file_mDeuEfficiency->Get("tpc_d");
  tof_d = (TH1D*)file_mDeuEfficiency->Get("tof_d");
  tofmatch = (TH2D*)file_mDeuEfficiency->Get("tofmatch");
}

float StMixEpManager::calcDeuEfficiency(float pT, float pMag, float etaLab, float yCms)
{
  const float dP_cut =  3.2;
  float d_eff = (float)tracking_d->GetBinContent(tracking_d->FindBin(yCms, pT));
  if(pMag < dP_cut) d_eff = d_eff*(float)tpc_d->GetBinContent(tpc_d->FindBin(pMag));
  if(pMag > dP_cut) d_eff = d_eff*(float)tof_d->GetBinContent(tof_d->FindBin(pMag));
  if(pMag > dP_cut) d_eff = d_eff*(float)tofmatch->GetBinContent(tofmatch->FindBin(etaLab, pT));
  if(d_eff>1) d_eff = 1.;

  return d_eff;
}

// deutron Directed Flow
void StMixEpManager::initMixSubEpFlow()
{
  p_mMixSubEpDeuV1Eff = new TProfile("p_mMixSubEpDeuV1Eff","p_mMixSubEpDeuV1Eff",20,-1.0,1.0);
  p_mMixSubEpDeuV1Eff->Sumw2();
  p_mMixSubEpDeuV1 = new TProfile("p_mMixSubEpDeuV1","p_mMixSubEpDeuV1",20,-1.0,1.0);
  p_mMixSubEpDeuV1->Sumw2();
}

void StMixEpManager::fillMixSubEpDeuV1(double yCms, double v1, double refWgt, double eff)
{
  p_mMixSubEpDeuV1Eff->Fill(yCms, v1, refWgt/eff);
  p_mMixSubEpDeuV1->Fill(yCms, v1, refWgt);
}

void StMixEpManager::writeMixSubEpFlow()
{
  p_mMixSubEpDeuV1Eff->Write();
  p_mMixSubEpDeuV1->Write();
}
//---------------------------------------------------------------------------------
// Event Plane Distribution
void StMixEpManager::initMixSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      std::string histName = Form("h_mMixEp1Grp%dRawCent%d",iGrp,iCent);
      h_mMixEp1RawCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StMixEpManager::fillMixSubEpRaw(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest)
{
  h_mMixEp1RawCorr[mCent9][0]->Fill(Psi1EpdGrp0, Psi1TpcEast); // 0: EpdEpGrp0 vs. TpcEpEast
  h_mMixEp1RawCorr[mCent9][1]->Fill(Psi1EpdGrp0, Psi1TpcWest); // 1: EpdEpGrp0 vs. TpcEpWest
  h_mMixEp1RawCorr[mCent9][2]->Fill(Psi1EpdGrp1, Psi1TpcEast); // 2: EpdEpGrp1 vs. TpcEpEast
  h_mMixEp1RawCorr[mCent9][3]->Fill(Psi1EpdGrp1, Psi1TpcWest); // 3: EpdEpGrp1 vs. TpcEpWest
  h_mMixEp1RawCorr[mCent9][4]->Fill(Psi1EpdGrp0, Psi1EpdGrp1); // 4: EpdEpGrp0 vs. EpdEpGrp1
  h_mMixEp1RawCorr[mCent9][5]->Fill(Psi1TpcEast, Psi1TpcWest); // 5: TpcEpEast vs. TpcEpWest
}

void StMixEpManager::writeMixSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      h_mMixEp1RawCorr[iCent][iGrp]->Write();
    }
  }
}

void StMixEpManager::initMixSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      std::string histName = Form("h_mMixEp1Grp%dReCtrCent%d",iGrp,iCent);
      h_mMixEp1ReCtrCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StMixEpManager::fillMixSubEpReCtr(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest)
{
  h_mMixEp1ReCtrCorr[mCent9][0]->Fill(Psi1EpdGrp0, Psi1TpcEast); // 0: EpdEpGrp0 vs. TpcEpEast
  h_mMixEp1ReCtrCorr[mCent9][1]->Fill(Psi1EpdGrp0, Psi1TpcWest); // 1: EpdEpGrp0 vs. TpcEpWest
  h_mMixEp1ReCtrCorr[mCent9][2]->Fill(Psi1EpdGrp1, Psi1TpcEast); // 2: EpdEpGrp1 vs. TpcEpEast
  h_mMixEp1ReCtrCorr[mCent9][3]->Fill(Psi1EpdGrp1, Psi1TpcWest); // 3: EpdEpGrp1 vs. TpcEpWest
  h_mMixEp1ReCtrCorr[mCent9][4]->Fill(Psi1EpdGrp0, Psi1EpdGrp1); // 4: EpdEpGrp0 vs. EpdEpGrp1
  h_mMixEp1ReCtrCorr[mCent9][5]->Fill(Psi1TpcEast, Psi1TpcWest); // 5: TpcEpEast vs. TpcEpWest
}

void StMixEpManager::writeMixSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      h_mMixEp1ReCtrCorr[iCent][iGrp]->Write();
    }
  }
}

void StMixEpManager::initMixSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      std::string histName = Form("h_mMixEp1Grp%dShiftCent%d",iGrp,iCent);
      h_mMixEp1ShiftCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StMixEpManager::fillMixSubEpShift(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest)
{
  h_mMixEp1ShiftCorr[mCent9][0]->Fill(Psi1EpdGrp0, Psi1TpcEast); // 0: EpdEpGrp0 vs. TpcEpEast
  h_mMixEp1ShiftCorr[mCent9][1]->Fill(Psi1EpdGrp0, Psi1TpcWest); // 1: EpdEpGrp0 vs. TpcEpWest
  h_mMixEp1ShiftCorr[mCent9][2]->Fill(Psi1EpdGrp1, Psi1TpcEast); // 2: EpdEpGrp1 vs. TpcEpEast
  h_mMixEp1ShiftCorr[mCent9][3]->Fill(Psi1EpdGrp1, Psi1TpcWest); // 3: EpdEpGrp1 vs. TpcEpWest
  h_mMixEp1ShiftCorr[mCent9][4]->Fill(Psi1EpdGrp0, Psi1EpdGrp1); // 4: EpdEpGrp0 vs. EpdEpGrp1
  h_mMixEp1ShiftCorr[mCent9][5]->Fill(Psi1TpcEast, Psi1TpcWest); // 5: TpcEpEast vs. TpcEpWest
}

void StMixEpManager::writeMixSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      h_mMixEp1ShiftCorr[iCent][iGrp]->Write();
    }
  }
}
//---------------------------------------------------------------------------------
