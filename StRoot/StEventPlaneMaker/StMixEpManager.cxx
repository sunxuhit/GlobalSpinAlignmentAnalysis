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
  mCent9 = -1;
}

void StMixEpManager::initMixEpManager(int cent9, int runIndex, int vzBin)
{
  mCent9 = cent9;
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
  double subRes1Grp0 = TMath::Cos(Psi1EpdGrp0-Psi1EpdGrp1)*TMath::Cos(Psi1EpdGrp0-Psi1TpcWest)/TMath::Cos(Psi1EpdGrp1-Psi1TpcWest);
  double subRes1Grp1 = TMath::Cos(Psi1EpdGrp0-Psi1EpdGrp1)*TMath::Cos(Psi1EpdGrp0-Psi1TpcEast)/TMath::Cos(Psi1EpdGrp1-Psi1TpcEast);
  double subRes1Grp2 = TMath::Cos(Psi1EpdGrp0-Psi1TpcEast)*TMath::Cos(Psi1EpdGrp0-Psi1TpcWest)/TMath::Cos(Psi1TpcEast-Psi1TpcWest);
  p_mMixSubEp1Res[0]->Fill((double)mCent9,subRes1Grp0); // 0: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpWest => default
  p_mMixSubEp1Res[1]->Fill((double)mCent9,subRes1Grp1); // 1: EpdEpGrp0 vs. EpdEpGrp1 && TpcEpEast
  p_mMixSubEp1Res[2]->Fill((double)mCent9,subRes1Grp2); // 2: EpdEpGrp0 vs. TpcEpEast && TpcEpWest

  double subRes1Grp3 = TMath::Cos(Psi1EpdGrp1-Psi1EpdGrp0)*TMath::Cos(Psi1EpdGrp1-Psi1TpcWest)/TMath::Cos(Psi1EpdGrp0-Psi1TpcWest);
  double subRes1Grp4 = TMath::Cos(Psi1EpdGrp1-Psi1EpdGrp0)*TMath::Cos(Psi1EpdGrp1-Psi1TpcEast)/TMath::Cos(Psi1EpdGrp0-Psi1TpcEast);
  double subRes1Grp5 = TMath::Cos(Psi1EpdGrp1-Psi1TpcEast)*TMath::Cos(Psi1EpdGrp1-Psi1TpcWest)/TMath::Cos(Psi1TpcEast-Psi1TpcWest);
  p_mMixSubEp1Res[3]->Fill((double)mCent9,subRes1Grp3); // 3: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpWest => main Systematic
  p_mMixSubEp1Res[4]->Fill((double)mCent9,subRes1Grp4); // 4: EpdEpGrp1 vs. EpdEpGrp0 && TpcEpEast
  p_mMixSubEp1Res[5]->Fill((double)mCent9,subRes1Grp5); // 5: EpdEpGrp1 vs. TpcEpEast && TpcEpWest
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

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      mMixSubEp1ResVal[iCent][iGrp]  = 0.0;
      mMixSubEp1ResErr[iCent][iGrp]  = 0.0;
    }
  }

  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      double valRes1Sub  = -999.9;
      double errRes1Sub  = 1.0;
      double valRes1Raw  = p_mMixSubEp1Res[iGrp]->GetBinContent(iCent+1);
      double errRes1Raw  = p_mMixSubEp1Res[iGrp]->GetBinError(iCent+1);
      if(valRes1Raw > 0)
      {
	valRes1Sub = TMath::Sqrt(valRes1Raw);
	errRes1Sub = errRes1Raw/(2.0*valRes1Sub);
      }
      mMixSubEp1ResVal[iCent][iGrp]  = valRes1Sub;
      mMixSubEp1ResErr[iCent][iGrp]  = errRes1Sub;
    }
  }

  file_mResolution->Close();
}

double StMixEpManager::getMixSubEp1ResVal(int cent9, int grpId)
{
  return mMixSubEp1ResVal[cent9][grpId];
}

double StMixEpManager::getMixSubEp1ResErr(int cent9, int grpId)
{
  return mMixSubEp1ResErr[cent9][grpId];
}
//---------------------------------------------------------------------------------
// deutron Directed Flow
void StMixEpManager::initMixSubEpFlow()
{
  p_mMixSubEpDeuV1 = new TProfile("p_mMixSubEpDeuV1","p_mMixSubEpDeuV1",20,-1.0,1.0);
  p_mMixSubEpDeuV1->Sumw2();
}

void StMixEpManager::fillMixSubEpDeuV1(double rap, double v1, double reweight)
{
  p_mMixSubEpDeuV1->Fill(rap, v1, reweight);
}

void StMixEpManager::writeMixSubEpFlow()
{
  p_mMixSubEpDeuV1->Write();
}
//---------------------------------------------------------------------------------
// Event Plane Distribution
void StEpdEpManager::initMixSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1Grp0RawEastCent%d",iCent); // EPD 1st EP
    h_mEpdEp1Grp0RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1Grp1RawEastCent%d",iCent);
    h_mEpdEp1Grp1RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());

    histName = Form("h_mTpcEp1RawEastCent%d",iCent); // TPC 1st EP
    h_mTpcEp1RawEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1RawWestCent%d",iCent);
    h_mTpcEp1RawWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      histName = Form("h_mMixEp1Grp%dRawCent%d",iGrp,iCent);
      h_mMixEp1RawCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StEpdEpManager::fillMixSubEpRaw(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest)
{
    h_mEpdEp1Grp0RawEast[mCent9]->Fill(mRunIndex, Psi1EpdGrp0);
    h_mEpdEp1Grp1RawEast[mCent9]->Fill(mRunIndex, Psi1EpdGrp1);
    h_mTpcEp1RawEast[mCent9]->Fill(mRunIndex, Psi1TpcEast);
    h_mTpcEp1RawWest[mCent9]->Fill(mRunIndex, Psi1TpcWest);

    h_mMixEp1RawCorr[mCent9][0]->Fill(Psi1EpdGrp0, Psi1TpcEast); // 0: EpdEpGrp0 vs. TpcEpEast
    h_mMixEp1RawCorr[mCent9][1]->Fill(Psi1EpdGrp0, Psi1TpcWest); // 1: EpdEpGrp0 vs. TpcEpWest
    h_mMixEp1RawCorr[mCent9][2]->Fill(Psi1EpdGrp1, Psi1TpcEast); // 2: EpdEpGrp1 vs. TpcEpEast
    h_mMixEp1RawCorr[mCent9][3]->Fill(Psi1EpdGrp1, Psi1TpcWest); // 3: EpdEpGrp1 vs. TpcEpWest
    h_mMixEp1RawCorr[mCent9][4]->Fill(Psi1EpdGrp0, Psi1EpdGrp1); // 4: EpdEpGrp0 vs. EpdEpGrp1
    h_mMixEp1RawCorr[mCent9][5]->Fill(Psi1TpcEast, Psi1TpcWest); // 5: TpcEpEast vs. TpcEpWest
  }
}

void StEpdEpManager::writeMixSubEpRaw()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1Grp0RawEast[iCent]->Write();
    h_mEpdEp1Grp1RawEast[iCent]->Write();
    h_mTpcEp1RawEast[iCent]->Write();
    h_mTpcEp1RawWest[iCent]->Write();
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      h_mMixEp1RawCorr[iCent][iGrp]->Write();
    }
  }
}

void StEpdEpManager::initMixSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1Grp0ReCtrEastCent%d",iCent); // EPD 1st EP
    h_mEpdEp1Grp0ReCtrEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1Grp1ReCtrEastCent%d",iCent);
    h_mEpdEp1Grp1ReCtrEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());

    histName = Form("h_mTpcEp1ReCtrEastCent%d",iCent); // TPC 1st EP
    h_mTpcEp1ReCtrEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1ReCtrWestCent%d",iCent);
    h_mTpcEp1ReCtrWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      histName = Form("h_mMixEp1Grp%dReCtrCent%d",iGrp,iCent);
      h_mMixEp1ReCtrCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StEpdEpManager::fillMixSubEpReCtr(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest)
{
    h_mEpdEp1Grp0ReCtrEast[mCent9]->Fill(mRunIndex, Psi1EpdGrp0);
    h_mEpdEp1Grp1ReCtrEast[mCent9]->Fill(mRunIndex, Psi1EpdGrp1);
    h_mTpcEp1ReCtrEast[mCent9]->Fill(mRunIndex, Psi1TpcEast);
    h_mTpcEp1ReCtrWest[mCent9]->Fill(mRunIndex, Psi1TpcWest);

    h_mMixEp1ReCtrCorr[mCent9][0]->Fill(Psi1EpdGrp0, Psi1TpcEast); // 0: EpdEpGrp0 vs. TpcEpEast
    h_mMixEp1ReCtrCorr[mCent9][1]->Fill(Psi1EpdGrp0, Psi1TpcWest); // 1: EpdEpGrp0 vs. TpcEpWest
    h_mMixEp1ReCtrCorr[mCent9][2]->Fill(Psi1EpdGrp1, Psi1TpcEast); // 2: EpdEpGrp1 vs. TpcEpEast
    h_mMixEp1ReCtrCorr[mCent9][3]->Fill(Psi1EpdGrp1, Psi1TpcWest); // 3: EpdEpGrp1 vs. TpcEpWest
    h_mMixEp1ReCtrCorr[mCent9][4]->Fill(Psi1EpdGrp0, Psi1EpdGrp1); // 4: EpdEpGrp0 vs. EpdEpGrp1
    h_mMixEp1ReCtrCorr[mCent9][5]->Fill(Psi1TpcEast, Psi1TpcWest); // 5: TpcEpEast vs. TpcEpWest
  }
}

void StEpdEpManager::writeMixSubEpReCtr()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1Grp0ReCtrEast[iCent]->Write();
    h_mEpdEp1Grp1ReCtrEast[iCent]->Write();
    h_mTpcEp1ReCtrEast[iCent]->Write();
    h_mTpcEp1ReCtrWest[iCent]->Write();
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      h_mMixEp1ReCtrCorr[iCent][iGrp]->Write();
    }
  }
}

void StEpdEpManager::initMixSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mEpdEp1Grp0ShiftEastCent%d",iCent); // EPD 1st EP
    h_mEpdEp1Grp0ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mEpdEp1Grp1ShiftEastCent%d",iCent);
    h_mEpdEp1Grp1ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());

    histName = Form("h_mTpcEp1ShiftEastCent%d",iCent); // TPC 1st EP
    h_mTpcEp1ShiftEast[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    histName = Form("h_mTpcEp1ShiftWestCent%d",iCent);
    h_mTpcEp1ShiftWest[iCent] = new TH2F(histName.c_str(),histName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5,540,-1.5*TMath::Pi(),1.5*TMath::Pi());
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      histName = Form("h_mMixEp1Grp%dShiftCent%d",iGrp,iCent);
      h_mMixEp1ShiftCorr[iCent][iGrp] = new TH2F(histName.c_str(),histName.c_str(),135,-1.5*TMath::Pi(),1.5*TMath::Pi(),135,-1.5*TMath::Pi(),1.5*TMath::Pi());
    }
  }
}

void StEpdEpManager::fillMixSubEpShift(double Psi1EpdGrp0, double Psi1EpdGrp1, double Psi1TpcEast, double Psi1TpcWest)
{
    h_mEpdEp1Grp0ShiftEast[mCent9]->Fill(mRunIndex, Psi1EpdGrp0);
    h_mEpdEp1Grp1ShiftEast[mCent9]->Fill(mRunIndex, Psi1EpdGrp1);
    h_mTpcEp1ShiftEast[mCent9]->Fill(mRunIndex, Psi1TpcEast);
    h_mTpcEp1ShiftWest[mCent9]->Fill(mRunIndex, Psi1TpcWest);

    h_mMixEp1ShiftCorr[mCent9][0]->Fill(Psi1EpdGrp0, Psi1TpcEast); // 0: EpdEpGrp0 vs. TpcEpEast
    h_mMixEp1ShiftCorr[mCent9][1]->Fill(Psi1EpdGrp0, Psi1TpcWest); // 1: EpdEpGrp0 vs. TpcEpWest
    h_mMixEp1ShiftCorr[mCent9][2]->Fill(Psi1EpdGrp1, Psi1TpcEast); // 2: EpdEpGrp1 vs. TpcEpEast
    h_mMixEp1ShiftCorr[mCent9][3]->Fill(Psi1EpdGrp1, Psi1TpcWest); // 3: EpdEpGrp1 vs. TpcEpWest
    h_mMixEp1ShiftCorr[mCent9][4]->Fill(Psi1EpdGrp0, Psi1EpdGrp1); // 4: EpdEpGrp0 vs. EpdEpGrp1
    h_mMixEp1ShiftCorr[mCent9][5]->Fill(Psi1TpcEast, Psi1TpcWest); // 5: TpcEpEast vs. TpcEpWest
  }
}

void StEpdEpManager::writeMixSubEpShift()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mEpdEp1Grp0ShiftEast[iCent]->Write();
    h_mEpdEp1Grp1ShiftEast[iCent]->Write();
    h_mTpcEp1ShiftEast[iCent]->Write();
    h_mTpcEp1ShiftWest[iCent]->Write();
    for(int iGrp = 0; iGrp < mNumEpGroup; ++iGrp)
    {
      h_mMixEp1ShiftCorr[iCent][iGrp]->Write();
    }
  }
}
//---------------------------------------------------------------------------------
