#include <string>
#include <iostream>

#include <TProfile.h>
#include <TMath.h>
#include <TString.h>

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StRunQAMaker/StRunQAProManager.h"

ClassImp(StRunQAProManager)

//---------------------------------------------------------------------------------

StRunQAProManager::StRunQAProManager(int beamType) : mType(beamType)
{
}

//---------------------------------------------------------------------------------

StRunQAProManager::~StRunQAProManager()
{
  /* */
}

//---------------------------------------------------------------------------------

void StRunQAProManager::initRunQA()
{
  for(int iCut = 0; iCut < mNumCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < mNumTriggerBins; ++iTrig)
    {
      std::string ProName = Form("p_mRefMult%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mRefMult[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mGRefMult%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mGRefMult[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mZdcX%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mZdcX[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mVz%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mVz[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mVr%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mVr[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mGDca%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mGDca[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mNHitsFit%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mNHitsFit[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mPrimPt%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mPrimPt[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mPrimEta%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mPrimEta[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mPrimPhi%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mPrimPhi[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mGlobPt%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mGlobPt[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mGlobEta%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mGlobEta[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);

      ProName = Form("p_mGlobPhi%sTrigger%d",str_mCutStatus[iCut].c_str(),iTrig);
      p_mGlobPhi[iCut][iTrig] = new TProfile(ProName.c_str(),ProName.c_str(),globCons::mNumRunIndex[mType],(double)globCons::mRunIndexLo[mType]-0.5,(double)globCons::mRunIndexHi[mType]-0.5);
    }
  }
}

void StRunQAProManager::fillRunQAEvt(int triggerBin, int runIndex, double refMult, double grefMult, double zdcX, double vx, double vy, double vz, int cutSelection)
{
  // for a specific triggerBin
  p_mRefMult[cutSelection][triggerBin]->Fill(runIndex, refMult);
  p_mGRefMult[cutSelection][triggerBin]->Fill(runIndex, grefMult);
  p_mZdcX[cutSelection][triggerBin]->Fill(runIndex, zdcX);
  p_mVz[cutSelection][triggerBin]->Fill(runIndex, vz);
  p_mVr[cutSelection][triggerBin]->Fill(runIndex, TMath::Sqrt(vx*vx+vy*vy));

  // for all triggers
  p_mRefMult[cutSelection][mNumTriggerBins-1]->Fill(runIndex, refMult);
  p_mGRefMult[cutSelection][mNumTriggerBins-1]->Fill(runIndex, grefMult);
  p_mZdcX[cutSelection][mNumTriggerBins-1]->Fill(runIndex, zdcX);
  p_mVz[cutSelection][mNumTriggerBins-1]->Fill(runIndex, vz);
  p_mVr[cutSelection][mNumTriggerBins-1]->Fill(runIndex, TMath::Sqrt(vx*vx+vy*vy));
}

void StRunQAProManager::fillRunQATrk(int triggerBin, int runIndex, double gDca, int nHitsFit, TVector3 pMom, TVector3 gMom, int cutSelection)
{
  // for a specific triggerBin
  p_mGDca[cutSelection][triggerBin]->Fill(runIndex, gDca);
  p_mNHitsFit[cutSelection][triggerBin]->Fill(runIndex, nHitsFit);
  p_mPrimPt[cutSelection][triggerBin]->Fill(runIndex, pMom.Pt());
  p_mPrimEta[cutSelection][triggerBin]->Fill(runIndex, pMom.Eta());
  p_mPrimPhi[cutSelection][triggerBin]->Fill(runIndex, pMom.Phi());
  p_mGlobPt[cutSelection][triggerBin]->Fill(runIndex, gMom.Pt());
  p_mGlobEta[cutSelection][triggerBin]->Fill(runIndex, gMom.Eta());
  p_mGlobPhi[cutSelection][triggerBin]->Fill(runIndex, gMom.Phi());

  // for all triggers
  p_mGDca[cutSelection][mNumTriggerBins-1]->Fill(runIndex, gDca);
  p_mNHitsFit[cutSelection][mNumTriggerBins-1]->Fill(runIndex, nHitsFit);
  p_mPrimPt[cutSelection][mNumTriggerBins-1]->Fill(runIndex, pMom.Pt());
  p_mPrimEta[cutSelection][mNumTriggerBins-1]->Fill(runIndex, pMom.Eta());
  p_mPrimPhi[cutSelection][mNumTriggerBins-1]->Fill(runIndex, pMom.Phi());
  p_mGlobPt[cutSelection][mNumTriggerBins-1]->Fill(runIndex, gMom.Pt());
  p_mGlobEta[cutSelection][mNumTriggerBins-1]->Fill(runIndex, gMom.Eta());
  p_mGlobPhi[cutSelection][mNumTriggerBins-1]->Fill(runIndex, gMom.Phi());
}

void StRunQAProManager::writeRunQA()
{
  for(int iCut = 0; iCut < mNumCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < mNumTriggerBins; ++iTrig)
    {
      p_mRefMult[iCut][iTrig]->Write();
      p_mGRefMult[iCut][iTrig]->Write();
      p_mZdcX[iCut][iTrig]->Write();
      p_mVz[iCut][iTrig]->Write();
      p_mVr[iCut][iTrig]->Write();
      p_mGDca[iCut][iTrig]->Write();
      p_mNHitsFit[iCut][iTrig]->Write();
      p_mPrimPt[iCut][iTrig]->Write();
      p_mPrimEta[iCut][iTrig]->Write();
      p_mPrimPhi[iCut][iTrig]->Write();
      p_mGlobPt[iCut][iTrig]->Write();
      p_mGlobEta[iCut][iTrig]->Write();
      p_mGlobPhi[iCut][iTrig]->Write();
    }
  }
}
