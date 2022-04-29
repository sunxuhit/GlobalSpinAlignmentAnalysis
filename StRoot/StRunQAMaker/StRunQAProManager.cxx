#include "StRoot/StRunQAMaker/StRunQAProManager.h"
#include "StRoot/StRunQAMaker/StRunQACons.h"
#include <TProfile.h>
#include <TMath.h>
#include <TString.h>

#include <string>

ClassImp(StRunQAProManager)

//---------------------------------------------------------------------------------

StRunQAProManager::StRunQAProManager()
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
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      std::string ProName = Form("p_mRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mRefMult[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5); 

      ProName = Form("p_mGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGRefMult[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5); 

      ProName = Form("p_mZdcX%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mZdcX[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mVz%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mVz[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mVr%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mVr[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mGDca%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGDca[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mNHitsFit%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mNHitsFit[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mPrimPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimPt[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mPrimEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimEta[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mPrimPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mPrimPhi[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mGlobPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobPt[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mGlobEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobEta[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);

      ProName = Form("p_mGlobPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      p_mGlobPhi[i_cut][i_trig] = new TProfile(ProName.c_str(),ProName.c_str(),runQA::mNumOfRunIndex,-0.5,(float)runQA::mNumOfRunIndex-0.5);
    }
  }
}

void StRunQAProManager::fillRunQA_Event(int triggerBin, int runIdenx, float refMult, float grefMult, float zdcX, float vx, float vy, float vz, int cutSelection)
{
  // for a specific triggerBin
  p_mRefMult[cutSelection][triggerBin]->Fill(runIdenx, refMult);
  p_mGRefMult[cutSelection][triggerBin]->Fill(runIdenx, grefMult);
  p_mZdcX[cutSelection][triggerBin]->Fill(runIdenx, zdcX);
  p_mVz[cutSelection][triggerBin]->Fill(runIdenx, vz);
  p_mVr[cutSelection][triggerBin]->Fill(runIdenx, TMath::Sqrt(vx*vx+vy*vy));

  // for all triggers
  p_mRefMult[cutSelection][9]->Fill(runIdenx, refMult);
  p_mGRefMult[cutSelection][9]->Fill(runIdenx, grefMult);
  p_mZdcX[cutSelection][9]->Fill(runIdenx, zdcX);
  p_mVz[cutSelection][9]->Fill(runIdenx, vz);
  p_mVr[cutSelection][9]->Fill(runIdenx, TMath::Sqrt(vx*vx+vy*vy));
}

void StRunQAProManager::fillRunQA_Track(int triggerBin, int runIdenx, float gDca, int nHitsFit, TVector3 pMom, TVector3 gMom, int cutSelection)
{
  // for a specific triggerBin
  p_mGDca[cutSelection][triggerBin]->Fill(runIdenx, gDca);
  p_mNHitsFit[cutSelection][triggerBin]->Fill(runIdenx, nHitsFit);
  p_mPrimPt[cutSelection][triggerBin]->Fill(runIdenx, pMom.Pt());
  p_mPrimEta[cutSelection][triggerBin]->Fill(runIdenx, pMom.Eta());
  p_mPrimPhi[cutSelection][triggerBin]->Fill(runIdenx, pMom.Phi());
  p_mGlobPt[cutSelection][triggerBin]->Fill(runIdenx, gMom.Pt());
  p_mGlobEta[cutSelection][triggerBin]->Fill(runIdenx, gMom.Eta());
  p_mGlobPhi[cutSelection][triggerBin]->Fill(runIdenx, gMom.Phi());

  // for all triggers
  p_mGDca[cutSelection][9]->Fill(runIdenx, gDca);
  p_mNHitsFit[cutSelection][9]->Fill(runIdenx, nHitsFit);
  p_mPrimPt[cutSelection][9]->Fill(runIdenx, pMom.Pt());
  p_mPrimEta[cutSelection][9]->Fill(runIdenx, pMom.Eta());
  p_mPrimPhi[cutSelection][9]->Fill(runIdenx, pMom.Phi());
  p_mGlobPt[cutSelection][9]->Fill(runIdenx, gMom.Pt());
  p_mGlobEta[cutSelection][9]->Fill(runIdenx, gMom.Eta());
  p_mGlobPhi[cutSelection][9]->Fill(runIdenx, gMom.Phi());
}

void StRunQAProManager::writeRunQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      p_mRefMult[i_cut][i_trig]->Write();
      p_mGRefMult[i_cut][i_trig]->Write();
      p_mZdcX[i_cut][i_trig]->Write();
      p_mVz[i_cut][i_trig]->Write();
      p_mVr[i_cut][i_trig]->Write();
      p_mGDca[i_cut][i_trig]->Write();
      p_mNHitsFit[i_cut][i_trig]->Write();
      p_mPrimPt[i_cut][i_trig]->Write();
      p_mPrimEta[i_cut][i_trig]->Write();
      p_mPrimPhi[i_cut][i_trig]->Write();
      p_mGlobPt[i_cut][i_trig]->Write();
      p_mGlobEta[i_cut][i_trig]->Write();
      p_mGlobPhi[i_cut][i_trig]->Write();
    }
  }
}
