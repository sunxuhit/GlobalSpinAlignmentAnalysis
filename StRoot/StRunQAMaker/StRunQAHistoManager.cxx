#include <iostream>

#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TString.h"

#include "StRoot/StRunQAMaker/StRunQAHistoManager.h"

ClassImp(StRunQAHistoManager)

//-------------------------------------------------------------------------------------------

StRunQAHistoManager::StRunQAHistoManager()
{
}

//-------------------------------------------------------------------------------------------

StRunQAHistoManager::~StRunQAHistoManager()
{
  /* */
}

//-------------------------------------------------------------------------------------------
//Event QA
void StRunQAHistoManager::initEventQA()
{
  for(int iCut = 0; iCut < mNumCuts; ++iCut)
  {
    std::string HistName = Form("h_mTriggerId%s",mCutStatus[iCut].c_str());
    h_mTriggerId[iCut] = new TH1F(HistName.c_str(),HistName.c_str(),10,-0.5,9.5);

    for(int iTrig = 0; iTrig < mNumTriggerBins; ++iTrig)
    {
      HistName = Form("h_mRefMult%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mRefMult[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5);

      HistName = Form("h_mGRefMult%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mGRefMult[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5);

      HistName = Form("h_mRefMultGRefMult%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mRefMultGRefMult[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5,2000,-0.5,1999.5);

      HistName = Form("h_mCentrality9%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mCentrality9[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),12,-1.5,10.5);

      HistName = Form("h_mTofMatchRefMult%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mTofMatchRefMult[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),4000,-0.5,3999.5,2000,-0.5,1999.5);

      HistName = Form("h_mTofHitsRefMult%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mTofHitsRefMult[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),4000,-0.5,3999.5,2000,-0.5,1999.5);

      HistName = Form("h_mTofMatchGRefMult%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mTofMatchGRefMult[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),4000,-0.5,3999.5,2000,-0.5,1999.5);

      HistName = Form("h_mTofHitsGRefMult%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mTofHitsGRefMult[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),4000,-0.5,3999.5,2000,-0.5,1999.5);

      HistName = Form("h_mVertexXY%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mVertexXY[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05,201,-10.05,10.05);

      HistName = Form("h_mVertexZ%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mVertexZ[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),301,-150.5,150.5);

      HistName = Form("h_mVzVzVpd%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mVzVzVpd[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),301,-150.5,150.5,301,-150.5,150.5);

      HistName = Form("h_mDiffVzVzVpd%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mDiffVzVzVpd[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05);
    }
  }
}

void StRunQAHistoManager::fillEventQA_RefMult(int triggerBin, int refMult, int grefMult, int cent9, double reweight, int tofHits, int tofMatch, int cutSelection)
{
  // for a specific triggerBin
  h_mRefMult[cutSelection][triggerBin]->Fill(refMult);
  h_mGRefMult[cutSelection][triggerBin]->Fill(grefMult);
  h_mRefMultGRefMult[cutSelection][triggerBin]->Fill(refMult,grefMult);
  h_mCentrality9[cutSelection][triggerBin]->Fill(cent9,reweight);
  h_mTofMatchRefMult[cutSelection][triggerBin]->Fill(tofMatch,refMult);
  h_mTofHitsRefMult[cutSelection][triggerBin]->Fill(tofHits,refMult);
  h_mTofMatchGRefMult[cutSelection][triggerBin]->Fill(tofMatch,grefMult);
  h_mTofHitsGRefMult[cutSelection][triggerBin]->Fill(tofHits,grefMult);

  // for all triggers
  h_mRefMult[cutSelection][mNumTriggerBins-1]->Fill(refMult);
  h_mGRefMult[cutSelection][mNumTriggerBins-1]->Fill(grefMult);
  h_mRefMultGRefMult[cutSelection][mNumTriggerBins-1]->Fill(refMult,grefMult);
  h_mCentrality9[cutSelection][mNumTriggerBins-1]->Fill(cent9,reweight);
  h_mTofMatchRefMult[cutSelection][mNumTriggerBins-1]->Fill(tofMatch,refMult);
  h_mTofHitsRefMult[cutSelection][mNumTriggerBins-1]->Fill(tofHits,refMult);
  h_mTofMatchGRefMult[cutSelection][mNumTriggerBins-1]->Fill(tofMatch,grefMult);
  h_mTofHitsGRefMult[cutSelection][mNumTriggerBins-1]->Fill(tofHits,grefMult);
}

void StRunQAHistoManager::fillEventQA_Vertex(int triggerBin, float vx, float vy, float vz, float vzVpd, int cutSelection)
{
  // for a specific triggerBin
  h_mVertexXY[cutSelection][triggerBin]->Fill(vx,vy);
  h_mVertexZ[cutSelection][triggerBin]->Fill(vz);
  h_mVzVzVpd[cutSelection][triggerBin]->Fill(vz,vzVpd);
  h_mDiffVzVzVpd[cutSelection][triggerBin]->Fill(vz-vzVpd);

  // for all triggers
  h_mVertexXY[cutSelection][mNumTriggerBins-1]->Fill(vx,vy);
  h_mVertexZ[cutSelection][mNumTriggerBins-1]->Fill(vz);
  h_mVzVzVpd[cutSelection][mNumTriggerBins-1]->Fill(vz,vzVpd);
  h_mDiffVzVzVpd[cutSelection][mNumTriggerBins-1]->Fill(vz-vzVpd);
}

void StRunQAHistoManager::fillEventQA_Trigger(int triggerBin, int cutSelection)
{
  h_mTriggerId[cutSelection]->Fill(triggerBin);
}

void StRunQAHistoManager::writeEventQA()
{
  for(int iCut = 0; iCut < mNumCuts; ++iCut)
  {
    h_mTriggerId[iCut]->Write();
    for(int iTrig = 0; iTrig < mNumTriggerBins; ++iTrig)
    {
      h_mRefMult[iCut][iTrig]->Write();
      h_mGRefMult[iCut][iTrig]->Write();
      h_mRefMultGRefMult[iCut][iTrig]->Write();
      h_mCentrality9[iCut][iTrig]->Write();
      h_mTofMatchRefMult[iCut][iTrig]->Write();
      h_mTofHitsRefMult[iCut][iTrig]->Write();
      h_mTofMatchGRefMult[iCut][iTrig]->Write();
      h_mTofHitsGRefMult[iCut][iTrig]->Write();

      h_mVertexXY[iCut][iTrig]->Write();
      h_mVertexZ[iCut][iTrig]->Write();
      h_mVzVzVpd[iCut][iTrig]->Write();
      h_mDiffVzVzVpd[iCut][iTrig]->Write();
    }
  }
}
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
//Track QA
void StRunQAHistoManager::initTrackQA()
{
  for(int iCut = 0; iCut < mNumCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < mNumTriggerBins; ++iTrig)
    {
      std::string HistName = Form("h_mPrimPt%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mPrimPt[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.5,49.5);

      HistName = Form("h_mPrimEta%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mPrimEta[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05);

      HistName = Form("h_mPrimPhi%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mPrimPhi[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),720,-TMath::TwoPi(),TMath::TwoPi());

      HistName = Form("h_mGlobPt%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mGlobPt[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.5,49.5);

      HistName = Form("h_mGlobEta%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mGlobEta[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05);

      HistName = Form("h_mGlobPhi%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mGlobPhi[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),720,-TMath::TwoPi(),TMath::TwoPi());

      HistName = Form("h_mDca%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mDca[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.1,9.9);

      HistName = Form("h_mNHitsFit%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mNHitsFit[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),100,-0.5,99.5);

      HistName = Form("h_mNHitsRatio%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mNHitsRatio[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),120,-0.05,1.15);

      HistName = Form("h_mNHitsDEdx%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mNHitsDEdx[iCut][iTrig] = new TH1F(HistName.c_str(),HistName.c_str(),100,-0.5,99.5);

      HistName = Form("h_mDEdxMom%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mMomDEdx[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),450,-4.5,4.5,400,0,40);

      HistName = Form("h_mMass2Mom%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mMomMass2[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),450,-4.5,4.5,200,-0.3,1.7);

      HistName = Form("h_mBetaMom%sTrigger%d",mCutStatus[iCut].c_str(),iTrig);
      h_mMomBeta[iCut][iTrig] = new TH2F(HistName.c_str(),HistName.c_str(),450,-4.5,4.5,300,0.0,3.0);
    }
  }
}

void StRunQAHistoManager::fillTrackQA_Kinematics(int triggerBin, TVector3 pMom, TVector3 gMom, int cutSelection)
{
  // for a specific triggerBin
  h_mPrimPt[cutSelection][triggerBin]->Fill(pMom.Pt()); 
  h_mPrimEta[cutSelection][triggerBin]->Fill(pMom.Eta()); 
  h_mPrimPhi[cutSelection][triggerBin]->Fill(pMom.Phi()); 
  h_mGlobPt[cutSelection][triggerBin]->Fill(gMom.Pt()); 
  h_mGlobEta[cutSelection][triggerBin]->Fill(gMom.Eta()); 
  h_mGlobPhi[cutSelection][triggerBin]->Fill(gMom.Phi()); 

  // for all triggers
  h_mPrimPt[cutSelection][mNumTriggerBins-1]->Fill(pMom.Pt()); 
  h_mPrimEta[cutSelection][mNumTriggerBins-1]->Fill(pMom.Eta()); 
  h_mPrimPhi[cutSelection][mNumTriggerBins-1]->Fill(pMom.Phi()); 
  h_mGlobPt[cutSelection][mNumTriggerBins-1]->Fill(gMom.Pt()); 
  h_mGlobEta[cutSelection][mNumTriggerBins-1]->Fill(gMom.Eta()); 
  h_mGlobPhi[cutSelection][mNumTriggerBins-1]->Fill(gMom.Phi()); 
}

void StRunQAHistoManager::fillTrackQA_Quliaty(int triggerBin, float gDca, int nHitsFit, int nHitsMax, int nHitsDEdx, int cutSelection)
{
  // for a specific triggerBin
  h_mDca[cutSelection][triggerBin]->Fill(gDca); 
  h_mNHitsDEdx[cutSelection][triggerBin]->Fill(nHitsDEdx); 
  h_mNHitsFit[cutSelection][triggerBin]->Fill(nHitsFit); 

  // for all triggers
  h_mDca[cutSelection][mNumTriggerBins-1]->Fill(gDca); 
  h_mNHitsDEdx[cutSelection][mNumTriggerBins-1]->Fill(nHitsDEdx); 
  h_mNHitsFit[cutSelection][mNumTriggerBins-1]->Fill(nHitsFit); 
  if(nHitsMax > 0)
  {
    float nHitsRatio = (float)nHitsFit/(float)nHitsMax;
    h_mNHitsRatio[cutSelection][triggerBin]->Fill(nHitsRatio); // for a specific triggerBin
    h_mNHitsRatio[cutSelection][mNumTriggerBins-1]->Fill(nHitsRatio); // for all triggers
  }
  else
  {
    h_mNHitsRatio[cutSelection][triggerBin]->Fill(-0.01); // for a specific triggerBin
    h_mNHitsRatio[cutSelection][mNumTriggerBins-1]->Fill(-0.01); // for all triggers
  }
}

void StRunQAHistoManager::fillTrackQA_PID(int triggerBin, float mom, short charge, float dEdx, float beta, float mass2, int cutSelection)
{
  // for a specific triggerBin
  h_mMomDEdx[cutSelection][triggerBin]->Fill(mom*charge, dEdx); 
  h_mMomMass2[cutSelection][triggerBin]->Fill(mom*charge, mass2); 
  h_mMomBeta[cutSelection][triggerBin]->Fill(mom*charge, 1.0/beta); 

  // for all triggers
  h_mMomDEdx[cutSelection][mNumTriggerBins-1]->Fill(mom*charge, dEdx); 
  h_mMomMass2[cutSelection][mNumTriggerBins-1]->Fill(mom*charge, mass2); 
  h_mMomBeta[cutSelection][mNumTriggerBins-1]->Fill(mom*charge, 1.0/beta); 
}

void StRunQAHistoManager::writeTrackQA()
{
  for(int iCut = 0; iCut < mNumCuts; ++iCut)
  {
    for(int iTrig = 0; iTrig < mNumTriggerBins; ++iTrig)
    {
      h_mPrimPt[iCut][iTrig]->Write(); 
      h_mPrimEta[iCut][iTrig]->Write(); 
      h_mPrimPhi[iCut][iTrig]->Write(); 
      h_mGlobPt[iCut][iTrig]->Write(); 
      h_mGlobEta[iCut][iTrig]->Write(); 
      h_mGlobPhi[iCut][iTrig]->Write(); 
      h_mDca[iCut][iTrig]->Write(); 
      h_mNHitsFit[iCut][iTrig]->Write(); 
      h_mNHitsRatio[iCut][iTrig]->Write(); 
      h_mNHitsDEdx[iCut][iTrig]->Write(); 
      h_mMomDEdx[iCut][iTrig]->Write(); 
      h_mMomMass2[iCut][iTrig]->Write(); 
      h_mMomBeta[iCut][iTrig]->Write(); 
    }
  }
}
//-------------------------------------------------------------------------------------------
