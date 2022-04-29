#include "StRoot/StRunQAMaker/StRunQAHistoManager.h"
#include "StRoot/StRunQAMaker/StRunQACons.h"

#include <TH2F.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>

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
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    std::string HistName = Form("h_mTriggerID_%s",mCutsQA[i_cut].c_str());
    h_mTriggerID[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),10,-0.5,9.5);

    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      HistName = Form("h_mRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mRefMult[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5);

      HistName = Form("h_mGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mGRefMult[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5);

      HistName = Form("h_mRefMultGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mRefMultGRefMult[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5,2000,-0.5,1999.5);

      HistName = Form("h_mCentrality9%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mCentrality9[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),12,-1.5,10.5);

      HistName = Form("h_mTofMatchRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofMatchRefMult[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),4000,-0.5,3999.5,2000,-0.5,1999.5);

      HistName = Form("h_mTofHitsRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofHitsRefMult[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),4000,-0.5,3999.5,2000,-0.5,1999.5);

      HistName = Form("h_mTofMatchGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofMatchGRefMult[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),4000,-0.5,3999.5,2000,-0.5,1999.5);

      HistName = Form("h_mTofHitsGRefMult%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mTofHitsGRefMult[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),4000,-0.5,3999.5,2000,-0.5,1999.5);

      HistName = Form("h_mVertexXY%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mVertexXY[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05,201,-10.05,10.05);

      HistName = Form("h_mVertexZ%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mVertexZ[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),301,-150.5,150.5);

      HistName = Form("h_mVzVzVpd%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mVzVzVpd[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),301,-150.5,150.5,301,-150.5,150.5);

      HistName = Form("h_mDiffVzVzVpd%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mDiffVzVzVpd[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05);
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
  h_mRefMult[cutSelection][9]->Fill(refMult);
  h_mGRefMult[cutSelection][9]->Fill(grefMult);
  h_mRefMultGRefMult[cutSelection][9]->Fill(refMult,grefMult);
  h_mCentrality9[cutSelection][9]->Fill(cent9,reweight);
  h_mTofMatchRefMult[cutSelection][9]->Fill(tofMatch,refMult);
  h_mTofHitsRefMult[cutSelection][9]->Fill(tofHits,refMult);
  h_mTofMatchGRefMult[cutSelection][9]->Fill(tofMatch,grefMult);
  h_mTofHitsGRefMult[cutSelection][9]->Fill(tofHits,grefMult);
}

void StRunQAHistoManager::fillEventQA_Vertex(int triggerBin, float vx, float vy, float vz, float vzVpd, int cutSelection)
{
  // for a specific triggerBin
  h_mVertexXY[cutSelection][triggerBin]->Fill(vx,vy);
  h_mVertexZ[cutSelection][triggerBin]->Fill(vz);
  h_mVzVzVpd[cutSelection][triggerBin]->Fill(vz,vzVpd);
  h_mDiffVzVzVpd[cutSelection][triggerBin]->Fill(vz-vzVpd);

  // for all triggers
  h_mVertexXY[cutSelection][9]->Fill(vx,vy);
  h_mVertexZ[cutSelection][9]->Fill(vz);
  h_mVzVzVpd[cutSelection][9]->Fill(vz,vzVpd);
  h_mDiffVzVzVpd[cutSelection][9]->Fill(vz-vzVpd);
}

void StRunQAHistoManager::fillEventQA_Trigger(int triggerBin, int cutSelection)
{
  h_mTriggerID[cutSelection]->Fill(triggerBin);
}

void StRunQAHistoManager::writeEventQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    h_mTriggerID[i_cut]->Write();
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      h_mRefMult[i_cut][i_trig]->Write();
      h_mGRefMult[i_cut][i_trig]->Write();
      h_mRefMultGRefMult[i_cut][i_trig]->Write();
      h_mCentrality9[i_cut][i_trig]->Write();
      h_mTofMatchRefMult[i_cut][i_trig]->Write();
      h_mTofHitsRefMult[i_cut][i_trig]->Write();
      h_mTofMatchGRefMult[i_cut][i_trig]->Write();
      h_mTofHitsGRefMult[i_cut][i_trig]->Write();

      h_mVertexXY[i_cut][i_trig]->Write();
      h_mVertexZ[i_cut][i_trig]->Write();
      h_mVzVzVpd[i_cut][i_trig]->Write();
      h_mDiffVzVzVpd[i_cut][i_trig]->Write();
    }
  }
}
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
//Track QA
void StRunQAHistoManager::initTrackQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      std::string HistName = Form("h_mPrimPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mPrimPt[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.5,49.5);

      HistName = Form("h_mPrimEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mPrimEta[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),205,-2.05,2.05);

      HistName = Form("h_mPrimPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mPrimPhi[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),720,-TMath::TwoPi(),TMath::TwoPi());

      HistName = Form("h_mGlobPt%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mGlobPt[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.5,49.5);

      HistName = Form("h_mGlobEta%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mGlobEta[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),205,-2.05,2.05);

      HistName = Form("h_mGlobPhi%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mGlobPhi[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),720,-TMath::TwoPi(),TMath::TwoPi());

      HistName = Form("h_mDca%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mDca[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.1,4.9);

      HistName = Form("h_mNHitsFit%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mNHitsFit[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),50,-0.5,49.5);

      HistName = Form("h_mNHitsRatio%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mNHitsRatio[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),120,-0.05,1.15);

      HistName = Form("h_mNHitsDEdx%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mNHitsDEdx[i_cut][i_trig] = new TH1F(HistName.c_str(),HistName.c_str(),50,-0.5,49.5);

      HistName = Form("h_mDEdxMom%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mDEdxMom[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),450,-4.5,4.5,400,0,40);

      HistName = Form("h_mMass2Mom%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mMass2Mom[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),450,-4.5,4.5,200,-0.3,1.7);

      HistName = Form("h_mBetaMom%s_trigger%d",mCutsQA[i_cut].c_str(),i_trig);
      h_mBetaMom[i_cut][i_trig] = new TH2F(HistName.c_str(),HistName.c_str(),450,-4.5,4.5,300,0.0,3.0);
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
  h_mPrimPt[cutSelection][9]->Fill(pMom.Pt()); 
  h_mPrimEta[cutSelection][9]->Fill(pMom.Eta()); 
  h_mPrimPhi[cutSelection][9]->Fill(pMom.Phi()); 
  h_mGlobPt[cutSelection][9]->Fill(gMom.Pt()); 
  h_mGlobEta[cutSelection][9]->Fill(gMom.Eta()); 
  h_mGlobPhi[cutSelection][9]->Fill(gMom.Phi()); 
}

void StRunQAHistoManager::fillTrackQA_Quliaty(int triggerBin, float gDca, int nHitsFit, int nHitsMax, int nHitsDEdx, int cutSelection)
{
  // for a specific triggerBin
  h_mDca[cutSelection][triggerBin]->Fill(gDca); 
  h_mNHitsDEdx[cutSelection][triggerBin]->Fill(nHitsDEdx); 
  h_mNHitsFit[cutSelection][triggerBin]->Fill(nHitsFit); 

  // for all triggers
  h_mDca[cutSelection][9]->Fill(gDca); 
  h_mNHitsDEdx[cutSelection][9]->Fill(nHitsDEdx); 
  h_mNHitsFit[cutSelection][9]->Fill(nHitsFit); 
  if(nHitsMax > 0)
  {
    float nHitsRatio = (float)nHitsFit/(float)nHitsMax;
    h_mNHitsRatio[cutSelection][triggerBin]->Fill(nHitsRatio); // for a specific triggerBin
    h_mNHitsRatio[cutSelection][9]->Fill(nHitsRatio); // for all triggers
  }
  else
  {
    h_mNHitsRatio[cutSelection][triggerBin]->Fill(-0.01); // for a specific triggerBin
    h_mNHitsRatio[cutSelection][9]->Fill(-0.01); // for all triggers
  }
}

void StRunQAHistoManager::fillTrackQA_PID(int triggerBin, float mom, short charge, float dEdx, float beta, float mass2, int cutSelection)
{
  // for a specific triggerBin
  h_mDEdxMom[cutSelection][triggerBin]->Fill(mom*charge, dEdx); 
  h_mMass2Mom[cutSelection][triggerBin]->Fill(mom*charge, mass2); 
  h_mBetaMom[cutSelection][triggerBin]->Fill(mom*charge, 1.0/beta); 

  // for all triggers
  h_mDEdxMom[cutSelection][9]->Fill(mom*charge, dEdx); 
  h_mMass2Mom[cutSelection][9]->Fill(mom*charge, mass2); 
  h_mBetaMom[cutSelection][9]->Fill(mom*charge, 1.0/beta); 
}

void StRunQAHistoManager::writeTrackQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    for(int i_trig = 0; i_trig < 10; ++i_trig)
    {
      h_mPrimPt[i_cut][i_trig]->Write(); 
      h_mPrimEta[i_cut][i_trig]->Write(); 
      h_mPrimPhi[i_cut][i_trig]->Write(); 
      h_mGlobPt[i_cut][i_trig]->Write(); 
      h_mGlobEta[i_cut][i_trig]->Write(); 
      h_mGlobPhi[i_cut][i_trig]->Write(); 
      h_mDca[i_cut][i_trig]->Write(); 
      h_mNHitsFit[i_cut][i_trig]->Write(); 
      h_mNHitsRatio[i_cut][i_trig]->Write(); 
      h_mNHitsDEdx[i_cut][i_trig]->Write(); 
      h_mDEdxMom[i_cut][i_trig]->Write(); 
      h_mMass2Mom[i_cut][i_trig]->Write(); 
      h_mBetaMom[i_cut][i_trig]->Write(); 
    }
  }
}
//-------------------------------------------------------------------------------------------
