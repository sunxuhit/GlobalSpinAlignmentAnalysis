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
    std::string HistName = Form("h_mRefMult_%s",mCutsQA[i_cut].c_str());
    h_mRefMult[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5);

    HistName = Form("h_mGRefMult_%s",mCutsQA[i_cut].c_str());
    h_mGRefMult[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5);

    HistName = Form("h_mCentrality9_%s",mCutsQA[i_cut].c_str());
    h_mCentrality9[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),12,-1.5,10.5);

    HistName = Form("h_mRefMultTofMatch_%s",mCutsQA[i_cut].c_str());
    h_mRefMultTofMatch[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5,2000,-0.5,1999.5);

    HistName = Form("h_mRefMultTofHits_%s",mCutsQA[i_cut].c_str());
    h_mRefMultTofHits[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5,2000,-0.5,1999.5);

    HistName = Form("h_mGRefMultTofMatch_%s",mCutsQA[i_cut].c_str());
    h_mGRefMultTofMatch[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5,2000,-0.5,1999.5);

    HistName = Form("h_mGRefMultTofHits_%s",mCutsQA[i_cut].c_str());
    h_mGRefMultTofHits[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),2000,-0.5,1999.5,2000,-0.5,1999.5);

    HistName = Form("h_mVertexXY_%s",mCutsQA[i_cut].c_str());
    h_mVertexXY[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05,201,-10.05,10.05);

    HistName = Form("h_mVertexZ_%s",mCutsQA[i_cut].c_str());
    h_mVertexZ[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),301,-150.5,150.5);

    HistName = Form("h_mVzVzVpd_%s",mCutsQA[i_cut].c_str());
    h_mVzVzVpd[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),301,-150.5,150.5,301,-150.5,150.5);

    HistName = Form("h_mDiffVzVzVpd_%s",mCutsQA[i_cut].c_str());
    h_mDiffVzVzVpd[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),201,-10.05,10.05);
  }
}

void StRunQAHistoManager::fillEventQA_RefMult(int refMult, int grefMult, int cent9, double reweight, int tofHits, int tofMatch, int cutSelection)
{
  h_mRefMult[cutSelection]->Fill(refMult);
  h_mGRefMult[cutSelection]->Fill(grefMult);
  h_mCentrality9[cutSelection]->Fill(cent9,reweight);
  h_mRefMultTofMatch[cutSelection]->Fill(refMult,tofMatch);
  h_mRefMultTofHits[cutSelection]->Fill(refMult,tofHits);
  h_mGRefMultTofMatch[cutSelection]->Fill(grefMult,tofMatch);
  h_mGRefMultTofHits[cutSelection]->Fill(grefMult,tofHits);
}

void StRunQAHistoManager::fillEventQA_Vertex(float vx, float vy, float vz, float vzVpd, int cutSelection)
{
  h_mVertexXY[cutSelection]->Fill(vx,vy);
  h_mVertexZ[cutSelection]->Fill(vz);
  h_mVzVzVpd[cutSelection]->Fill(vz,vzVpd);
  h_mDiffVzVzVpd[cutSelection]->Fill(vz-vzVpd);
}

void StRunQAHistoManager::writeEventQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    h_mRefMult[i_cut]->Write();
    h_mGRefMult[i_cut]->Write();
    h_mCentrality9[i_cut]->Write();
    h_mRefMultTofMatch[i_cut]->Write();
    h_mRefMultTofHits[i_cut]->Write();
    h_mGRefMultTofMatch[i_cut]->Write();
    h_mGRefMultTofHits[i_cut]->Write();

    h_mVertexXY[i_cut]->Write();
    h_mVertexZ[i_cut]->Write();
    h_mVzVzVpd[i_cut]->Write();
    h_mDiffVzVzVpd[i_cut]->Write();
  }
}
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
//Track QA
void StRunQAHistoManager::initTrackQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    std::string HistName = Form("h_mPrimPt_%s",mCutsQA[i_cut].c_str());
    h_mPrimPt[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.5,49.5);

    HistName = Form("h_mPrimEta_%s",mCutsQA[i_cut].c_str());
    h_mPrimEta[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),205,-2.05,2.05);

    HistName = Form("h_mPrimPhi_%s",mCutsQA[i_cut].c_str());
    h_mPrimPhi[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),720,-TMath::TwoPi(),TMath::TwoPi());

    HistName = Form("h_mGlobPt_%s",mCutsQA[i_cut].c_str());
    h_mGlobPt[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.5,49.5);

    HistName = Form("h_mGlobEta_%s",mCutsQA[i_cut].c_str());
    h_mGlobEta[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),205,-2.05,2.05);

    HistName = Form("h_mGlobPhi_%s",mCutsQA[i_cut].c_str());
    h_mGlobPhi[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),720,-TMath::TwoPi(),TMath::TwoPi());

    HistName = Form("h_mDca_%s",mCutsQA[i_cut].c_str());
    h_mDca[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),1000,-0.1,4.9);

    HistName = Form("h_mNHitsFit_%s",mCutsQA[i_cut].c_str());
    h_mNHitsFit[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),50,-0.5,49.5);

    HistName = Form("h_mNHitsRatio_%s",mCutsQA[i_cut].c_str());
    h_mNHitsRatio[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),120,-0.05,1.15);

    HistName = Form("h_mNHitsDEdx_%s",mCutsQA[i_cut].c_str());
    h_mNHitsDEdx[i_cut] = new TH1F(HistName.c_str(),HistName.c_str(),50,-0.5,49.5);

    HistName = Form("h_mDEdxMom_%s",mCutsQA[i_cut].c_str());
    h_mDEdxMom[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),450,-4.5,4.5,400,0,40);

    HistName = Form("h_mMass2Mom_%s",mCutsQA[i_cut].c_str());
    h_mMass2Mom[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),450,-4.5,4.5,200,-0.3,1.7);

    HistName = Form("h_mBetaMom_%s",mCutsQA[i_cut].c_str());
    h_mBetaMom[i_cut] = new TH2F(HistName.c_str(),HistName.c_str(),450,-4.5,4.5,300,0.0,3.0);
  }
}

void StRunQAHistoManager::fillTrackQA_Kinematics(TVector3 pMom, TVector3 gMom, int cutSelection)
{
  h_mPrimPt[cutSelection]->Fill(pMom.Pt()); 
  h_mPrimEta[cutSelection]->Fill(pMom.Eta()); 
  h_mPrimPhi[cutSelection]->Fill(pMom.Phi()); 

  h_mGlobPt[cutSelection]->Fill(gMom.Pt()); 
  h_mGlobEta[cutSelection]->Fill(gMom.Eta()); 
  h_mGlobPhi[cutSelection]->Fill(gMom.Phi()); 
}

void StRunQAHistoManager::fillTrackQA_Quliaty(float gDca, int nHitsFit, int nHitsMax, int nHitsDEdx, int cutSelection)
{
  h_mDca[cutSelection]->Fill(gDca); 
  h_mNHitsFit[cutSelection]->Fill(nHitsFit); 
  if(nHitsMax > 0)
  {
    float nHitsRatio = (float)nHitsFit/(float)nHitsMax;
    h_mNHitsRatio[cutSelection]->Fill(nHitsRatio); 
  }
  else
  {
    h_mNHitsRatio[cutSelection]->Fill(-0.01); 
  }
  h_mNHitsDEdx[cutSelection]->Fill(nHitsDEdx); 
}

void StRunQAHistoManager::fillTrackQA_PID(float mom, short charge, float dEdx, float beta, float mass2, int cutSelection)
{
  h_mDEdxMom[cutSelection]->Fill(mom*charge, dEdx); 
  h_mMass2Mom[cutSelection]->Fill(mom*charge, mass2); 
  h_mBetaMom[cutSelection]->Fill(mom*charge, 1.0/beta); 
}

void StRunQAHistoManager::writeTrackQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    h_mPrimPt[i_cut]->Write(); 
    h_mPrimEta[i_cut]->Write(); 
    h_mPrimPhi[i_cut]->Write(); 
    h_mGlobPt[i_cut]->Write(); 
    h_mGlobEta[i_cut]->Write(); 
    h_mGlobPhi[i_cut]->Write(); 
    h_mDca[i_cut]->Write(); 
    h_mNHitsFit[i_cut]->Write(); 
    h_mNHitsRatio[i_cut]->Write(); 
    h_mNHitsDEdx[i_cut]->Write(); 
    h_mDEdxMom[i_cut]->Write(); 
    h_mMass2Mom[i_cut]->Write(); 
    h_mBetaMom[i_cut]->Write(); 
  }
}
//-------------------------------------------------------------------------------------------
