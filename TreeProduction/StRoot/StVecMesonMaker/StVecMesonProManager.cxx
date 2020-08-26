#include "StRoot/StVecMesonMaker/StVecMesonProManager.h"
#include <TProfile.h>
#include <TMath.h>
#include <TString.h>

#include <string>

ClassImp(StVecMesonProManager)

//---------------------------------------------------------------------------------

StVecMesonProManager::StVecMesonProManager()
{
}

//---------------------------------------------------------------------------------

StVecMesonProManager::~StVecMesonProManager()
{
  /* */
}

//---------------------------------------------------------------------------------

void StVecMesonProManager::initRunQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    std::string ProName = Form("p_mQA_RefMult_%s",mCutsQA[i_cut].c_str());
    p_mQA_RefMult[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5); 

    ProName = Form("p_mQA_ZdcX_%s",mCutsQA[i_cut].c_str());
    p_mQA_ZdcX[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_Vz_%s",mCutsQA[i_cut].c_str());
    p_mQA_Vz[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_Vr_%s",mCutsQA[i_cut].c_str());
    p_mQA_Vr[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_gDca_%s",mCutsQA[i_cut].c_str());
    p_mQA_gDca[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_nHitsFit_%s",mCutsQA[i_cut].c_str());
    p_mQA_nHitsFit[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_PrimPt_%s",mCutsQA[i_cut].c_str());
    p_mQA_PrimPt[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_PrimEta_%s",mCutsQA[i_cut].c_str());
    p_mQA_PrimEta[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_PrimPhi_%s",mCutsQA[i_cut].c_str());
    p_mQA_PrimPhi[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_GlobPt_%s",mCutsQA[i_cut].c_str());
    p_mQA_GlobPt[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_GlobEta_%s",mCutsQA[i_cut].c_str());
    p_mQA_GlobEta[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);

    ProName = Form("p_mQA_GlobPhi_%s",mCutsQA[i_cut].c_str());
    p_mQA_GlobPhi[i_cut] = new TProfile(ProName.c_str(),ProName.c_str(),2000,-0.5,1999.5);
  }
}

void StVecMesonProManager::fillRunQA_Event(int runIdenx, float refMult, float zdcX, float vx, float vy, float vz, int cutSelection)
{
    p_mQA_RefMult[cutSelection]->Fill(runIdenx, refMult);
    p_mQA_ZdcX[cutSelection]->Fill(runIdenx, zdcX);
    p_mQA_Vz[cutSelection]->Fill(runIdenx, vz);
    p_mQA_Vr[cutSelection]->Fill(runIdenx, TMath::Sqrt(vx*vx+vy*vy));
}

void StVecMesonProManager::fillRunQA_Track(int runIdenx, float gDca, int nHitsFit, TVector3 pMom, TVector3 gMom, int cutSelection)
{
    p_mQA_gDca[cutSelection]->Fill(runIdenx, gDca);
    p_mQA_nHitsFit[cutSelection]->Fill(runIdenx, nHitsFit);
    p_mQA_PrimPt[cutSelection]->Fill(runIdenx, pMom.Pt());
    p_mQA_PrimEta[cutSelection]->Fill(runIdenx, pMom.Eta());
    p_mQA_PrimPhi[cutSelection]->Fill(runIdenx, pMom.Phi());
    p_mQA_GlobPt[cutSelection]->Fill(runIdenx, gMom.Pt());
    p_mQA_GlobEta[cutSelection]->Fill(runIdenx, gMom.Eta());
    p_mQA_GlobPhi[cutSelection]->Fill(runIdenx, gMom.Phi());
}

void StVecMesonProManager::writeRunQA()
{
  for(int i_cut = 0; i_cut < 2; ++i_cut)
  {
    p_mQA_RefMult[i_cut]->Write();
    p_mQA_ZdcX[i_cut]->Write();
    p_mQA_Vz[i_cut]->Write();
    p_mQA_Vr[i_cut]->Write();
    p_mQA_gDca[i_cut]->Write();
    p_mQA_nHitsFit[i_cut]->Write();
    p_mQA_PrimPt[i_cut]->Write();
    p_mQA_PrimEta[i_cut]->Write();
    p_mQA_PrimPhi[i_cut]->Write();
    p_mQA_GlobPt[i_cut]->Write();
    p_mQA_GlobEta[i_cut]->Write();
    p_mQA_GlobPhi[i_cut]->Write();
  }
}
