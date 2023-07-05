#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"

#include "StMessMgr.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StPhiMesonAnalyzer/StPhiMesonHistoManger.h"

ClassImp(StPhiMesonHistoManger)

//-------------------------------------------------------------
StPhiMesonHistoManger::StPhiMesonHistoManger(const int beamType, const int flagME) : mType(beamType), mFlagME(flagME)
{
}

StPhiMesonHistoManger::~StPhiMesonHistoManger()
{
}
//-------------------------------------------------------------
void StPhiMesonHistoManger::initPhiQA()
{
  for(int iCent = 0; iCent < mNumCentBinQA; ++iCent)
  {
    std::string histName = Form("h_mInvMassPhiQA%sCent%d",str_mMixEvt[mFlagME].c_str(),iCent);
    h_mInvMassPhiQA[iCent] = new TH2F(histName.c_str(),histName.c_str(),mNumPtBinQA,0.0,5.0,anaUtils::mNumInvMassPhi,anaUtils::mMassPhiMin,anaUtils::mMassPhiMax);
    for(int iPt = 0; iPt < mNumPtBinQA; ++iPt)
    {
      for(int iRap = 0; iRap < mNumRapBinQA; ++iRap)
      {
	std::string acptPhiLabKey = Form("h_mAcptPhi%sLabCent%dPt%dRap%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iRap);
	h_mAcptPhiLab[acptPhiLabKey] = new TH1F(acptPhiLabKey.c_str(),acptPhiLabKey.c_str(),anaUtils::mNumInvMassPhi,anaUtils::mMassPhiMin,anaUtils::mMassPhiMax);
	std::string acptPhiCmsKey = Form("h_mAcptPhi%sCmsCent%dPt%dRap%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iRap);
	h_mAcptPhiCms[acptPhiCmsKey] = new TH1F(acptPhiCmsKey.c_str(),acptPhiCmsKey.c_str(),anaUtils::mNumInvMassPhi,anaUtils::mMassPhiMin,anaUtils::mMassPhiMax);
      }
    }
  }
}

void StPhiMesonHistoManger::fillPhiQA(int cent9, double pt, double yLab, double yCms, double invMass, double refWgt)
{
  int ptBin = getPtBinPhiQA(pt);
  int yLabBin = getRapBinPhiQA(yLab);
  int yCmsBin = getRapBinPhiQA(yCms);
  h_mInvMassPhiQA[cent9]->Fill(pt,invMass,refWgt);
  if(is2060(cent9)) h_mInvMassPhiQA[9]->Fill(pt,invMass,refWgt);

  if(ptBin >= 0 && yLabBin >= 0 && yCmsBin >= 0)
  {
    std::string acptPhiLabKey = Form("h_mAcptPhi%sLabCent%dPt%dRap%d",str_mMixEvt[mFlagME].c_str(),cent9,ptBin,yLabBin);
    h_mAcptPhiLab[acptPhiLabKey]->Fill(invMass,refWgt);
    std::string acptPhiCmsKey = Form("h_mAcptPhi%sCmsCent%dPt%dRap%d",str_mMixEvt[mFlagME].c_str(),cent9,ptBin,yCmsBin);
    h_mAcptPhiCms[acptPhiCmsKey]->Fill(invMass,refWgt);

    if( is2060(cent9) )
    { // fill 20-60%
      acptPhiLabKey = Form("h_mAcptPhi%sLabCent9Pt%dRap%d",str_mMixEvt[mFlagME].c_str(),ptBin,yLabBin);
      h_mAcptPhiLab[acptPhiLabKey]->Fill(invMass,refWgt);
      acptPhiCmsKey = Form("h_mAcptPhi%sCmsCent9Pt%dRap%d",str_mMixEvt[mFlagME].c_str(),ptBin,yCmsBin);
      h_mAcptPhiCms[acptPhiCmsKey]->Fill(invMass,refWgt);
    }
  }
}

void StPhiMesonHistoManger::writePhiQA()
{
  for(int iCent = 0; iCent < mNumCentBinQA; ++iCent)
  {
    h_mInvMassPhiQA[iCent]->Write();
    for(int iPt = 0; iPt < mNumPtBinQA; ++iPt)
    {
      for(int iRap = 0; iRap < mNumRapBinQA; ++iRap)
      {
	std::string acptPhiLabKey = Form("h_mAcptPhi%sLabCent%dPt%dRap%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iRap);
	h_mAcptPhiLab[acptPhiLabKey]->Write();
	std::string acptPhiCmsKey = Form("h_mAcptPhi%sCmsCent%dPt%dRap%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iRap);
	h_mAcptPhiCms[acptPhiCmsKey]->Write();
      }
    }
  }
}

int StPhiMesonHistoManger::getPtBinPhiQA(double pt)
{
  int ptBin = -1;

  const double ptMin      = 0.0;
  const double ptMax      = 5.0;
  const double ptBinWidth = (ptMax-ptMin)/(double)mNumPtBinQA; // 0.2

  if(std::abs(pt-ptMin) < std::numeric_limits<double>::epsilon()) ptBin = 0;
  for(int iPt = 0; iPt < mNumPtBinQA; ++iPt)
  {
    if((pt > ptMin+iPt*ptBinWidth) && (pt <= ptMin+(iPt+1)*ptBinWidth))
    {
      ptBin = iPt;
    }
  }

  return ptBin;
}

int StPhiMesonHistoManger::getRapBinPhiQA(double y)
{
  int yBin = -1;

  const double yMin      = -1.25;
  const double yMax      = 1.25;
  const double yBinWidth = (yMax-yMin)/(double)mNumRapBinQA; // 0.1

  if(std::abs(y-yMin) < std::numeric_limits<double>::epsilon()) yBin = 0;
  for(int iPt = 0; iPt < mNumRapBinQA; ++iPt)
  {
    if((y > yMin+iPt*yBinWidth) && (y <= yMin+(iPt+1)*yBinWidth))
    {
      yBin = iPt;
    }
  }

  return yBin;
}
//-------------------------------------------------------------
bool StPhiMesonHistoManger::is0080(int cent9)
{
  if(cent9 < 0 || cent9 > 8) return false;

  return true;
}

bool StPhiMesonHistoManger::is2060(int cent9)
{
  if(cent9 < 2 || cent9 > 5) return false;

  return true;
}
//-------------------------------------------------------------
void StPhiMesonHistoManger::initIsoPhiFlow()
{
  for(int iCent = 0; iCent < mNumCentBinIsoPhiFlow; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinIsoPhiFlow; ++iPt)
    {
      for(int iPsi = 0; iPsi < mNumPsiBinIsoPhiFlow; ++iPsi)
      {
	std::string isoPhiV2Key = Form("h_mInvMassIsoPhiV2%sCent%dPt%dPsi%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iPsi);
	h_mInvMassIsoPhiV2[isoPhiV2Key] = new TH1F(isoPhiV2Key.c_str(),isoPhiV2Key.c_str(),anaUtils::mNumInvMassPhi,anaUtils::mMassPhiMin,anaUtils::mMassPhiMax);
	std::string isoPhiV3Key = Form("h_mInvMassIsoPhiV3%sCent%dPt%dPsi%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iPsi);
	h_mInvMassIsoPhiV3[isoPhiV3Key] = new TH1F(isoPhiV3Key.c_str(),isoPhiV3Key.c_str(),anaUtils::mNumInvMassPhi,anaUtils::mMassPhiMin,anaUtils::mMassPhiMax);
      }
    }
  }

  for(int iCent = 0; iCent < mNumCentBinIsoPhiYileds; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinIsoPhiFlow; ++iPt)
    {
      std::string isoPhiYieldsKey = Form("h_mInvMassIsoPhiYields%sCent%dPt%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt);
      h_mInvMassIsoPhiYields[isoPhiYieldsKey] = new TH1F(isoPhiYieldsKey.c_str(),isoPhiYieldsKey.c_str(),anaUtils::mNumInvMassPhi,anaUtils::mMassPhiMin,anaUtils::mMassPhiMax);
    }
  }
}

void StPhiMesonHistoManger::fillIsoPhiV2(int cent9, double pt, double yCms, double phi, double Psi2, double invMass, double res2, double refWgt)
{
  int centBin = getCentBinIsoPhiFlow(cent9);
  int ptBin   = getPtBinIsoPhiFlow(pt);
  int PsiBin  = getPsi2BinIsoPhiFlow(phi,Psi2);

  if(centBin > 0 && ptBin >= 0 && PsiBin >= 0)
  {
    std::string isoPhiV2CentKey = Form("h_mInvMassIsoPhiV2%sCent%dPt%dPsi%d",str_mMixEvt[mFlagME].c_str(),centBin,ptBin,PsiBin);
    h_mInvMassIsoPhiV2[isoPhiV2CentKey]->Fill(invMass,refWgt/res2); // centrality dependence

    std::string isoPhiV2MinBiasKey = Form("h_mInvMassIsoPhiV2%sCent0Pt%dPsi%d",str_mMixEvt[mFlagME].c_str(),ptBin,PsiBin);
    h_mInvMassIsoPhiV2[isoPhiV2MinBiasKey]->Fill(invMass,refWgt/res2); // minimum bias
  }
}

void StPhiMesonHistoManger::fillIsoPhiV3(int cent9, double pt, double yCms, double phi, double Psi3, double invMass, double res3, double refWgt)
{
  int centBin = getCentBinIsoPhiFlow(cent9);
  int ptBin   = getPtBinIsoPhiFlow(pt);
  int PsiBin  = getPsi3BinIsoPhiFlow(phi,Psi3);

  if(centBin > 0 && ptBin >= 0 && PsiBin >= 0)
  {
    std::string isoPhiV3CentKey = Form("h_mInvMassIsoPhiV3%sCent%dPt%dPsi%d",str_mMixEvt[mFlagME].c_str(),centBin,ptBin,PsiBin);
    h_mInvMassIsoPhiV3[isoPhiV3CentKey]->Fill(invMass,refWgt/res3); // centrality dependence

    std::string isoPhiV3MinBiasKey = Form("h_mInvMassIsoPhiV3%sCent0Pt%dPsi%d",str_mMixEvt[mFlagME].c_str(),ptBin,PsiBin);
    h_mInvMassIsoPhiV3[isoPhiV3MinBiasKey]->Fill(invMass,refWgt/res3); // minimum bias
  }
}

void StPhiMesonHistoManger::fillIsoPhiYields(int cent9, double pt, double yCms, double invMass, double refWgt)
{
  int ptBin = getPtBinIsoPhiFlow(pt);

  if(ptBin >= 0)
  {
    std::string isoPhiYieldsKey = Form("h_mInvMassIsoPhiYields%sCent%dPt%d",str_mMixEvt[mFlagME].c_str(),cent9,ptBin);
    h_mInvMassIsoPhiYields[isoPhiYieldsKey]->Fill(invMass,refWgt);
  }
}

void StPhiMesonHistoManger::writeIsoPhiFlow()
{
  for(int iCent = 0; iCent < mNumCentBinIsoPhiFlow; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinIsoPhiFlow; ++iPt)
    {
      for(int iPsi = 0; iPsi < mNumPsiBinIsoPhiFlow; ++iPsi)
      {
	std::string isoPhiV2Key = Form("h_mInvMassIsoPhiV2%sCent%dPt%dPsi%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iPsi);
	h_mInvMassIsoPhiV2[isoPhiV2Key]->Write();
	std::string isoPhiV3Key = Form("h_mInvMassIsoPhiV3%sCent%dPt%dPsi%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iPsi);
	h_mInvMassIsoPhiV3[isoPhiV3Key]->Write();
      }
    }
  }

  for(int iCent = 0; iCent < mNumCentBinIsoPhiYileds; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinIsoPhiFlow; ++iPt)
    {
      std::string isoPhiYieldsKey = Form("h_mInvMassIsoPhiYields%sCent%dPt%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt);
      h_mInvMassIsoPhiYields[isoPhiYieldsKey]->Write();
    }
  }
}

int StPhiMesonHistoManger::getCentBinIsoPhiFlow(int cent9)
{
  int centBin = -1;

  if(cent9 >= 7 && cent9 <= 8) return 1; //  0-10%
  if(cent9 >= 4 && cent9 <= 6) return 2; // 10-40%
  if(cent9 >= 0 && cent9 <= 3) return 3; // 40-80%

  return centBin;
}

int StPhiMesonHistoManger::getPtBinIsoPhiFlow(double pt)
{
  int ptBin = -1;

  const double ptLo[mNumPtBinIsoPhiFlow] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.4};
  const double ptHi[mNumPtBinIsoPhiFlow] = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.8,4.2,4.6,5.4,7.2};

  if(std::abs(pt-ptLo[0]) < std::numeric_limits<double>::epsilon()) ptBin = 0;
  for(int iPt = 0; iPt < mNumPtBinIsoPhiFlow; ++iPt)
  {
    if((pt > ptLo[iPt]) && (pt <= ptHi[iPt]))
    {
      ptBin = iPt;
    }
  }

  return ptBin;
}

int StPhiMesonHistoManger::getPsi2BinIsoPhiFlow(double phi, double Psi2)
{
  int PsiBin = -1;

  const double phiPsi2 = transPsi2(phi,Psi2); // [-pi/2,pi/2]
  if(!isPsi2InRange(phiPsi2)) return PsiBin;
  
  const double phiPsiMin      = 0.0;
  const double phiPsiMax      = TMath::Pi()/2.0;
  const double phiPsiBinWidth = (phiPsiMax-phiPsiMin)/(double)mNumPsiBinIsoPhiFlow; // pi/14

  if(std::abs(std::abs(phiPsi2)-phiPsiMin) < std::numeric_limits<double>::epsilon()) PsiBin = 0;
  for(int iPsi = 0; iPsi < mNumPsiBinIsoPhiFlow; ++iPsi)
  {
    if((std::abs(phiPsi2) > phiPsiMin+iPsi*phiPsiBinWidth) && (std::abs(phiPsi2) <= phiPsiMin+(iPsi+1)*phiPsiBinWidth))
    {
      PsiBin = iPsi;
    }
  }

  return PsiBin;
}

double StPhiMesonHistoManger::transPsi2(double phi, double Psi2)
{
  const double phiPsi2 = phi - Psi2; // [-3pi/2,3pi/2]
  double phiPsi2Corr = phiPsi2;
  if(phiPsi2 >  TMath::Pi()/2.0) phiPsi2Corr = phiPsi2 - TMath::Pi();
  if(phiPsi2 < -TMath::Pi()/2.0) phiPsi2Corr = phiPsi2 + TMath::Pi();

  return phiPsi2Corr;
}

bool StPhiMesonHistoManger::isPsi2InRange(double phiPsi2)
{
  if(phiPsi2 < -TMath::Pi()/2.0 || phiPsi2 > TMath::Pi()/2.0)
  {
    return false;
  }

  return true;
}

int StPhiMesonHistoManger::getPsi3BinIsoPhiFlow(double phi, double Psi3)
{
  int PsiBin = -1;

  const double phiPsi3 = transPsi3(phi,Psi3); // [-pi/2,pi/2]
  if(!isPsi3InRange(phiPsi3)) return PsiBin;
  
  const double phiPsiMin      = 0.0;
  const double phiPsiMax      = TMath::Pi()/3.0;
  const double phiPsiBinWidth = (phiPsiMax-phiPsiMin)/(double)mNumPsiBinIsoPhiFlow; // pi/21

  if(std::abs(std::abs(phiPsi3)-phiPsiMin) < std::numeric_limits<double>::epsilon()) PsiBin = 0;
  for(int iPsi = 0; iPsi < mNumPsiBinIsoPhiFlow; ++iPsi)
  {
    if((std::abs(phiPsi3) > phiPsiMin+iPsi*phiPsiBinWidth) && (std::abs(phiPsi3) <= phiPsiMin+(iPsi+1)*phiPsiBinWidth))
    {
      PsiBin = iPsi;
    }
  }

  return PsiBin;
}

double StPhiMesonHistoManger::transPsi3(double phi, double Psi3)
{
  const double phiPsi3 = phi - Psi3; // [-4pi/3,4pi/3]
  double phiPsi3Corr = phiPsi3;
  if(phiPsi3 >  TMath::Pi()) phiPsi3Corr = phiPsi3 - 2.0*TMath::TwoPi()/3.0;
  if(phiPsi3 < -TMath::Pi()) phiPsi3Corr = phiPsi3 + 2.0*TMath::TwoPi()/3.0;
  if(phiPsi3 >  TMath::Pi()/3.0 && phiPsi3 <=  TMath::Pi()) phiPsi3Corr = phiPsi3 - TMath::TwoPi()/3.0;
  if(phiPsi3 < -TMath::Pi()/3.0 && phiPsi3 >= -TMath::Pi()) phiPsi3Corr = phiPsi3 + TMath::TwoPi()/3.0;

  return phiPsi3Corr;
}

bool StPhiMesonHistoManger::isPsi3InRange(double phiPsi3)
{
  if(phiPsi3 < -TMath::Pi()/3.0 || phiPsi3 > TMath::Pi()/3.0)
  {
    return false;
  }

  return true;
}
//-------------------------------------------------------------
void StPhiMesonHistoManger::initFxtPhiFlow()
{
  for(int iCent = 0; iCent < mNumCentBinFxtPhiFlow; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinFxtPhiFlow; ++iPt)
    {
      for(int iPsi = 0; iPsi < mNumPsiBinFxtPhiFlow; ++iPsi)
      {
	std::string fxtPhiV2Key = Form("h_mInvMassFxtPhiV2%sCent%dPt%dPsi%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iPsi);
	h_mInvMassFxtPhiV2[fxtPhiV2Key] = new TH1F(fxtPhiV2Key.c_str(),fxtPhiV2Key.c_str(),anaUtils::mNumInvMassPhi,anaUtils::mMassPhiMin,anaUtils::mMassPhiMax);
      }
    }
  }

  for(int iCent = 0; iCent < mNumCentBinFxtPhiYileds; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinFxtPhiFlow; ++iPt)
    {
      std::string fxtPhiYieldsKey = Form("h_mInvMassFxtPhiYields%sCent%dPt%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt);
      h_mInvMassFxtPhiYields[fxtPhiYieldsKey] = new TH1F(fxtPhiYieldsKey.c_str(),fxtPhiYieldsKey.c_str(),anaUtils::mNumInvMassPhi,anaUtils::mMassPhiMin,anaUtils::mMassPhiMax);
    }
  }
}

void StPhiMesonHistoManger::fillFxtPhiV2(int cent9, double pt, double yCms, double phi, double Psi1, double invMass, double res12, double refWgt)
{
  int centBin = getCentBinFxtPhiFlow(cent9);
  int ptBin   = getPtBinFxtPhiFlow(pt);
  int PsiBin  = getPsi12BinFxtPhiFlow(phi,Psi1);

  if(centBin > 0 && ptBin >= 0 && PsiBin >= 0)
  {
    std::string fxtPhiV2CentKey = Form("h_mInvMassFxtPhiV2%sCent%dPt%dPsi%d",str_mMixEvt[mFlagME].c_str(),centBin,ptBin,PsiBin);
    h_mInvMassFxtPhiV2[fxtPhiV2CentKey]->Fill(invMass,refWgt/res12); // centrality dependence

    std::string fxtPhiV2MinBiasKey = Form("h_mInvMassFxtPhiV2%sCent0Pt%dPsi%d",str_mMixEvt[mFlagME].c_str(),ptBin,PsiBin);
    h_mInvMassFxtPhiV2[fxtPhiV2MinBiasKey]->Fill(invMass,refWgt/res12); // minimum bias
  }
}

void StPhiMesonHistoManger::fillFxtPhiYields(int cent9, double pt, double yCms, double invMass, double refWgt)
{
  int ptBin = getPtBinFxtPhiFlow(pt);

  if(ptBin >= 0)
  {
    std::string fxtPhiYieldsKey = Form("h_mInvMassFxtPhiYields%sCent%dPt%d",str_mMixEvt[mFlagME].c_str(),cent9,ptBin);
    h_mInvMassFxtPhiYields[fxtPhiYieldsKey]->Fill(invMass,refWgt);
  }
}

void StPhiMesonHistoManger::writeFxtPhiFlow()
{
  for(int iCent = 0; iCent < mNumCentBinFxtPhiFlow; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinFxtPhiFlow; ++iPt)
    {
      for(int iPsi = 0; iPsi < mNumPsiBinFxtPhiFlow; ++iPsi)
      {
	std::string fxtPhiV2Key = Form("h_mInvMassFxtPhiV2%sCent%dPt%dPsi%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt,iPsi);
	h_mInvMassFxtPhiV2[fxtPhiV2Key]->Write();
      }
    }
  }

  for(int iCent = 0; iCent < mNumCentBinFxtPhiYileds; ++iCent)
  {
    for(int iPt = 0; iPt < mNumPtBinFxtPhiFlow; ++iPt)
    {
      std::string fxtPhiYieldsKey = Form("h_mInvMassFxtPhiYields%sCent%dPt%d",str_mMixEvt[mFlagME].c_str(),iCent,iPt);
      h_mInvMassFxtPhiYields[fxtPhiYieldsKey]->Write();
    }
  }
}

int StPhiMesonHistoManger::getCentBinFxtPhiFlow(int cent9)
{
  int centBin = -1;

  if(cent9 >= 7 && cent9 <= 8) return 1; //  0-10%
  if(cent9 >= 4 && cent9 <= 6) return 2; // 10-40%
  if(cent9 >= 2 && cent9 <= 3) return 3; // 40-60%

  return centBin;
}

int StPhiMesonHistoManger::getPtBinFxtPhiFlow(double pt)
{
  int ptBin = -1;

  const double ptLo[mNumPtBinFxtPhiFlow] = {0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
  const double ptHi[mNumPtBinFxtPhiFlow] = {0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,3.0};

  if(std::abs(pt-ptLo[0]) < std::numeric_limits<double>::epsilon()) ptBin = 0;
  for(int iPt = 0; iPt < mNumPtBinFxtPhiFlow; ++iPt)
  {
    if((pt > ptLo[iPt]) && (pt <= ptHi[iPt]))
    {
      ptBin = iPt;
    }
  }

  return ptBin;
}

int StPhiMesonHistoManger::getPsi12BinFxtPhiFlow(double phi, double Psi1)
{
  int PsiBin = -1;

  const double phiPsi12 = transPsi12(phi,Psi1); // [-pi/2,pi/2]
  if(!isPsi12InRange(phiPsi12)) return PsiBin;
  
  const double phiPsiMin      = 0.0;
  const double phiPsiMax      = TMath::Pi()/2.0;
  const double phiPsiBinWidth = (phiPsiMax-phiPsiMin)/(double)mNumPsiBinFxtPhiFlow; // pi/10

  if(std::abs(std::abs(phiPsi12)-phiPsiMin) < std::numeric_limits<double>::epsilon()) PsiBin = 0;
  for(int iPsi = 0; iPsi < mNumPsiBinFxtPhiFlow; ++iPsi)
  {
    if((std::abs(phiPsi12) > phiPsiMin+iPsi*phiPsiBinWidth) && (std::abs(phiPsi12) <= phiPsiMin+(iPsi+1)*phiPsiBinWidth))
    {
      PsiBin = iPsi;
    }
  }

  return PsiBin;
}

double StPhiMesonHistoManger::transPsi12(double phi, double Psi1)
{
  const double phiPsi12 = phi - Psi1; // [-2pi,2pi]
  double phiPsi12Corr = phiPsi12;
  if(phiPsi12 >  3.0*TMath::Pi()/2.0) phiPsi12Corr = phiPsi12 - 2.0*TMath::Pi();
  if(phiPsi12 < -3.0*TMath::Pi()/2.0) phiPsi12Corr = phiPsi12 + 2.0*TMath::Pi();
  if(phiPsi12 >  TMath::Pi()/2.0 && phiPsi12 <=  3.0*TMath::Pi()/2.0) phiPsi12Corr = phiPsi12 - TMath::Pi();
  if(phiPsi12 < -TMath::Pi()/2.0 && phiPsi12 >= -3.0*TMath::Pi()/2.0) phiPsi12Corr = phiPsi12 + TMath::Pi();

  return phiPsi12Corr;
}

bool StPhiMesonHistoManger::isPsi12InRange(double phiPsi12)
{
  if(phiPsi12 < -TMath::Pi()/2.0 || phiPsi12 > TMath::Pi()/2.0)
  {
    return false;
  }

  return true;
}
