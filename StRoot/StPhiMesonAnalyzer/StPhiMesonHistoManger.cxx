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
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    std::string histName = Form("h_mInvMassPhi%sCent%d",str_mMixEvt[mFlagME].c_str(),iCent);
    h_mInvMassPhi[iCent] = new TH2F(histName.c_str(),histName.c_str(),mNumPtBinQA,0.0,5.0,anaUtils::mNumInvMassPhi,anaUtils::mMassPhiMin,anaUtils::mMassPhiMax);
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
  int ptBin = getPtBinQA(pt);
  int yLabBin = getRapBinQA(yLab);
  int yCmsBin = getRapBinQA(yCms);
  h_mInvMassPhi[cent9]->Fill(pt,invMass,refWgt);
  if(is2060(cent9)) h_mInvMassPhi[9]->Fill(pt,invMass,refWgt);

  if(ptBin >= 0 && yLabBin >= 0 && yCmsBin >= 0)
  {
    std::string acptPhiLabKey = Form("h_mAcptPhi%sLabCent%dPt%dRap%d",str_mMixEvt[mFlagME].c_str(),cent9,ptBin,yLabBin);
    h_mAcptPhiLab[acptPhiLabKey]->Fill(invMass,refWgt);
    std::string acptPhiCmsKey = Form("h_mAcptPhi%sCmsCent%dPt%dRap%d",str_mMixEvt[mFlagME].c_str(),cent9,ptBin,yCmsBin);
    h_mAcptPhiCms[acptPhiCmsKey]->Fill(invMass,refWgt);
  }
}

void StPhiMesonHistoManger::writePhiQA()
{
  for(int iCent = 0; iCent < mNumCentrality; ++iCent)
  {
    h_mInvMassPhi[iCent]->Write();
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
//-------------------------------------------------------------
int StPhiMesonHistoManger::getPtBinQA(double pt)
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

int StPhiMesonHistoManger::getRapBinQA(double y)
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

bool StPhiMesonHistoManger::is2060(int cent9)
{
  if(cent9 < 2 || cent9 >5) return false;

  return true;
}
