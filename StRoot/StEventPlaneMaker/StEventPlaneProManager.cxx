#include <string>

#include <TProfile.h>
#include <TString.h>
#include <TMath.h>

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneProManager.h"
// #include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

ClassImp(StEventPlaneProManager)

//---------------------------------------------------------------------------------

StEventPlaneProManager::StEventPlaneProManager(int beamType) : mType(beamType)
{
}

StEventPlaneProManager::~StEventPlaneProManager()
{
  /* */
}

//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// ZDC-SMD ReCenter Correction
void StEventPlaneProManager::initZdcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    std::string ProName = Form("p_mZdcQEastVertVz%d",iVz);
    p_mZdcQEastVertical[iVz] = new TProfile2D(ProName.c_str(),ProName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    ProName = Form("p_mZdcQEastHoriVz%d",iVz);
    p_mZdcQEastHorizontal[iVz] = new TProfile2D(ProName.c_str(),ProName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);

    ProName = Form("p_mZdcQWestVertVz%d",iVz);
    p_mZdcQWestVertical[iVz] = new TProfile2D(ProName.c_str(),ProName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
    ProName = Form("p_mZdcQWestHoriVz%d",iVz);
    p_mZdcQWestHorizontal[iVz] = new TProfile2D(ProName.c_str(),ProName.c_str(),globCons::mMaxRunIndex[mType],-0.5,(double)globCons::mMaxRunIndex[mType]-0.5,9,-0.5,8.5);
  }
}

void StEventPlaneProManager::fillZdcReCenterEast(TVector2 qVector, int Cent9, int RunIndex, int vzBin)
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  // Event Plane method
  p_mZdcQEastVertical[vzBin]->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mZdcQEastHorizontal[vzBin]->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StEventPlaneProManager::fillZdcReCenterWest(TVector2 qVector, int Cent9, int RunIndex, int vzBin)
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  p_mZdcQWestVertical[vzBin]->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mZdcQWestHorizontal[vzBin]->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StEventPlaneProManager::writeZdcReCenter()
{
  for(int iVz = 0; iVz < mNumVzBin; ++iVz)
  {
    p_mZdcQEastVertical[iVz]->Write();
    p_mZdcQEastHorizontal[iVz]->Write();
    p_mZdcQWestVertical[iVz]->Write();
    p_mZdcQWestHorizontal[iVz]->Write();
  }
}
//---------------------------------------------------------------------------------

