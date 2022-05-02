#include <string>

#include <TProfile.h>
#include <TMath.h>
#include <TString.h>

#include "StRoot/StEventPlaneMaker/StEventPlaneProManager.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"

ClassImp(StEventPlaneProManager)

//---------------------------------------------------------------------------------

StEventPlaneProManager::StEventPlaneProManager()
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
  string ProName;

  ProName = "p_mZdcQEastVertical";
  p_mZdcQEastVertical = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  ProName = "p_mZdcQEastHorizontal";
  p_mZdcQEastHorizontal = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);

  ProName = "p_mZdcQWestVertical";
  p_mZdcQWestVertical = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
  ProName = "p_mZdcQWestHorizontal";
  p_mZdcQWestHorizontal = new TProfile2D(ProName.c_str(),ProName.c_str(),recoEP::mNumOfRunIndex,-0.5,(float)recoEP::mNumOfRunIndex-0.5,9,-0.5,8.5);
}

void StEventPlaneProManager::fillZdcReCenterEast(TVector2 qVector, int Cent9, int RunIndex) // vz_sign = vertex pos/neg 
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  // Event Plane method
  p_mZdcQEastVertical->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mZdcQEastHorizontal->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StEventPlaneProManager::fillZdcReCenterWest(TVector2 qVector, int Cent9, int RunIndex) // vz_sign = vertex pos/neg 
{
  const double qx = qVector.X();
  const double qy = qVector.Y();

  p_mZdcQWestVertical->Fill((double)RunIndex,(double)Cent9,(double)qx);
  p_mZdcQWestHorizontal->Fill((double)RunIndex,(double)Cent9,(double)qy);
}

void StEventPlaneProManager::writeZdcReCenter()
{
  p_mZdcQEastVertical->Write();
  p_mZdcQEastHorizontal->Write();
  p_mZdcQWestVertical->Write();
  p_mZdcQWestHorizontal->Write();
}
//---------------------------------------------------------------------------------

