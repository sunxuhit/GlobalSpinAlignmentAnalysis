#ifndef StEventPlaneProManager_h
#define StEventPlaneProManager_h

#include <TString.h>
#include <TVector2.h>
#include <TProfile2D.h>
#include "StMessMgr.h"

class TProfile;

class StEventPlaneProManager
{
  public:
    StEventPlaneProManager();
    virtual ~StEventPlaneProManager();

    // ZDC-SMD ReCenter Correction
    void initZdcReCenter();
    void fillZdcReCenterEast(TVector2 qVector, int Cent9, int RunIndex); // vz_sign = vertex pos/neg
    void fillZdcReCenterWest(TVector2 qVector, int Cent9, int RunIndex);
    void writeZdcReCenter();

  private:
    // ZDC-SMD ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mZdcQEastVertical; // 0 = vertex pos/neg
    TProfile2D *p_mZdcQEastHorizontal;
    TProfile2D *p_mZdcQWestVertical;
    TProfile2D *p_mZdcQWestHorizontal;

    ClassDef(StEventPlaneProManager,1)
};

#endif
