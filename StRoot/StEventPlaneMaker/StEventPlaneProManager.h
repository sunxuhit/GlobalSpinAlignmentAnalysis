#ifndef StEventPlaneProManager_h
#define StEventPlaneProManager_h

#include <TVector2.h>
#include <TProfile2D.h>
// #include <TString.h>
// #include "StMessMgr.h"

class TProfile;

class StEventPlaneProManager
{
  public:
    StEventPlaneProManager(int beamType);
    virtual ~StEventPlaneProManager();

    // ZDC-SMD ReCenter Correction
    void initZdcReCenter();
    void fillZdcReCenterEast(TVector2 qVector, int Cent9, int RunIndex, int vzBin);
    void fillZdcReCenterWest(TVector2 qVector, int Cent9, int RunIndex, int vzBin);
    void writeZdcReCenter();

  private:
    static const int mNumVzBin = 2; // 0: vz < 0 | 1: vz >= 0

    // ZDC-SMD ReCenter Correction | x axis is RunIndex, y axis is Centrality
    TProfile2D *p_mZdcQEastVertical[mNumVzBin]; // 0: vz < 0 | 1: vz >= 0
    TProfile2D *p_mZdcQEastHorizontal[mNumVzBin];
    TProfile2D *p_mZdcQWestVertical[mNumVzBin];
    TProfile2D *p_mZdcQWestHorizontal[mNumVzBin];

    const int mType;

    ClassDef(StEventPlaneProManager,1)
};

#endif
