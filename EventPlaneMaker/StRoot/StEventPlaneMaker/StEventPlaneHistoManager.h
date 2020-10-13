#ifndef StEventPlaneHistoManager_h
#define StEventPlaneHistoManager_h

#include "StMessMgr.h"
#include "TVector2.h"

class TH1F;
class TH2F;

class StEventPlaneHistoManager
{
  public:
    StEventPlaneHistoManager();
    virtual ~StEventPlaneHistoManager();

    //--------------ZDC EP---------------
    void initZdcGainCorr();
    void fillZdcGainCorr(int i_eastwest, int i_verthori, int i_slat, int runIndex, float zdcsmd);
    void writeZdcGainCorr();

    void initZdcRawEP(); // raw EP
    void fillZdcRawEP(TVector2 QEast, TVector2 QWest, TVector2 QFull, int Cent9, int runIndex);
    void writeZdcRawEP();
    //--------------ZDC EP---------------
    
  private:
    //--------------ZDC EP---------------
    TH2F *h_mZdcGainCorr[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y); | x-axis: runIndex | y-axis: ADC

    TH2F *h_mZdcRawEast[9]; // raw EP
    TH2F *h_mZdcRawWest[9];
    TH2F *h_mZdcRawFull[9]; // Qwest-QEast
    //--------------ZDC EP---------------

  ClassDef(StEventPlaneHistoManager,1)
};
#endif
