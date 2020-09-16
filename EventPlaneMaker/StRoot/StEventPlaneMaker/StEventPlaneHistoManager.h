#ifndef StEventPlaneHistoManager_h
#define StEventPlaneHistoManager_h

#include "StMessMgr.h"
#include "TVector3.h"

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
    //--------------ZDC EP---------------

    /*
    void InitEP();
    void FillEP_Sub(Float_t Psi2East_ReCenter, Float_t Psi2East_Shift, Float_t Psi2West_ReCenter, Float_t Psi2West_Shift);
    void FillEP_Ran(Float_t Psi2RanA_ReCenter, Float_t Psi2RanA_Shift, Float_t Psi2RanB_ReCenter, Float_t Psi2RanB_Shift, Float_t Psi2Full_ReCenter, Float_t Psi2Full_Shift);
    void WriteEP();
    */
    
  private:
    //--------------ZDC EP---------------
    TH2F *h_mZdcGainCorr[2][2][8]; // 0: east/west | 1: vertical(x)/horizontal(y) | 2: 7 slats(x)/8 slats(y); | x-axis: runIndex | y-axis: ADC
    //--------------ZDC EP---------------

  ClassDef(StEventPlaneHistoManager,1)
};
#endif
