#ifndef StEventPlaneCut_h
#define StEventPlaneCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;

class StEventPlaneCut : public TObject
{
  public:
    StEventPlaneCut(int energy);
    virtual ~StEventPlaneCut();

    bool isMinBias(StPicoEvent*);
    bool isBES();
    bool isPileUpEvent(int, int, int); // refmult/grefmult & nTofMatch & nTofHits
    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackQA(StPicoTrack*, StPicoEvent*);
    float getBeta(StPicoDst*, int); // return beta of i-th track (tof || -999)
    float getPrimaryMass2(StPicoDst*, int); // return m^2 of i-th track (primary || -999)
    float getGlobalMass2(StPicoDst*, int); // return m^2 of i-th track (global || -999)
    int getTriggerBin(StPicoEvent*); // return trigger bin for event QA

  private:
    // int mMatchedToF;
    // int mN_prim;
    // int mN_non_prim;
    int mEnergy;

    ClassDef(StEventPlaneCut,1)
};
#endif
