#ifndef StVecMesonCut_h
#define StVecMesonCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;

class StVecMesonCut : public TObject
{
  public:
    StVecMesonCut(int energy);
    virtual ~StVecMesonCut();

    bool isMinBias(StPicoEvent*);
    bool isBES(int energy);
    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackEP(StPicoTrack*);
    bool passSigPionCut(StPicoTrack*, float);
    bool passSigKaonCut(StPicoTrack*, float);
    bool passSigProntonCut(StPicoTrack*, float);
    bool passTrackPhi(StPicoTrack*);
    int getMatchedToF();
    int getNpirm();
    int getNnonprim();
    float getPrimaryMass2(StPicoDst*, int); // return m^2 of i-th track (primary || -999)
    float getGlobalMass2(StPicoDst*, int); // return m^2 of i-th track (global || -999)

  private:
    int mMatchedToF;
    int mN_prim;
    int mN_non_prim;
    int mEnergy;

    ClassDef(StVecMesonCut,1)
};
#endif
