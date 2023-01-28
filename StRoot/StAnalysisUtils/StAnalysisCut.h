#ifndef StAnalysisCut_h
#define StAnalysisCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;

class StAnalysisCut : public TObject
{
  public:
    StAnalysisCut(int beamType);
    virtual ~StAnalysisCut();

    // Run Cuts
    bool isFixedTarget();
    bool isIsobar();

    // Event Cuts
    bool isMinBias(StPicoEvent*);
    bool isPileUpEvent(double, double, double); // refmult & nTofMatch & vz
    bool passEventCut(StPicoDst*);

    // Track Cuts
    bool passTrackBasic(StPicoTrack*);
    bool passTrackQA(StPicoTrack*, StPicoEvent*);
    bool passTrackEP(StPicoTrack*, StPicoEvent*);

  private:
    const int mType;

    ClassDef(StAnalysisCut,1)
};
#endif
