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

    bool isBES();
    bool isIsobar();
    bool isMinBias(StPicoEvent*);
    bool isPileUpEvent(double, double, double); // refmult & nTofMatch & vz
    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackQA(StPicoTrack*, StPicoEvent*);

  private:
    int mType;

    ClassDef(StAnalysisCut,1)
};
#endif
