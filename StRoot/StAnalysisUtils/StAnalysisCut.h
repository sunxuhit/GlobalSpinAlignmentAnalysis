#ifndef StAnalysisCut_h
#define StAnalysisCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class TVector3;

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
    bool passEventCut(StPicoEvent*);

    // Track Cuts
    bool passTrackBasic(StPicoTrack*);
    bool passTrackQA(StPicoTrack*, TVector3 primVtx);
    // TPC EP
    bool passTrackTpcEpFull(StPicoTrack*, TVector3 primVtx);
    bool passTrackTpcEpEast(StPicoTrack*, TVector3 primVtx);
    bool passTrackTpcEpWest(StPicoTrack*, TVector3 primVtx);
    bool passNumTrackTpcSubEpRaw(int numTrackEast, int numTrackWest);
    bool passNumTrackTpcSubEpReCenter(int numTrackEast, int numTrackWest);

  private:
    const int mType;

    ClassDef(StAnalysisCut,1)
};
#endif
