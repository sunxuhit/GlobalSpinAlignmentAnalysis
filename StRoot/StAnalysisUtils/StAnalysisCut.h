#ifndef StAnalysisCut_h
#define StAnalysisCut_h

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoEpdHit;

class StAnalysisCut : public TObject
{
  public:
    StAnalysisCut(int beamType);
    virtual ~StAnalysisCut();

    // Run Cuts
    bool isFixedTarget();
    bool isIsobar();

    // Event Cuts
    bool isMinBias(StPicoEvent *picoEvent);
    bool isPileUpEvent(double refMult, double numOfBTofMatch, double vz);
    bool passEventCut(StPicoEvent *picoEvent);

    // Track Cuts
    bool passTrackBasic(StPicoTrack *picoTrack);
    bool passTrackQA(StPicoTrack *picoTrack, TVector3 primVtx);
    // TPC EP
    bool passTrackTpcEpFull(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passTrackTpcEpEast(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passTrackTpcEpWest(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passNumTrackTpcSubEpRaw(int numTrackEast, int numTrackWest);
    bool passNumTrackTpcSubEpReCenter(int numTrackEast, int numTrackWest);
    // EPD EP
    bool passHitEpdEpFull(StPicoEpdHit *picoEpdHit);
    bool passHitEpdEpEast(StPicoEpdHit *picoEpdHit);
    bool passHitEpdEpWest(StPicoEpdHit *picoEpdHit);

  private:
    const int mType;

    ClassDef(StAnalysisCut,1)
};
#endif
