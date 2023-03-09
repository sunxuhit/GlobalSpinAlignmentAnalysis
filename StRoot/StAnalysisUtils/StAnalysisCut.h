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
    bool isGoodCent9(int cent9);
    bool passEventCut(StPicoEvent *picoEvent);

    // Track Cuts
    bool passTrkBasic(StPicoTrack *picoTrack);
    bool passTrkQA(StPicoTrack *picoTrack, TVector3 primVtx);
    // TPC EP
    bool passTrkTpcEpFull(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passTrkTpcEpEast(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passTrkTpcEpWest(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passNumTrkTpcSubEpRaw(int numTrackEast, int numTrackWest);
    bool passNumTrkTpcSubEpReCtr(int numTrackEast, int numTrackWest);
    // TPC Flow
    bool passTrkTpcFlowFull(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passTrkTpcFlowEast(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passTrkTpcFlowWest(StPicoTrack *picoTrack, TVector3 primVtx);
    // Kaon Candidate
    bool passTrkKaonFull(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passTrkKaonEast(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passTrkKaonWest(StPicoTrack *picoTrack, TVector3 primVtx);

    // EPD Hit Cuts for EPD EP
    bool passHitEpdEpFull(StPicoEpdHit *picoEpdHit);
    bool passHitEpdEpEast(StPicoEpdHit *picoEpdHit);
    bool passHitEpdEpWest(StPicoEpdHit *picoEpdHit);
    bool passHitEpdFlowEast(StPicoEpdHit *picoEpdHit);
    bool passHitEpdFlowWest(StPicoEpdHit *picoEpdHit);

  private:
    const int mType;

    ClassDef(StAnalysisCut,1)
};
#endif
