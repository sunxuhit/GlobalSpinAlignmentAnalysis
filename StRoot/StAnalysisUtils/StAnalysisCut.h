#ifndef StAnalysisCut_h
#define StAnalysisCut_h

#include "TObject.h"
#include "TString.h"
#include "TVector2.h"
#include "TVector3.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoEpdHit;
class StPhiMesonTrack;

class StAnalysisCut : public TObject
{
  public:
    StAnalysisCut(int beamType);
    virtual ~StAnalysisCut();

    // Run Cuts
    bool isIsobar();
    bool isFxt3p85GeV_2018();

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
    // Kaon Candidate: used in StPhiMesonMaker && loose cuts to fill the TTree
    bool passTrkTpcKaonFull(StPicoTrack *picoTrack, TVector3 primVtx); // apply to single K+/K- track
    bool passTrkTpcKaonEast(StPicoTrack *picoTrack, TVector3 primVtx); 
    bool passTrkTpcKaonWest(StPicoTrack *picoTrack, TVector3 primVtx); 
    bool passTrkTofKaonMass(TVector3 primMom, int charge, double mass2);
    bool passTrkTofKaonBeta(TVector3 primMom, int charge, double beta); // cross-check with Guannan Xie
    // Kaon Candidate: used in StPhiMesonAnalyzer && strict cuts for flow/alignment analysis
    bool passTrkTpcKaonFull(StPhiMesonTrack *phiTrk); // apply to K+&K- track pairs
    bool passTrkTpcKaonEast(StPhiMesonTrack *phiTrk, int charge); // apply to single K+/K- track && mostly to K+
    bool passTrkTpcKaonWest(StPhiMesonTrack *phiTrk, int charge);
    bool passTrkTofKaonMass(StPhiMesonTrack *phiTrk); // apply to K+&K- track pairs
    bool passTrkTofKaonBeta(StPhiMesonTrack *phiTrk); // cross-check with Guannan Xie

    // EPD Hit Cuts for EPD EP
    bool passHitEpdEpFull(StPicoEpdHit *picoEpdHit);
    bool passHitEpdEpEast(StPicoEpdHit *picoEpdHit);
    bool passHitEpdEpWest(StPicoEpdHit *picoEpdHit);
    bool passHitEpdFlowEast(StPicoEpdHit *picoEpdHit);
    bool passHitEpdFlowWest(StPicoEpdHit *picoEpdHit);
    bool passQVecEpdSide(TVector2 Q1VecEast, TVector2 Q1VecWest, TVector2 Q1VecFull);
    bool passQVecEpdGrp(TVector2 Q1VecEast, TVector2 Q1VecWest, TVector2 Q1VecFull, int grpId);

    // ZDC Hit Cuts for ZDC EP
    bool passQVecZdc(TVector2 Q1VecEast, TVector2 Q1VecWest, TVector2 Q1VecFull);

    // only used for proton & deuteron flow comparison in Fxt3p85GeV_2018
    bool passTrkTpcFlow(StPicoTrack *picoTrack, TVector3 primVtx);
    bool passTrkProFlow(StPicoTrack *picoTrack);
    bool passTrkDeuFlow(double pMag, double deuteronZ, double mass2);

    // TPC EP: used in StPhiMesonAnalyzer
    bool passTrkTpcEpFull(TVector3 primMom, double gDca);
    bool passTrkTpcEpEast(TVector3 primMom, double gDca);
    bool passTrkTpcEpWest(TVector3 primMom, double gDca);
    // phi Flow: used in StPhiMesonAnalyzer
    bool passTrkPhiFlowEast(double yPhiCms);
    bool passTrkPhiFlowWest(double yPhiCms);

  private:
    const int mType;

    ClassDef(StAnalysisCut,1)
};
#endif
