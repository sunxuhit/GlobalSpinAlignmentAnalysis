#ifndef StPhiMesonEvent_h
#define StPhiMesonEvent_h

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include "TVector3.h"
// #include "TLorentzVector.h"

// A. Schmah 19.12.2011

class StPhiMesonTrack : public TObject
{
  private:
    // Track properties
    // TrackA: for phi => Kplus
    // TrackB: for phi => Kminus
    double mMass2A; // mass2 of track A
    double mMass2B;
    double mNSigKaonA; // nsigma dE/dx of particle A
    double mNSigKaonB;
    double mDcaA; // distance of closest approach of particle A
    double mDcaB;
    int mChargeA; // charge of particle A
    int mChargeB;
    TVector3 v_mTrkMomA; // Momentum for track A (px, py, pz)
    TVector3 v_mTrkMomB;
    int mFlagA; // Flag for Event A: 0 for Same Event, others for Mixed Event
    int mFlagB;

  public:
    StPhiMesonTrack()
    {
      mMass2A    = -999.9;
      mMass2B    = -999.9;
      mNSigKaonA = -999.9;
      mNSigKaonB = -999.9;
      mDcaA      = -999.9;
      mDcaB      = -999.9;
      mChargeA   = -999;
      mChargeB   = -999;
      v_mTrkMomA.SetXYZ(0.0,0.0,0.0);
      v_mTrkMomB.SetXYZ(0.0,0.0,0.0);
      mFlagA = -1;
      mFlagB = -1;
    }
    ~StPhiMesonTrack() {}

    // setters
    void setMass2A(double f)           { mMass2A    = f; }
    void setMass2B(double f)           { mMass2B    = f; }
    void setNSigKaonA(double f)        { mNSigKaonA = f; }
    void setNSigKaonB(double f)        { mNSigKaonB = f; }
    void setDcaA(double f)             { mDcaA      = f; }
    void setDcaB(double f)             { mDcaB      = f; }
    void setTrkMomA(TVector3 f)        { v_mTrkMomA = f; }
    void setTrkMomB(TVector3 f)        { v_mTrkMomB = f; }
    void setFlagA(int f)               { mFlagA     = f; }
    void setFlagB(int f)               { mFlagB     = f; }

    // getters
    double getMass2A() const           { return mMass2A;    }
    double getMass2B() const           { return mMass2B;    }
    double getNSigKaonA() const        { return mNSigKaonA; }
    double getNSigKaonB() const        { return mNSigKaonB; }
    double getDcaA() const             { return mDcaA;      }
    double getDcaB() const             { return mDcaB;      }
    TVector3 getTrkMomA() const        { return v_mTrkMomA; }
    TVector3 getTrkMomB() const        { return v_mTrkMomB; }
    int getFlagA() const               { return mFlagA;     }
    int getFlagB() const               { return mFlagB;     }

    ClassDef(StPhiMesonTrack,1)  // A simple track of a particle
};

class StPhiMesonEvent : public TObject
{
  private:
    int    mRunId;
    int    mEventId;
    int    mRefMult;
    int    mCent9;
    int    mCent16;
    int    mNumTofMatch;
    double mRefWgt;
    double mZDCx;
    double mBBCx;
    double mVzVpd;
    TVector3 mPrimVtx;

    bool mFlagZdcEp // 0 for no ZDC EP | 1 for ZDC EP
    TVector2 v_mQ1ZdcShiftEast; // Q1Vector
    TVector2 v_mQ1ZdcShiftWest; 
    TVector2 v_mQ1ZdcShiftFull; // shift corrected v_mQ1ZdcShiftWest-v_mQ1ZdcShiftEast 

    bool mFlagEpdEp // 0 for no EPD EP | 1 for EPD EP
    TVector2 v_mQ1EpdShiftEast; // Q1Vector
    TVector2 v_mQ1EpdShiftWest; 
    TVector2 v_mQ1EpdShiftFull; // shift corrected v_mQ1EpdShiftWest-v_mQ1EpdShiftEast 

    bool mFlagTpcEp // 0 for no TPC EP | 1 for TPC EP
    TVector2 v_mQ2TpcReCtrEast; // Q2Vector
    TVector2 v_mQ2TpcReCtrWest;
    TVector2 v_mQ3TpcReCtrEast; // Q3Vector
    TVector2 v_mQ3TpcReCtrWest;
    int mNumTrkReCtrEast;
    int mNumTrkReCtrWest;

    unsigned short fNumTracks;

    TClonesArray* fTracks; //->

  public:
    StPhiMesonEvent()
    {
      mRunId       = -1;
      mEventId     = -1;
      mRefMult     = -1;
      mCent9       = -1;
      mCent16      = -1;
      mNumTofMatch = -1;
      mRefWgt      = -1.0;
      mZDCx        = -1.0;
      mBBCx        = -1.0;
      mVzVpd       = -1.0;
      mPrimVtx.Set(0.0,0.0,0.0);

      mFlagZdcEp = false; // ZDC EP
      v_mQ1ZdcShiftEast.Set(0.0,0.0);
      v_mQ1ZdcShiftWest.Set(0.0,0.0);
      v_mQ1ZdcShiftFull.Set(0.0,0.0);

      mFlagEpdEp = false; // EPD EP
      v_mQ1EpdShiftEast.Set(0.0,0.0);
      v_mQ1EpdShiftWest.Set(0.0,0.0);
      v_mQ1EpdShiftFull.Set(0.0,0.0);

      mFlagTpcEp = false; // TPC EP
      v_mQ2TpcReCtrEast.Set(0.0,0.0);
      v_mQ2TpcReCtrWest.Set(0.0,0.0);
      v_mQ3TpcReCtrEast.Set(0.0,0.0);
      v_mQ3TpcReCtrWest.Set(0.0,0.0);
      mNumTrkReCtrEast = -1;
      mNumTrkReCtrWest = -1;

      fNumTracks = 0;

      fTracks = new TClonesArray( "StPhiMesonTrack", 10 );
    }

    ~StPhiMesonEvent()
    {
      delete fTracks;
      fTracks = NULL;
    }

    void setRunId(int r)               { mRunId            = r; }
    void setEventId(int r)             { mEventId          = r; }
    void setRefMult(int r)             { mRefMult          = r; }
    void setCentrality9(int r)         { mCent9            = r; }
    void setCentrality16(int r)        { mCent16           = r; }
    void setNumTofMatch(int r)         { mNumTofMatch      = r; }
    void setRefWgt(double r)           { mRefWgt           = r; }
    void setZDCx(double r)             { mZDCx             = r; }
    void setBBCx(double r)             { mBBCx             = r; }
    void setVzVpd(double r)            { mVzVpd            = r; }
    void setPrimVtx(TVector3 r)        { mPrimVtx          = r; }
    void setFlagZdcEp(bool r)          { mFlagZdcEp        = r; }
    void setQ1VecZdcEast(TVector2 r)   { v_mQ1ZdcShiftEast = r; }
    void setQ1VecZdcWest(TVector2 r)   { v_mQ1ZdcShiftWest = r; }
    void setQ1VecZdcFull(TVector2 r)   { v_mQ1ZdcShiftFull = r; }
    void setFlagEpdEp(bool r)          { mFlagEpdEp        = r; }
    void setQ1VecEpdEast(TVector2 r)   { v_mQ1EpdShiftEast = r; }
    void setQ1VecEpdWest(TVector2 r)   { v_mQ1EpdShiftWest = r; }
    void setQ1VecEpdFull(TVector2 r)   { v_mQ1EpdShiftFull = r; }
    void setFlagTpcEp(bool r)          { mFlagTpcEp        = r; }
    void setQ2VecTpcEast(TVector2 r)   { v_mQ2TpcReCtrEast = r; }
    void setQ2VecTpcWest(TVector2 r)   { v_mQ2TpcReCtrWest = r; }
    void setQ3VecTpcEast(TVector2 r)   { v_mQ3TpcReCtrEast = r; }
    void setQ3VecTpcWest(TVector2 r)   { v_mQ3TpcReCtrWest = r; }
    void setNumTrkReCtrEast(int r)     { mNumTrkReCtrEast  = r; }
    void setNumTrkReCtrWest(int r)     { mNumTrkReCtrWest  = r; }

    int getRunId() const               { return mRunId;            }
    int getEventId() const             { return mEventId;          }
    int getRefMult() const             { return mRefMult;          }
    int getCentrality9() const         { return mCent9;            }
    int getCentrality16() const        { return mCent16;           }
    int getNumTofMatch() const         { return mNumTofMatch;      }
    double getRefWgt() const           { return mRefWgt;           }
    double getZDCx() const             { return mZDCx;             }
    double getBBCx() const             { return mBBCx;             }
    double getVzVpd() const            { return mVzVpd;            }
    TVector3 getPrimVtx() const        { return mPrimVtx;          }
    bool getFlagZdcEp() const          { return mFlagZdcEp;        }
    TVector2 getQ1VecZdcEast() const   { return v_mQ1ZdcShiftEast; }
    TVector2 getQ1VecZdcWest() const   { return v_mQ1ZdcShiftWest; }
    TVector2 getQ1VecZdcFull() const   { return v_mQ1ZdcShiftFull; }
    bool getFlagEpdEp() const          { return mFlagEpdEp;        }
    TVector2 getQ1VecEpdEast() const   { return v_mQ1EpdShiftEast; }
    TVector2 getQ1VecEpdWest() const   { return v_mQ1EpdShiftWest; }
    TVector2 getQ1VecEpdFull() const   { return v_mQ1EpdShiftFull; }
    bool getFlagTpcEp() const          { return mFlagTpcEp;        }
    TVector2 getQ2VecTpcEast() const   { return v_mQ2TpcReCtrEast; }
    TVector2 getQ2VecTpcWest() const   { return v_mQ2TpcReCtrWest; }
    TVector2 getQ3VecTpcEast() const   { return v_mQ3TpcReCtrEast; }
    TVector2 getQ3VecTpcWest() const   { return v_mQ3TpcReCtrWest; }
    int getNumTrkReCtrEast() const     { return mNumTrkReCtrEast;  }
    int getNumTrkReCtrWest() const     { return mNumTrkReCtrWest;  }
    // -----------------------------------Number of Tracks----------------------------------------
    StPhiMesonTrack* createTrack()
    {
      if (fNumTracks == fTracks->GetSize()) { fTracks->Expand( fNumTracks + 10 ); }
      if (fNumTracks >= 10000)
      {
	Fatal( "StPhiMesonEvent::createTrack()", "ERROR: Too many tracks (>10000)!" );
	exit( 2 );
      }

      new((*fTracks)[fNumTracks++]) StPhiMesonTrack;
      return (StPhiMesonTrack*)((*fTracks)[fNumTracks - 1]);
    }

    void clearTrackList()
    {
      fNumTracks = 0;
      fTracks->Clear();
    }

    UShort_t getNumTracks() const
    {
      return fNumTracks;
    }

    StPhiMesonTrack* getTrack(UShort_t i) const
    {
      return i < fNumTracks ? (StPhiMesonTrack*)((*fTracks)[i]) : NULL;
    }

    ClassDef(StPhiMesonEvent,1)  // A simple event compiled of tracks
};

#endif // StPhiMesonEvent_h
