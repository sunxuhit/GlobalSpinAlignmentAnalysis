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
    TVector3 v_mTrkMomKp; // Momentum for Kplus (px, py, pz)
    TVector3 v_mTrkMomKm;
    double mMass2Kp; // mass2 of Kplus
    double mMass2Km;
    double mNSigKp; // nsigma dE/dx of Kplus
    double mNSigKm;
    double mDcaKp; // distance of closest approach of Kplus
    double mDcaKm;
    int mChargeKp; // charge of Kplus
    int mChargeKm;
    int mNHitsFitKp; // nHitsFit of Kplus
    int mNHitsFitKm;
    int mFlagKp; // Flag for Kplus: 0 for Same Event, others for Mixed Event
    int mFlagKm;

  public:
    StPhiMesonTrack()
    {
      v_mTrkMomKp.SetXYZ(0.0,0.0,0.0);
      v_mTrkMomKm.SetXYZ(0.0,0.0,0.0);
      mMass2Kp    = -999.9;
      mMass2Km    = -999.9;
      mNSigKp     = -999.9;
      mNSigKm     = -999.9;
      mDcaKp      = -999.9;
      mDcaKm      = -999.9;
      mChargeKp   = -999;
      mChargeKm   = -999;
      mNHitsFitKp = -999;
      mNHitsFitKm = -999;
      mFlagKp     = -1;
      mFlagKm     = -1;
    }
    ~StPhiMesonTrack() {}

    // setters
    void setTrkMomKp(TVector3 f)        { v_mTrkMomKp = f; }
    void setTrkMomKm(TVector3 f)        { v_mTrkMomKm = f; }
    void setMass2Kp(double f)           { mMass2Kp    = f; }
    void setMass2Km(double f)           { mMass2Km    = f; }
    void setNSigKp(double f)            { mNSigKp     = f; }
    void setNSigKm(double f)            { mNSigKm     = f; }
    void setDcaKp(double f)             { mDcaKp      = f; }
    void setDcaKm(double f)             { mDcaKm      = f; }
    void setChargeKp(int f)             { mChargeKp   = f; }
    void setChargeKm(int f)             { mChargeKm   = f; }
    void setNHitsFitKp(int f)           { mNHitsFitKp = f; }
    void setNHitsFitKm(int f)           { mNHitsFitKm = f; }
    void setFlagKp(int f)               { mFlagKp     = f; }
    void setFlagKm(int f)               { mFlagKm     = f; }

    // getters
    TVector3 getTrkMomKp() const        { return v_mTrkMomKp; }
    TVector3 getTrkMomKm() const        { return v_mTrkMomKm; }
    double getMass2Kp() const           { return mMass2Kp;    }
    double getMass2Km() const           { return mMass2Km;    }
    double getNSigKp() const            { return mNSigKp;     }
    double getNSigKm() const            { return mNSigKm;     }
    double getDcaKp() const             { return mDcaKp;      }
    double getDcaKm() const             { return mDcaKm;      }
    int getChargeKp()                   { return mChargeKp;   }
    int getChargeKm()                   { return mChargeKm;   }
    int getNHitsFitKp()                 { return mNHitsFitKp; }
    int getNHitsFitKm()                 { return mNHitsFitKm; }
    int getFlagKp() const               { return mFlagKp;     }
    int getFlagKm() const               { return mFlagKm;     }

    ClassDef(StPhiMesonTrack,1)  // A simple track of a particle
};

class StPhiMesonEvent : public TObject
{
  private:
    int    mRunId;
    int    mEvtId;
    int    mRefMult;
    int    mNumTofMatch;
    int    mCent9;
    int    mCent16;
    double mRefWgt;
    double mZDCx;
    double mBBCx;
    double mVzVpd;
    TVector3 mPrimVtx;

    int mFlagZdcEp; // 0 for no ZDC EP | 1 for ZDC EP
    TVector2 v_mQ1ZdcShiftEast; // Q1Vector
    TVector2 v_mQ1ZdcShiftWest; 
    TVector2 v_mQ1ZdcShiftFull; // shift corrected v_mQ1ZdcShiftWest-v_mQ1ZdcShiftEast 

    int mFlagEpdEp; // 0 for no EPD EP | 1 for EPD EP
    TVector2 v_mQ1EpdShiftEast; // Q1Vector
    TVector2 v_mQ1EpdShiftWest; 
    TVector2 v_mQ1EpdShiftFull; // shift corrected v_mQ1EpdShiftWest-v_mQ1EpdShiftEast 

    int mFlagTpcEp; // 0 for no TPC EP | 1 for TPC EP
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
      mEvtId       = -1;
      mRefMult     = -1;
      mNumTofMatch = -1;
      mCent9       = -1;
      mCent16      = -1;
      mRefWgt      = -1.0;
      mZDCx        = -1.0;
      mBBCx        = -1.0;
      mVzVpd       = -1.0;
      mPrimVtx.Set(-999.9,-999.9,-999.9);

      mFlagZdcEp = 0; // ZDC EP
      v_mQ1ZdcShiftEast.Set(0.0,0.0);
      v_mQ1ZdcShiftWest.Set(0.0,0.0);
      v_mQ1ZdcShiftFull.Set(0.0,0.0);

      mFlagEpdEp = 0; // EPD EP
      v_mQ1EpdShiftEast.Set(0.0,0.0);
      v_mQ1EpdShiftWest.Set(0.0,0.0);
      v_mQ1EpdShiftFull.Set(0.0,0.0);

      mFlagTpcEp = 0; // TPC EP
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
    void setEvtId(int r)               { mEvtId            = r; }
    void setRefMult(int r)             { mRefMult          = r; }
    void setNumTofMatch(int r)         { mNumTofMatch      = r; }
    void setCentrality9(int r)         { mCent9            = r; }
    void setCentrality16(int r)        { mCent16           = r; }
    void setRefWgt(double r)           { mRefWgt           = r; }
    void setZDCx(double r)             { mZDCx             = r; }
    void setBBCx(double r)             { mBBCx             = r; }
    void setVzVpd(double r)            { mVzVpd            = r; }
    void setPrimVtx(TVector3 r)        { mPrimVtx          = r; }
    void setFlagZdcEp(int r)           { mFlagZdcEp        = r; }
    void setQ1VecZdcEast(TVector2 r)   { v_mQ1ZdcShiftEast = r; }
    void setQ1VecZdcWest(TVector2 r)   { v_mQ1ZdcShiftWest = r; }
    void setQ1VecZdcFull(TVector2 r)   { v_mQ1ZdcShiftFull = r; }
    void setFlagEpdEp(int r)           { mFlagEpdEp        = r; }
    void setQ1VecEpdEast(TVector2 r)   { v_mQ1EpdShiftEast = r; }
    void setQ1VecEpdWest(TVector2 r)   { v_mQ1EpdShiftWest = r; }
    void setQ1VecEpdFull(TVector2 r)   { v_mQ1EpdShiftFull = r; }
    void setFlagTpcEp(int r)           { mFlagTpcEp        = r; }
    void setQ2VecTpcEast(TVector2 r)   { v_mQ2TpcReCtrEast = r; }
    void setQ2VecTpcWest(TVector2 r)   { v_mQ2TpcReCtrWest = r; }
    void setQ3VecTpcEast(TVector2 r)   { v_mQ3TpcReCtrEast = r; }
    void setQ3VecTpcWest(TVector2 r)   { v_mQ3TpcReCtrWest = r; }
    void setNumTrkReCtrEast(int r)     { mNumTrkReCtrEast  = r; }
    void setNumTrkReCtrWest(int r)     { mNumTrkReCtrWest  = r; }

    int getRunId() const               { return mRunId;            }
    int getEvtId() const               { return mEvtId;            }
    int getRefMult() const             { return mRefMult;          }
    int getNumTofMatch() const         { return mNumTofMatch;      }
    int getCentrality9() const         { return mCent9;            }
    int getCentrality16() const        { return mCent16;           }
    double getRefWgt() const           { return mRefWgt;           }
    double getZDCx() const             { return mZDCx;             }
    double getBBCx() const             { return mBBCx;             }
    double getVzVpd() const            { return mVzVpd;            }
    TVector3 getPrimVtx() const        { return mPrimVtx;          }
    int getFlagZdcEp() const           { return mFlagZdcEp;        }
    TVector2 getQ1VecZdcEast() const   { return v_mQ1ZdcShiftEast; }
    TVector2 getQ1VecZdcWest() const   { return v_mQ1ZdcShiftWest; }
    TVector2 getQ1VecZdcFull() const   { return v_mQ1ZdcShiftFull; }
    int getFlagEpdEp() const           { return mFlagEpdEp;        }
    TVector2 getQ1VecEpdEast() const   { return v_mQ1EpdShiftEast; }
    TVector2 getQ1VecEpdWest() const   { return v_mQ1EpdShiftWest; }
    TVector2 getQ1VecEpdFull() const   { return v_mQ1EpdShiftFull; }
    int getFlagTpcEp() const           { return mFlagTpcEp;        }
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
