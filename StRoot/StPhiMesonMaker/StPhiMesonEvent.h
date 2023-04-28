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
    double mBetaKp; // 1/beta of Kplus
    double mBetaKm;
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
      mBetaKp     = -999.9;
      mBetaKm     = -999.9;
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
    void setBetaKp(double f)            { mBetaKp     = f; }
    void setBetaKm(double f)            { mBetaKm     = f; }
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
    double getBetaKp() const            { return mBetaKp;     }
    double getBetaKm() const            { return mBetaKm;     }
    double getNSigKp() const            { return mNSigKp;     }
    double getNSigKm() const            { return mNSigKm;     }
    double getDcaKp() const             { return mDcaKp;      }
    double getDcaKm() const             { return mDcaKm;      }
    int getChargeKp() const             { return mChargeKp;   }
    int getChargeKm() const             { return mChargeKm;   }
    int getNHitsFitKp() const           { return mNHitsFitKp; }
    int getNHitsFitKm() const           { return mNHitsFitKm; }
    int getFlagKp() const               { return mFlagKp;     }
    int getFlagKm() const               { return mFlagKm;     }

    ClassDef(StPhiMesonTrack,1)  // A simple track of a particle
};

class StPhiMesonEvent : public TObject
{
  private:
    int mRunId;
    int mRunIdx;
    int mEvtId;
    int mRefMult;
    int mNumTofMatch;
    int mCent9;
    int mCent16;
    double mRefWgt;
    double mZDCx;
    double mBBCx;
    double mVzVpd;
    TVector3 mPrimVtx;

    int mFlagZdcEp; // -1 for NOT Used | 0 for no ZDC EP | 1 for ZDC EP
    TVector2 v_mQ1ZdcShiftEast; // Q1Vector
    TVector2 v_mQ1ZdcShiftWest; 
    TVector2 v_mQ1ZdcShiftFull; // shift corrected v_mQ1ZdcShiftWest-v_mQ1ZdcShiftEast 

    int mFlagEpdSideEp; // -1 for NOT Used | 0 for no EPD Side EP | 1 for EPD Side EP
    TVector2 v_mQ1EpdSideShiftEast; // Q1Vector
    TVector2 v_mQ1EpdSideShiftWest; 
    TVector2 v_mQ1EpdSideShiftFull; // shift corrected v_mQ1EpdShiftWest-v_mQ1EpdShiftEast 

    int mFlagEpdGrp0Ep; // -1 for NOT Used | 0 for no EPD Grp0 EP | 1 for Grp0 Side EP
    TVector2 v_mQ1EpdGrp0ShiftEast; // Q1Vector
    TVector2 v_mQ1EpdGrp0ShiftWest; 
    TVector2 v_mQ1EpdGrp0ShiftFull; // shift corrected v_mQ1EpdShiftWest-v_mQ1EpdShiftEast 
    int mFlagEpdGrp1Ep; // -1 for NOT Used | 0 for no EPD Grp1 EP | 1 for EPD Grp1 EP
    TVector2 v_mQ1EpdGrp1ShiftEast; // Q1Vector
    TVector2 v_mQ1EpdGrp1ShiftWest; 
    TVector2 v_mQ1EpdGrp1ShiftFull; // shift corrected v_mQ1EpdShiftWest-v_mQ1EpdShiftEast 

    int mFlagTpcEp; // -1 for NOT Used | 0 for no TPC EP | 1 for TPC EP
    TVector2 v_mQ1TpcReCtrEast; // Q2Vector
    TVector2 v_mQ1TpcReCtrWest;
    TVector2 v_mQ2TpcReCtrEast; // Q2Vector
    TVector2 v_mQ2TpcReCtrWest;
    TVector2 v_mQ3TpcReCtrEast; // Q3Vector
    TVector2 v_mQ3TpcReCtrWest;
    int mNumTrkReCtrEast;
    int mNumTrkReCtrWest;

    unsigned short mNumTracks;

    TClonesArray* phiTracks; //->

  public:
    StPhiMesonEvent()
    {
      mRunId       = -1;
      mRunIdx      = -1;
      mEvtId       = -1;
      mRefMult     = -1;
      mNumTofMatch = -1;
      mCent9       = -1;
      mCent16      = -1;
      mRefWgt      = -1.0;
      mZDCx        = -1.0;
      mBBCx        = -1.0;
      mVzVpd       = -1.0;
      mPrimVtx.SetXYZ(-999.9,-999.9,-999.9);

      mFlagZdcEp = -1; // ZDC EP
      v_mQ1ZdcShiftEast.Set(0.0,0.0);
      v_mQ1ZdcShiftWest.Set(0.0,0.0);
      v_mQ1ZdcShiftFull.Set(0.0,0.0);

      mFlagEpdSideEp = -1; // EPD Side EP
      v_mQ1EpdSideShiftEast.Set(0.0,0.0);
      v_mQ1EpdSideShiftWest.Set(0.0,0.0);
      v_mQ1EpdSideShiftFull.Set(0.0,0.0);

      mFlagEpdGrp0Ep = -1; // EPD Grp0 EP
      v_mQ1EpdGrp0ShiftEast.Set(0.0,0.0);
      v_mQ1EpdGrp0ShiftWest.Set(0.0,0.0);
      v_mQ1EpdGrp0ShiftFull.Set(0.0,0.0);
      mFlagEpdGrp1Ep = -1; // EPD Grp1 EP
      v_mQ1EpdGrp1ShiftEast.Set(0.0,0.0);
      v_mQ1EpdGrp1ShiftWest.Set(0.0,0.0);
      v_mQ1EpdGrp1ShiftFull.Set(0.0,0.0);

      mFlagTpcEp = -1; // TPC EP
      v_mQ1TpcReCtrEast.Set(0.0,0.0);
      v_mQ1TpcReCtrWest.Set(0.0,0.0);
      v_mQ2TpcReCtrEast.Set(0.0,0.0);
      v_mQ2TpcReCtrWest.Set(0.0,0.0);
      v_mQ3TpcReCtrEast.Set(0.0,0.0);
      v_mQ3TpcReCtrWest.Set(0.0,0.0);
      mNumTrkReCtrEast = -1;
      mNumTrkReCtrWest = -1;

      mNumTracks = 0;

      phiTracks = new TClonesArray( "StPhiMesonTrack", 10 );
    }

    ~StPhiMesonEvent()
    {
      delete phiTracks;
      phiTracks = NULL;
    }

    void setRunId(int r)                   { mRunId                = r; } // Evt Header
    void setRunIdx(int r)                  { mRunIdx               = r; }
    void setEvtId(int r)                   { mEvtId                = r; }
    void setRefMult(int r)                 { mRefMult              = r; }
    void setNumTofMatch(int r)             { mNumTofMatch          = r; }
    void setCentrality9(int r)             { mCent9                = r; }
    void setCentrality16(int r)            { mCent16               = r; }
    void setRefWgt(double r)               { mRefWgt               = r; }
    void setZDCx(double r)                 { mZDCx                 = r; }
    void setBBCx(double r)                 { mBBCx                 = r; }
    void setVzVpd(double r)                { mVzVpd                = r; }
    void setPrimVtx(TVector3 r)            { mPrimVtx              = r; }

    void setFlagZdcEp(int r)               { mFlagZdcEp            = r; } // ZDC EP
    void setQ1VecZdcEast(TVector2 r)       { v_mQ1ZdcShiftEast     = r; }
    void setQ1VecZdcWest(TVector2 r)       { v_mQ1ZdcShiftWest     = r; }
    void setQ1VecZdcFull(TVector2 r)       { v_mQ1ZdcShiftFull     = r; }

    void setFlagEpdSideEp(int r)           { mFlagEpdSideEp        = r; } // EPD Side EP
    void setQ1VecEpdSideEast(TVector2 r)   { v_mQ1EpdSideShiftEast = r; }
    void setQ1VecEpdSideWest(TVector2 r)   { v_mQ1EpdSideShiftWest = r; }
    void setQ1VecEpdSideFull(TVector2 r)   { v_mQ1EpdSideShiftFull = r; }

    void setFlagEpdGrp0Ep(int r)           { mFlagEpdGrp0Ep        = r; } // EPD Grp0 EP
    void setQ1VecEpdGrp0East(TVector2 r)   { v_mQ1EpdGrp0ShiftEast = r; }
    void setQ1VecEpdGrp0West(TVector2 r)   { v_mQ1EpdGrp0ShiftWest = r; }
    void setQ1VecEpdGrp0Full(TVector2 r)   { v_mQ1EpdGrp0ShiftFull = r; }
    void setFlagEpdGrp1Ep(int r)           { mFlagEpdGrp1Ep        = r; } // EPD Grp1 EP
    void setQ1VecEpdGrp1East(TVector2 r)   { v_mQ1EpdGrp1ShiftEast = r; }
    void setQ1VecEpdGrp1West(TVector2 r)   { v_mQ1EpdGrp1ShiftWest = r; }
    void setQ1VecEpdGrp1Full(TVector2 r)   { v_mQ1EpdGrp1ShiftFull = r; }

    void setFlagTpcEp(int r)               { mFlagTpcEp            = r; } // TPC EP
    void setQ1VecTpcEast(TVector2 r)       { v_mQ1TpcReCtrEast     = r; }
    void setQ1VecTpcWest(TVector2 r)       { v_mQ1TpcReCtrWest     = r; }
    void setQ2VecTpcEast(TVector2 r)       { v_mQ2TpcReCtrEast     = r; }
    void setQ2VecTpcWest(TVector2 r)       { v_mQ2TpcReCtrWest     = r; }
    void setQ3VecTpcEast(TVector2 r)       { v_mQ3TpcReCtrEast     = r; }
    void setQ3VecTpcWest(TVector2 r)       { v_mQ3TpcReCtrWest     = r; }
    void setNumTrkReCtrEast(int r)         { mNumTrkReCtrEast      = r; }
    void setNumTrkReCtrWest(int r)         { mNumTrkReCtrWest      = r; }

    int getRunId() const                   { return mRunId;                } // Evt Header
    int getRunIdx() const                  { return mRunIdx;               }
    int getEvtId() const                   { return mEvtId;                }
    int getRefMult() const                 { return mRefMult;              }
    int getNumTofMatch() const             { return mNumTofMatch;          }
    int getCentrality9() const             { return mCent9;                }
    int getCentrality16() const            { return mCent16;               }
    double getRefWgt() const               { return mRefWgt;               }
    double getZDCx() const                 { return mZDCx;                 }
    double getBBCx() const                 { return mBBCx;                 }
    double getVzVpd() const                { return mVzVpd;                }
    TVector3 getPrimVtx() const            { return mPrimVtx;              }

    int getFlagZdcEp() const               { return mFlagZdcEp;            } // ZDC EP
    TVector2 getQ1VecZdcEast() const       { return v_mQ1ZdcShiftEast;     }
    TVector2 getQ1VecZdcWest() const       { return v_mQ1ZdcShiftWest;     }
    TVector2 getQ1VecZdcFull() const       { return v_mQ1ZdcShiftFull;     }

    int getFlagEpdSideEp() const           { return mFlagEpdSideEp;        } // EPD Side EP
    TVector2 getQ1VecEpdSideEast() const   { return v_mQ1EpdSideShiftEast; }
    TVector2 getQ1VecEpdSideWest() const   { return v_mQ1EpdSideShiftWest; }
    TVector2 getQ1VecEpdSideFull() const   { return v_mQ1EpdSideShiftFull; }

    int getFlagEpdGrp0Ep() const           { return mFlagEpdGrp0Ep;        } // EPD Grp0 EP
    TVector2 getQ1VecEpdGrp0East() const   { return v_mQ1EpdGrp0ShiftEast; }
    TVector2 getQ1VecEpdGrp0West() const   { return v_mQ1EpdGrp0ShiftWest; }
    TVector2 getQ1VecEpdGrp0Full() const   { return v_mQ1EpdGrp0ShiftFull; }
    int getFlagEpdGrp1Ep() const           { return mFlagEpdGrp1Ep;        } // EPD Grp1 EP
    TVector2 getQ1VecEpdGrp1East() const   { return v_mQ1EpdGrp1ShiftEast; }
    TVector2 getQ1VecEpdGrp1West() const   { return v_mQ1EpdGrp1ShiftWest; }
    TVector2 getQ1VecEpdGrp1Full() const   { return v_mQ1EpdGrp1ShiftFull; }

    int getFlagTpcEp() const               { return mFlagTpcEp;            } // TPC EP
    TVector2 getQ1VecTpcEast() const       { return v_mQ1TpcReCtrEast;     }
    TVector2 getQ1VecTpcWest() const       { return v_mQ1TpcReCtrWest;     }
    TVector2 getQ2VecTpcEast() const       { return v_mQ2TpcReCtrEast;     }
    TVector2 getQ2VecTpcWest() const       { return v_mQ2TpcReCtrWest;     }
    TVector2 getQ3VecTpcEast() const       { return v_mQ3TpcReCtrEast;     }
    TVector2 getQ3VecTpcWest() const       { return v_mQ3TpcReCtrWest;     }
    int getNumTrkReCtrEast() const         { return mNumTrkReCtrEast;      }
    int getNumTrkReCtrWest() const         { return mNumTrkReCtrWest;      }
    // -----------------------------------Number of Tracks----------------------------------------
    StPhiMesonTrack* createTrack()
    {
      if (mNumTracks == phiTracks->GetSize()) { phiTracks->Expand( mNumTracks + 10 ); }
      if (mNumTracks >= 10000)
      {
	Fatal( "StPhiMesonEvent::createTrack()", "ERROR: Too many tracks (>10000)!" );
	exit( 2 );
      }

      new((*phiTracks)[mNumTracks++]) StPhiMesonTrack;
      return (StPhiMesonTrack*)((*phiTracks)[mNumTracks - 1]);
    }

    void clearEvtHeader()
    {
      mRunId       = -1;
      mRunIdx      = -1;
      mEvtId       = -1;
      mRefMult     = -1;
      mNumTofMatch = -1;
      mCent9       = -1;
      mCent16      = -1;
      mRefWgt      = -1.0;
      mZDCx        = -1.0;
      mBBCx        = -1.0;
      mVzVpd       = -1.0;
      mPrimVtx.SetXYZ(-999.9,-999.9,-999.9);

      mFlagZdcEp = -1; // ZDC EP
      v_mQ1ZdcShiftEast.Set(0.0,0.0);
      v_mQ1ZdcShiftWest.Set(0.0,0.0);
      v_mQ1ZdcShiftFull.Set(0.0,0.0);

      mFlagEpdSideEp = -1; // EPD Side EP
      v_mQ1EpdSideShiftEast.Set(0.0,0.0);
      v_mQ1EpdSideShiftWest.Set(0.0,0.0);
      v_mQ1EpdSideShiftFull.Set(0.0,0.0);

      mFlagEpdGrp0Ep = -1; // EPD Grp0 EP
      v_mQ1EpdGrp0ShiftEast.Set(0.0,0.0);
      v_mQ1EpdGrp0ShiftWest.Set(0.0,0.0);
      v_mQ1EpdGrp0ShiftFull.Set(0.0,0.0);
      mFlagEpdGrp1Ep = -1; // EPD Grp1 EP
      v_mQ1EpdGrp1ShiftEast.Set(0.0,0.0);
      v_mQ1EpdGrp1ShiftWest.Set(0.0,0.0);
      v_mQ1EpdGrp1ShiftFull.Set(0.0,0.0);

      mFlagTpcEp = -1; // TPC EP
      v_mQ1TpcReCtrEast.Set(0.0,0.0);
      v_mQ1TpcReCtrWest.Set(0.0,0.0);
      v_mQ2TpcReCtrEast.Set(0.0,0.0);
      v_mQ2TpcReCtrWest.Set(0.0,0.0);
      v_mQ3TpcReCtrEast.Set(0.0,0.0);
      v_mQ3TpcReCtrWest.Set(0.0,0.0);
      mNumTrkReCtrEast = -1;
      mNumTrkReCtrWest = -1;
    }

    void clearTrackList()
    {
      mNumTracks = 0;
      phiTracks->Clear();
    }

    UShort_t getNumTracks() const
    {
      return mNumTracks;
    }

    StPhiMesonTrack* getTrack(UShort_t iTrack) const
    {
      return iTrack < mNumTracks ? (StPhiMesonTrack*)((*phiTracks)[iTrack]) : NULL;
    }

    ClassDef(StPhiMesonEvent,1)  // A simple event compiled of tracks
};

#endif // StPhiMesonEvent_h
