#ifndef StPhiMesonEvent_h
#define StPhiMesonEvent_h

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

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
    TLorentzVector mTrackA; // Lorentz Vector for track A (px, py, pz, mass2)
    TLorentzVector mTrackB;
    int mFlagA; // Flag for Event A: 0 for Same Event, others for Mixed Event
    int mFlagB;

  public:
    StPhiMesonTrack() :
      mMass2A(-999.9),mMass2B(-999.9),mNSigKaonA(-999.9),mNSigKaonB(-999.9),mDcaA(-999.9),mDcaB(-999.9),mTrackA(0,0,0,0),mTrackB(0,0,0,0),mFlagA(-1),mFlagB(-1)
  {
  }
    ~StPhiMesonTrack() {}

    // setters
    void setMass2A(double f)           { mMass2A    = f; }
    void setMass2B(double f)           { mMass2B    = f; }
    void setNSigKaonA(double f)        { mNSigKaonA = f; }
    void setNSigKaonB(double f)        { mNSigKaonB = f; }
    void setDcaA(double f)             { mDcaA      = f; }
    void setDcaB(double f)             { mDcaB      = f; }
    void setTrackA(TLorentzVector f)   { mTrackA    = f; }
    void setTrackB(TLorentzVector f)   { mTrackB    = f; }
    void setFlagA(int f)               { mFlagA     = f; }
    void setFlagB(int f)               { mFlagB     = f; }

    // getters
    double getMass2A() const           { return mMass2A;    }
    double getMass2B() const           { return mMass2B;    }
    double getNSigKaonA() const        { return mNSigKaonA; }
    double getNSigKaonB() const        { return mNSigKaonB; }
    double getDcaA() const             { return mDcaA;      }
    double getDcaB() const             { return mDcaB;      }
    TLorentzVector getTrackA() const   { return mTrackA;    }
    TLorentzVector getTrackB() const   { return mTrackB;    }
    int getFlagA() const               { return mFlagA;     }
    int getFlagB() const               { return mFlagB;     }

    ClassDef(StPhiMesonTrack,1)  // A simple track of a particle
};

class StPhiMesonEvent : public TObject
{
  private:
    TVector3 mPrimVtx;
    int   mRunId;
    int   mEventId;
    int   mRefMult;
    int   mCent9;
    int   mCent16;
    int   mNumTofMatch;
    double mRefWgt;
    double mZDCx;
    double mBBCx;
    double mVzVpd;

    unsigned short fNumTracks;

    TVector2 mQ2TpcReCtrEast; // QVector
    TVector2 mQ2TpcReCtrWest;
    int   mNumTpcTrackEast;
    int   mNumTpcTrackWest;

    TClonesArray* fTracks; //->

  public:
    StPhiMesonEvent() :
      mPrimVtx(-1.0,-1.0,-1.0),mRunId(-1),mEventId(-1),mRefMult(-1),mCent9(-1),mN_prim(-1),mN_non_prim(-1),mN_Tof_match(-1),mZDCx(-1),mBBCx(-1),mVzVpd(-1),fNumTracks(0)
  {
    mQ2East.Set(-999.9,-999.9); // QVector2 East
    mQ2West.Set(-999.9,-999.9); // QVector2 West
    mQ2Full.Set(-999.9,-999.9); // QVector2 West

    mNumTrackEast = 0;
    mNumTrackWest = 0;
    mNumTrackFull = 0;
    mNumTrackFullEast = 0;
    mNumTrackFullWest = 0;

    fTracks      = new TClonesArray( "StPhiMesonTrack", 10 );
  }
    ~StPhiMesonEvent()
    {
      delete fTracks;
      fTracks = NULL;
    }

    void       setPrimaryVertex(TVector3 r)      { mPrimVtx = r;       }
    TVector3 getPrimaryVertex() const         { return mPrimVtx;    }

    void       setRunId(int  r)                      { mRunId = r;               }
    int      getRunId() const                        { return mRunId;            }

    void       setEventId(int  r)                    { mEventId = r;             }
    int      getEventId() const                      { return mEventId;          }

    void       setRefMult(int r)                     { mRefMult = r;             }
    int      getRefMult() const                      { return mRefMult;          }

    void       setCentrality(int r)                  { mCent9 = r;          }
    int      getCentrality() const                   { return mCent9;       }

    void       setN_prim(int r)                      { mN_prim = r;              }
    int      getN_prim() const                       { return mN_prim;           }

    void       setN_non_prim(int r)                  { mN_non_prim = r;          }
    int      getN_non_prim() const                   { return mN_non_prim;       }

    void       setN_Tof_match(int r)                 { mN_Tof_match = r;         }
    int      getN_Tof_match() const                  { return mN_Tof_match;      }


    void       setZDCx(double r)                      { mZDCx = r;                }
    double    getZDCx() const                         { return mZDCx;             }

    void       setBBCx(double r)                      { mBBCx = r;                }
    double    getBBCx() const                         { return mBBCx;             }

    void       setVzVpd(double r)                     { mVzVpd = r;               }
    double    getVzVpd() const                        { return mVzVpd;            }

    // ---------------------------------------QVector---------------------------------------------
    // QVector2 East
    void       setQ2East(TVector2 r)                   { mQ2East = r;              }
    TVector2   getQ2East() const                       { return mQ2East;           }
    // QVector2 West
    void       setQ2West(TVector2 r)                   { mQ2West = r;              }
    TVector2   getQ2West() const                       { return mQ2West;           }
    // QVector2 Full 
    void       setQ2Full(TVector2 r)                   { mQ2Full = r;              }
    TVector2   getQ2Full() const                       { return mQ2Full;           }
    // ---------------------------------------QVector---------------------------------------------

    // -----------------------------------Number of Tracks----------------------------------------
    // East
    void       setNumTrackEast(int r)                { mNumTrackEast = r;        }
    int      getNumTrackEast() const                 { return mNumTrackEast;     }
    // West
    void       setNumTrackWest(int r)                { mNumTrackWest = r;        }
    int      getNumTrackWest() const                 { return mNumTrackWest;     }
    // Full 
    void       setNumTrackFull(int r)                { mNumTrackFull = r;        }
    int      getNumTrackFull() const                 { return mNumTrackFull;     }
    // Full East
    void       setNumTrackFullEast(int r)            { mNumTrackFullEast = r;    }
    int      getNumTrackFullEast() const             { return mNumTrackFullEast; }
    // Full West
    void       setNumTrackFullWest(int r)            { mNumTrackFullWest = r;    }
    int      getNumTrackFullWest() const             { return mNumTrackFullWest; }
    // -----------------------------------Number of Tracks----------------------------------------
    StPhiMesonTrack* createTrack()
    {
      if (fNumTracks == fTracks->GetSize())
	fTracks->Expand( fNumTracks + 10 );
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
      fNumTracks   = 0;
      fTracks      ->Clear();
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


#endif // __STALEXPHIMESONEVENT_H__
