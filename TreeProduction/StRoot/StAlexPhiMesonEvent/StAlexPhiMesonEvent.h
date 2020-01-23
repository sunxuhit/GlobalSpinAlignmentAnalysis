#ifndef __STALEXPHIMESONEVENT_H__
#define __STALEXPHIMESONEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include "TLorentzVector.h"

// A. Schmah 19.12.2011

class StAlexPhiMesonTrack : public TObject
{
  private:
    // Track properties
    // TrackA: for phi => K_plus
    // TrackB: for phi => K_minus
    Float_t mMass2A; // mass2 of track A
    Float_t mMass2B;
    Float_t mNSigKaonA; // nsigma dE/dx of particle A
    Float_t mNSigKaonB;
    Float_t mDcaA; // distance of closest approach of particle A * charge
    Float_t mDcaB;
    TLorentzVector mTrackA; // Lorentz Vector for track A (px, py, pz, mass2)
    TLorentzVector mTrackB;
    Int_t mFlagA; // Flag for Event A: 0 for Same Event, others for Mixed Event
    Int_t mFlagB;

  public:
    StAlexPhiMesonTrack() :
      mMass2A(-999.9),mMass2B(-999.9),mNSigKaonA(-999.9),mNSigKaonB(-999.9),mDcaA(-999.9),mDcaB(-999.9),mTrackA(0,0,0,0),mTrackB(0,0,0,0),mFlagA(-1),mFlagB(-1)
  {
  }
    ~StAlexPhiMesonTrack() {}

    // setters
    void setMass2A(Float_t f)                  { mMass2A = f;       }
    void setMass2B(Float_t f)                  { mMass2B = f;       }
    void setNSigKaonA(Float_t f)               { mNSigKaonA = f;    }
    void setNSigKaonB(Float_t f)               { mNSigKaonB = f;    }
    void setDcaA(Float_t f)                    { mDcaA = f;         }
    void setDcaB(Float_t f)                    { mDcaB = f;         }
    void setTrackA(TLorentzVector f)           { mTrackA = f;       }
    void setTrackB(TLorentzVector f)           { mTrackB = f;       }
    void setFlagA(Int_t f)                     { mFlagA = f;        }
    void setFlagB(Int_t f)                     { mFlagB = f;        }

    // getters
    Float_t getMass2A() const                  { return mMass2A;    }
    Float_t getMass2B() const                  { return mMass2B;    }
    Float_t getNSigKaonA() const               { return mNSigKaonA; }
    Float_t getNSigKaonB() const               { return mNSigKaonB; }
    Float_t getDcaA() const                    { return mDcaA;      }
    Float_t getDcaB() const                    { return mDcaB;      }
    TLorentzVector getTrackA() const           { return mTrackA;    }
    TLorentzVector getTrackB() const           { return mTrackB;    }
    Int_t getFlagA() const                     { return mFlagA;     }
    Int_t getFlagB() const                     { return mFlagB;     }

    ClassDef(StAlexPhiMesonTrack,1)  // A simple track of a particle
};

class StAlexPhiMesonEvent : public TObject
{
  private:
    StThreeVectorF mPrimaryvertex;
    Int_t   mRunId;
    Int_t   mEventId;
    Int_t   mRefMult;
    Int_t   mCentrality;
    Int_t   mN_prim;
    Int_t   mN_non_prim;
    Int_t   mN_Tof_match;

    Float_t mZDCx;
    Float_t mBBCx;
    Float_t mVzVpd;

    UShort_t      fNumTracks;

    TVector2 mQ2East;
    TVector2 mQ2West;
    TVector2 mQ2Full;
    Int_t   mNumTrackEast;
    Int_t   mNumTrackWest;
    Int_t   mNumTrackFull;
    Int_t   mNumTrackFullEast;
    Int_t   mNumTrackFullWest;

    TClonesArray* fTracks;      //->

  public:
    StAlexPhiMesonEvent() :
      mPrimaryvertex(-1.0,-1.0,-1.0),mRunId(-1),mEventId(-1),mRefMult(-1),mCentrality(-1),mN_prim(-1),mN_non_prim(-1),mN_Tof_match(-1),mZDCx(-1),mBBCx(-1),mVzVpd(-1),fNumTracks(0)
  {
    mQ2East.Set(-999.9,-999.9); // QVector2 East
    mQ2West.Set(-999.9,-999.9); // QVector2 West
    mQ2Full.Set(-999.9,-999.9); // QVector2 West

    mNumTrackEast = 0;
    mNumTrackWest = 0;
    mNumTrackFull = 0;
    mNumTrackFullEast = 0;
    mNumTrackFullWest = 0;

    fTracks      = new TClonesArray( "StAlexPhiMesonTrack", 10 );
  }
    ~StAlexPhiMesonEvent()
    {
      delete fTracks;
      fTracks = NULL;
    }

    void       setPrimaryVertex(StThreeVectorF r)      { mPrimaryvertex = r;       }
    StThreeVectorF    getPrimaryVertex() const         { return mPrimaryvertex;    }

    void       setRunId(Int_t  r)                      { mRunId = r;               }
    Int_t      getRunId() const                        { return mRunId;            }

    void       setEventId(Int_t  r)                    { mEventId = r;             }
    Int_t      getEventId() const                      { return mEventId;          }

    void       setRefMult(Int_t r)                     { mRefMult = r;             }
    Int_t      getRefMult() const                      { return mRefMult;          }

    void       setCentrality(Int_t r)                  { mCentrality = r;          }
    Int_t      getCentrality() const                   { return mCentrality;       }

    void       setN_prim(Int_t r)                      { mN_prim = r;              }
    Int_t      getN_prim() const                       { return mN_prim;           }

    void       setN_non_prim(Int_t r)                  { mN_non_prim = r;          }
    Int_t      getN_non_prim() const                   { return mN_non_prim;       }

    void       setN_Tof_match(Int_t r)                 { mN_Tof_match = r;         }
    Int_t      getN_Tof_match() const                  { return mN_Tof_match;      }


    void       setZDCx(Float_t r)                      { mZDCx = r;                }
    Float_t    getZDCx() const                         { return mZDCx;             }

    void       setBBCx(Float_t r)                      { mBBCx = r;                }
    Float_t    getBBCx() const                         { return mBBCx;             }

    void       setVzVpd(Float_t r)                     { mVzVpd = r;               }
    Float_t    getVzVpd() const                        { return mVzVpd;            }

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
    void       setNumTrackEast(Int_t r)                { mNumTrackEast = r;        }
    Int_t      getNumTrackEast() const                 { return mNumTrackEast;     }
    // West
    void       setNumTrackWest(Int_t r)                { mNumTrackWest = r;        }
    Int_t      getNumTrackWest() const                 { return mNumTrackWest;     }
    // Full 
    void       setNumTrackFull(Int_t r)                { mNumTrackFull = r;        }
    Int_t      getNumTrackFull() const                 { return mNumTrackFull;     }
    // Full East
    void       setNumTrackFullEast(Int_t r)            { mNumTrackFullEast = r;    }
    Int_t      getNumTrackFullEast() const             { return mNumTrackFullEast; }
    // Full West
    void       setNumTrackFullWest(Int_t r)            { mNumTrackFullWest = r;    }
    Int_t      getNumTrackFullWest() const             { return mNumTrackFullWest; }
    // -----------------------------------Number of Tracks----------------------------------------
    StAlexPhiMesonTrack* createTrack()
    {
      if (fNumTracks == fTracks->GetSize())
	fTracks->Expand( fNumTracks + 10 );
      if (fNumTracks >= 10000)
      {
	Fatal( "StAlexPhiMesonEvent::createTrack()", "ERROR: Too many tracks (>10000)!" );
	exit( 2 );
      }

      new((*fTracks)[fNumTracks++]) StAlexPhiMesonTrack;
      return (StAlexPhiMesonTrack*)((*fTracks)[fNumTracks - 1]);
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
    StAlexPhiMesonTrack* getTrack(UShort_t i) const
    {
      return i < fNumTracks ? (StAlexPhiMesonTrack*)((*fTracks)[i]) : NULL;
    }

    ClassDef(StAlexPhiMesonEvent,1)  // A simple event compiled of tracks
};


#endif // __STALEXPHIMESONEVENT_H__
