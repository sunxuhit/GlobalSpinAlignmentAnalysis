#ifndef StRunQAUtility_h
#define StRunQAUtility_h

#include "StMessMgr.h"
#include <map>
#include <vector>

class StPicoDst;
class StPicoEvent;
class StPicoTrack;

class StRunQAUtility
{
  public:
    StRunQAUtility(int beamType);
    virtual ~StRunQAUtility();

    void initRunIndex();

    bool read_in_runIndex();
    int findRunIndex(int runId);

    bool read_in_badRunList();
    bool isBadRun(int runId);

    float getBeta(StPicoDst*, int); // return beta of i-th track (tof || -999)
    float getPrimaryMass2(StPicoDst*, int); // return m^2 of i-th track (primary || -999)
    float getGlobalMass2(StPicoDst*, int); // return m^2 of i-th track (global || -999)
    int getTriggerBin(StPicoEvent*); // return trigger bin for event QA

  private:
    int mType;
    std::map<int,int> map_runIndex;
    std::vector<int> vec_badRunId;

    ClassDef(StRunQAUtility,1)
};

#endif
