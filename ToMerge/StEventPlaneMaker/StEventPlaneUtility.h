#ifndef StEventPlaneUtility_h
#define StEventPlaneUtility_h

#include "StMessMgr.h"
#include <map>

class StEventPlaneUtility
{
  public:
    StEventPlaneUtility(int energy);
    virtual ~StEventPlaneUtility();

    void initRunIndex();
    bool read_in_runIndex();
    int findRunIndex(int runId);

  private:
    int mEnergy;
    std::map<int,int> map_runIndex;

    ClassDef(StEventPlaneUtility,1)
};

#endif
