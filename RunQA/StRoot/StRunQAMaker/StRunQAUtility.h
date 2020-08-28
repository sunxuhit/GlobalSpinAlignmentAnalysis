#ifndef StRunQAUtility_h
#define StRunQAUtility_h

#include "StMessMgr.h"
#include <map>

class StRunQAUtility
{
  public:
    StRunQAUtility(int energy);
    virtual ~StRunQAUtility();

    void initRunIndex();
    bool read_in_runIndex();
    int findRunIndex(int runId);

  private:
    int mEnergy;
    std::map<int,int> map_runIndex;

    ClassDef(StRunQAUtility,1)
};

#endif
