#ifndef StRunQAUtility_h
#define StRunQAUtility_h

#include "StMessMgr.h"
#include <map>
#include <vector>

class StRunQAUtility
{
  public:
    StRunQAUtility(int energy);
    virtual ~StRunQAUtility();

    void initRunIndex();

    bool read_in_runIndex();
    int findRunIndex(int runId);

    bool read_in_badRunList();
    bool isBadRun(int runId);

  private:
    int mEnergy;
    std::map<int,int> map_runIndex;
    std::vector<int> vec_badRunId;

    ClassDef(StRunQAUtility,1)
};

#endif
