#ifndef StVecMesonUtility_h
#define StVecMesonUtility_h

#include "StMessMgr.h"
#include <map>

class StVecMesonUtility
{
  public:
    StVecMesonUtility(int energy);
    virtual ~StVecMesonUtility();

    void initRunIndex();
    bool read_in_runIndex();
    int findRunIndex(int runId);

  private:
    int mEnergy;
    std::map<int,int> map_runIndex;

    ClassDef(StVecMesonUtility,1)
};

#endif
