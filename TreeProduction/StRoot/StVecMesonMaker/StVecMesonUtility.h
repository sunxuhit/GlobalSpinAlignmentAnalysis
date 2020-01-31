#ifndef StVecMesonUtility_h
#define StVecMesonUtility_h

#include "StMessMgr.h"
#include <map>

class StVecMesonUtility
{
  public:
    StVecMesonUtility(int energy);
    virtual ~StVecMesonUtility();

    void Init_RunIndex();
    bool read_in_runIndex();
    int find_runIndex(int runId);

  private:
    int mEnergy;
    std::map<int,int> map_runIndex;

    ClassDef(StVecMesonUtility,1)
};

#endif
