#ifndef StAnalysisUtils_h
#define StAnalysisUtils_h

#include <map>
#include <vector>

class StPicoDst;
class StPicoEvent;
class StPicoTrack;

class StAnalysisUtils
{
  public:
    StAnalysisUtils(int beamType);
    virtual ~StAnalysisUtils();

    void initRunIndex();

    bool readRunIndex();
    int findRunIndex(int runId);

    bool readBadRunList();
    bool isBadRun(int runId);

    double getBeta(StPicoDst*, int); // return beta of i-th track (tof || -999)
    double getPrimaryMass2(StPicoDst*, int); // return m^2 of i-th track (primary || -999)
    double getGlobalMass2(StPicoDst*, int); // return m^2 of i-th track (global || -999)
    int getTriggerBin(StPicoEvent*); // return trigger bin for event QA
    int getVzBin(double); // return vz bin

  private:
    const int mType;
    std::map<int,int> map_mRunIndex;
    std::vector<int> vec_mBadRunId;

    ClassDef(StAnalysisUtils,1)
};

#endif
