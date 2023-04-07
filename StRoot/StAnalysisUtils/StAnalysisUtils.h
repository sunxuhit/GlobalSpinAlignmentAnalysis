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

    int getTriggerBin(StPicoEvent *picoEvent); // return trigger bin for event QA
    int getVzBin(double vz); // return vz bin
    double getVxReCtr(double vx); // return vx - anaUtils::mVxCtr[mType]
    double getVyReCtr(double vy); // return vy - anaUtils::mVyCtr[mType]
    double getBeta(StPicoDst *picoDst, int iTrack); // return beta of i-th track (tof || -999)
    double getPrimMass2(StPicoDst *picoDst, int iTrack); // return m^2 of i-th track (primary || -999)
    double getGlobMass2(StPicoDst *picoDst, int iTrack); // return m^2 of i-th track (global || -999)
    double getRapidityLab(StPicoTrack *picoTrack, int pid); // return particle rapidity in Lab frame
    double getRapidityCMS(double rapLab); // return particle rapidity in CMS: rapLab - anaUtils::mRapCtrM[mType]
    double calcNSigmaZ(int charge, double mass, double mom, double dEdx); // calculate nSigmaZ for deuteron

  private:
    const int mType;
    std::map<int,int> map_mRunIndex;
    std::vector<int> vec_mBadRunId;

    ClassDef(StAnalysisUtils,1)
};

#endif
