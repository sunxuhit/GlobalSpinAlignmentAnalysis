#include <iostream>

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StMessMgr.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"

ClassImp(StAnalysisUtils)

//---------------------------------------------------------------------------------

StAnalysisUtils::StAnalysisUtils(int beamType) : mType(beamType)
{
  // mType = beamType;
}

//---------------------------------------------------------------------------------

StAnalysisUtils::~StAnalysisUtils()
{
  /* */
}

//---------------------------------------------------------------------------------

void StAnalysisUtils::initRunIndex()
{
  bool isOpen_runIndex = readRunIndex(); // read in runId vs. runIndex
  if(isOpen_runIndex) std::cout << "Run Index read in!" << std::endl;

  bool isOpen_badRunList = readBadRunList(); // read in Bad Run List
  if(isOpen_badRunList) std::cout << "Bad Run List read in!" << std::endl;
}

bool StAnalysisUtils::readRunIndex()
{
  map_mRunIndex.clear();

  std::string inputfile = Form("Utility/RunIndex/%s/runIndex_%s.txt",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  std::cout << "inputfile = " << inputfile.c_str() << std::endl;
  std::ifstream file_runIndex ( inputfile.c_str() );
  if ( !file_runIndex.is_open() )
  {
    std::cout << "Abort. Fail to read in run Index file: " << inputfile << std::endl;
    return false;
  }

  int temp_runId = 0, temp_runIndex = 0;
  std::cout << "reading run Index: " << inputfile.c_str() << std::endl;
  while (file_runIndex >> temp_runId >> temp_runIndex)
  {
    // std::cout << "runId = " << temp_runId << ", runIndex = " << temp_runIndex << std::endl;
    map_mRunIndex[temp_runId] = temp_runIndex;
  }
  file_runIndex.close();

  // print map_mRunIndex content:
  // for (std::map<int,int>::iterator it=map_mRunIndex.begin(); it!=map_mRunIndex.end(); ++it)
  // {
  //   std::cout << it->first << " => " << it->second << '\n';
  // }

  return true;
}

int StAnalysisUtils::findRunIndex(int runId)
{
  // print map_mRunIndex content:
  // for (std::map<int,int>::iterator it=map_mRunIndex.begin(); it!=map_mRunIndex.end(); ++it)
  // {
  //   std::cout << it->first << " => " << it->second << '\n';
  // }

  std::map<int,int>::iterator it_runId = map_mRunIndex.find(runId);
  if(it_runId == map_mRunIndex.end())
  {
    // std::cout << "StAnalysisUtils -> could not find in full run list! & send signal to kill the run!" << std::endl;
    return -999;
  }
  else
  {
    // std::cout << "StAnalysisUtils -> runId: " << it_runId->first << " => runIndex: " << it_runId->second << std::endl;
    return it_runId->second;
  }

  return -999;
}

bool StAnalysisUtils::readBadRunList()
{
  vec_mBadRunId.clear();

  std::string inputfile = Form("Utility/FileList/%s/badRunList_%s.txt",globCons::str_mBeamType[mType].c_str(),globCons::str_mBeamType[mType].c_str());
  std::cout << "inputfile = " << inputfile.c_str() << std::endl;
  std::ifstream file_badRunList ( inputfile.c_str() );
  if ( !file_badRunList.is_open() )
  {
    std::cout << "Abort. Fail to read in bad Run List file: " << inputfile << std::endl;
    return false;
  }

  int runId = 0;
  std::cout << "reading bad Bun List: " << inputfile.c_str() << std::endl;
  while (file_badRunList >> runId)
  {
    vec_mBadRunId.push_back(runId);
  }
  file_badRunList.close();

  // print vec_mBadRunId content:
  // for (std::vector<int>::iterator it=vec_mBadRunId.begin(); it!=vec_mBadRunId.end(); ++it)
  // {
  //   std::cout << (*it) << std::endl;;
  // }

  return true;
}

bool StAnalysisUtils::isBadRun(int runId)
{
  // print vec_mBadRunId content:
  // for (std::vector<int>::iterator it=vec_mBadRunId.begin(); it!=vec_mBadRunId.end(); ++it)
  // {
  //   std::cout << (*it) << std::endl;;
  // }

  std::vector<int>::iterator it_runId = std::find(vec_mBadRunId.begin(), vec_mBadRunId.end(), runId);

  return ( it_runId != vec_mBadRunId.end() ); // true if can be found in bad run list
}

//---------------------------------------------------------------------------------
double StAnalysisUtils::getBeta(StPicoDst *picoDst, int i_track)
{
  double beta = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    beta = tofTrack->btofBeta();
  }

  return beta;
}

double StAnalysisUtils::getPrimaryMass2(StPicoDst *picoDst, int i_track)
{
  double mass2 = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    double beta = tofTrack->btofBeta();
    const TVector3 primMom = picoTrack->pMom(); // primary Momentum
    double primMomMag = primMom.Mag(); // primary momentum magnitude

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && beta != 0)
    {
      mass2 = primMomMag*primMomMag*(1.0/(beta*beta) - 1.0);
    }
  }

  return mass2;
}

double StAnalysisUtils::getGlobalMass2(StPicoDst *picoDst, int i_track)
{
  double mass2 = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    double beta = tofTrack->btofBeta();
    const TVector3 globMom = picoTrack->gMom(); // global Momentum
    double globMomMag = globMom.Mag(); // global momentum magnitude

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && beta != 0)
    {
      mass2 = globMomMag*globMomMag*(1.0/(beta*beta) - 1.0);
    }
  }

  return mass2;
}

int StAnalysisUtils::getTriggerBin(StPicoEvent *picoEvent)
{
  // std::cout << "year: " << picoEvent->year() << std::endl;
  // std::cout << "day: " << picoEvent->day() << std::endl;
  // std::cout << "triggerIds: " << picoEvent->triggerIds()[0] << std::endl;
  if( (mType == 0 || mType == 1) && globCons::mBeamYear[mType] == picoEvent->year() )
  { // ZrZr200GeV_2018 || RuRu200GeV_2018
    if( picoEvent->isTrigger(600001) ) return 0; // VPDMB-30
    if( picoEvent->isTrigger(600011) ) return 1; // VPDMB-30
    if( picoEvent->isTrigger(600021) ) return 2; // VPDMB-30
    if( picoEvent->isTrigger(600031) ) return 3; // VPDMB-30
  }
  if( (mType == 2) && globCons::mBeamYear[mType] == picoEvent->year() )
  { // Fxt3p85GeV_2018
    if( picoEvent->isTrigger(620052) ) return 0; // bbce_tofmult1
  }

  return -1;
}

int StAnalysisUtils::getVzBin(double vz)
{
  if(vz >= anaUtils::mVzMin[mType] && vz < anaUtils::mVzCtr[mType]) 
  {
    return 0;
  }
  if(vz >= anaUtils::mVzCtr[mType] && vz <= anaUtils::mVzMax[mType]) 
  {
    return 1;
  }

  return -1;
}
