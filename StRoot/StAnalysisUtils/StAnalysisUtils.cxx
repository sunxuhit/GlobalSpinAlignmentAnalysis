#include <iostream>

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StMessMgr.h"

#include "Utility/include/StSpinAlignmentCons.h"
#include "StRoot/StAnalysisUtils/StAnalysisUtils.h"

ClassImp(StAnalysisUtils)

//---------------------------------------------------------------------------------

StAnalysisUtils::StAnalysisUtils(int beamType)
{
  mType = beamType;
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
  map_runIndex.clear();

  std::string inputfile = Form("Utility/RunIndex/runIndex_%s.txt",globCons::mBeamType[mType].c_str());
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
    map_runIndex[temp_runId] = temp_runIndex;
  }
  file_runIndex.close();

  // print map_runIndex content:
  // for (std::map<int,int>::iterator it=map_runIndex.begin(); it!=map_runIndex.end(); ++it)
  // {
  //   std::cout << it->first << " => " << it->second << '\n';
  // }

  return true;
}

int StAnalysisUtils::findRunIndex(int runId)
{
  // print map_runIndex content:
  // for (std::map<int,int>::iterator it=map_runIndex.begin(); it!=map_runIndex.end(); ++it)
  // {
  //   std::cout << it->first << " => " << it->second << '\n';
  // }

  std::map<int,int>::iterator it_runId = map_runIndex.find(runId);
  if(it_runId == map_runIndex.end())
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
  vec_badRunId.clear();

  std::string inputfile = Form("Utility/RunIndex/badRunList_%s.txt",globCons::mBeamType[mType].c_str());
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
    vec_badRunId.push_back(runId);
  }
  file_badRunList.close();

  // print vec_badRunId content:
  // for (std::vector<int>::iterator it=vec_badRunId.begin(); it!=vec_badRunId.end(); ++it)
  // {
  //   std::cout << (*it) << std::endl;;
  // }

  return true;
}

bool StAnalysisUtils::isBadRun(int runId)
{
  // print vec_badRunId content:
  // for (std::vector<int>::iterator it=vec_badRunId.begin(); it!=vec_badRunId.end(); ++it)
  // {
  //   std::cout << (*it) << std::endl;;
  // }

  std::vector<int>::iterator it_runId = std::find(vec_badRunId.begin(), vec_badRunId.end(), runId);

  return ( it_runId != vec_badRunId.end() ); // true if can be found in bad run list
}

//---------------------------------------------------------------------------------
float StAnalysisUtils::getBeta(StPicoDst *picoDst, int i_track)
{
  float beta = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    beta = tofTrack->btofBeta();
  }

  return beta;
}

float StAnalysisUtils::getPrimaryMass2(StPicoDst *picoDst, int i_track)
{
  float Mass2 = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    float beta = tofTrack->btofBeta();
    // float Momentum = picoTrack->pMom().mag(); // primary momentum for 54GeV_2017
    // float Momentum = picoTrack->pMom().Mag(); // primary momentum for 27GeV_2018
    TVector3 primMom; // temp fix for StThreeVectorF & TVector3
    float primPx    = picoTrack->pMom().x(); // x works for both TVector3 and StThreeVectorF
    float primPy    = picoTrack->pMom().y();
    float primPz    = picoTrack->pMom().z();
    primMom.SetXYZ(primPx,primPy,primPz);
    float Momentum = primMom.Mag(); // primary momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && beta != 0)
    {
      Mass2 = Momentum*Momentum*(1.0/(beta*beta) - 1.0);
    }
  }

  return Mass2;
}

float StAnalysisUtils::getGlobalMass2(StPicoDst *picoDst, int i_track)
{
  float Mass2 = -999.9;
  StPicoTrack *picoTrack = (StPicoTrack*)picoDst->track(i_track); // return ith track
  int tofIndex = picoTrack->bTofPidTraitsIndex(); // return ToF PID traits
  if(tofIndex >= 0)
  {
    StPicoBTofPidTraits *tofTrack = picoDst->btofPidTraits(tofIndex);
    float beta = tofTrack->btofBeta();
    // float Momentum = picoTrack->gMom().mag(); // global momentum for 54GeV_2017
    // float Momentum = picoTrack->gMom().Mag(); // global momentum for 27GeV_2018
    TVector3 globMom; // temp fix for StThreeVectorF & TVector3
    float globPx     = picoTrack->gMom().x(); // x works for both TVector3 and StThreeVectorF
    float globPy     = picoTrack->gMom().y();
    float globPz     = picoTrack->gMom().z();
    globMom.SetXYZ(globPx,globPy,globPz);
    float Momentum = globMom.Mag(); // global momentum

    if(tofTrack->btofMatchFlag() > 0 && tofTrack->btof() != 0 && beta != 0)
    {
      Mass2 = Momentum*Momentum*(1.0/(beta*beta) - 1.0);
    }
  }

  return Mass2;
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

  return -1;
}
