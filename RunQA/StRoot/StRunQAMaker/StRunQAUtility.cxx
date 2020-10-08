#include "StRoot/StRunQAMaker/StRunQAUtility.h"
#include "StRoot/StRunQAMaker/StRunQACons.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TMath.h"

ClassImp(StRunQAUtility)

//---------------------------------------------------------------------------------

StRunQAUtility::StRunQAUtility(int energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StRunQAUtility::~StRunQAUtility()
{
  /* */
}

//---------------------------------------------------------------------------------


void StRunQAUtility::initRunIndex()
{
  bool isOpen_runIndex = read_in_runIndex(); // read in runId vs. runIndex
  if(isOpen_runIndex) std::cout << "Run Index read in!" << std::endl;

  bool isOpen_badRunList = read_in_badRunList(); // read in Bad Run List
  if(isOpen_badRunList) std::cout << "Bad Run List read in!" << std::endl;
}

bool StRunQAUtility::read_in_runIndex()
{
  map_runIndex.clear();

  // std::string inputfile = Form("StRoot/StRunQAMaker/runIndex_%s.txt",runQA::mBeamEnergy[mEnergy].c_str());
  std::string inputfile = Form("StRoot/StRunQAUtility/RunIndex/runIndex_%s.txt",runQA::mBeamEnergy[mEnergy].c_str());
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

int StRunQAUtility::findRunIndex(int runId)
{
  // print map_runIndex content:
  // for (std::map<int,int>::iterator it=map_runIndex.begin(); it!=map_runIndex.end(); ++it)
  // {
  //   std::cout << it->first << " => " << it->second << '\n';
  // }

  std::map<int,int>::iterator it_runId = map_runIndex.find(runId);
  if(it_runId == map_runIndex.end())
  {
    // std::cout << "StRunQAUtility -> could not find in full run list! & send signal to kill the run!" << std::endl;
    return -999;
  }
  else
  {
    // std::cout << "StRunQAUtility -> runId: " << it_runId->first << " => runIndex: " << it_runId->second << std::endl;
    return it_runId->second;
  }

  return -999;
}

bool StRunQAUtility::read_in_badRunList()
{
  vec_badRunId.clear();

  // std::string inputfile = Form("StRoot/StRunQAMaker/runIndex_%s.txt",runQA::mBeamEnergy[mEnergy].c_str());
  std::string inputfile = Form("StRoot/StRunQAUtility/RunIndex/badRunList_%s.txt",runQA::mBeamEnergy[mEnergy].c_str());
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

bool StRunQAUtility::isBadRun(int runId)
{
  // print vec_badRunId content:
  // for (std::vector<int>::iterator it=vec_badRunId.begin(); it!=vec_badRunId.end(); ++it)
  // {
  //   std::cout << (*it) << std::endl;;
  // }

  std::vector<int>::iterator it_runId = std::find(vec_badRunId.begin(), vec_badRunId.end(), runId);

  return ( it_runId != vec_badRunId.end() ); // true if can be found in bad run list
}
