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
  map_runIndex.clear();
  bool isOpen_runIndex = read_in_runIndex();
  if(isOpen_runIndex) std::cout << "Run Index read in!" << std::endl;
}

int StRunQAUtility::findRunIndex(int runId)
{
  // print map_runIndex content:
  /*
     for (std::map<int,int>::iterator it=map_runIndex.begin(); it!=map_runIndex.end(); ++it)
     {
     std::cout << it->first << " => " << it->second << '\n';
     }
     */

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

bool StRunQAUtility::read_in_runIndex()
{
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
  std::cout << "reading run Index: " << std::endl;
  while (file_runIndex >> temp_runId >> temp_runIndex)
  {
    // std::cout << "runId = " << temp_runId << ", runIndex = " << temp_runIndex << std::endl;
    map_runIndex[temp_runId] = temp_runIndex;
  }
  file_runIndex.close();

  return true;
}
