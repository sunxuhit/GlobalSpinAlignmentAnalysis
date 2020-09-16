#include "StRoot/StEventPlaneMaker/StEventPlaneUtility.h"
#include "StRoot/StEventPlaneMaker/StEventPlaneCons.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TMath.h"

ClassImp(StEventPlaneUtility)

//---------------------------------------------------------------------------------

StEventPlaneUtility::StEventPlaneUtility(int energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StEventPlaneUtility::~StEventPlaneUtility()
{
  /* */
}

//---------------------------------------------------------------------------------


void StEventPlaneUtility::initRunIndex()
{
  map_runIndex.clear();
  bool isOpen_runIndex = read_in_runIndex();
  if(isOpen_runIndex) std::cout << "Run Index read in!" << std::endl;
}

int StEventPlaneUtility::findRunIndex(int runId)
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
    // std::cout << "StEventPlaneUtility -> could not find in full run list! & send signal to kill the run!" << std::endl;
    return -999;
  }
  else
  {
    // std::cout << "StEventPlaneUtility -> runId: " << it_runId->first << " => runIndex: " << it_runId->second << std::endl;
    return it_runId->second;
  }

  return -999;
}

bool StEventPlaneUtility::read_in_runIndex()
{
  // std::string inputfile = Form("StRoot/StEventPlaneMaker/runIndex_%s.txt",recoEP::mBeamEnergy[mEnergy].c_str());
  std::string inputfile = Form("StRoot/StEventPlaneUtility/RunIndex/runIndex_%s.txt",recoEP::mBeamEnergy[mEnergy].c_str());
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
