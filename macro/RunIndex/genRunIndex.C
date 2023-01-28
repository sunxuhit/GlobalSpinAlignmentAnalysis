#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <iterator>
#include "../../Utility/include/StSpinAlignmentCons.h"

int genRunIndex(int beamType = 0)
{
  const int numOfRuns = globCons::mMaxRunIndex[beamType];
  int runId[numOfRuns];
  int runIndex[numOfRuns];
  for(int i_run = 0; i_run < numOfRuns; ++i_run)
  {
    runId[i_run] = -999;
    runIndex[i_run] = -999;
  }

  // const string mBeamSpec[2] = {"ZrZr200GeV_2018","RuRu200GeV_2018"};

  string inputfile = Form("../../Utility/FileList/%s/runNumberRunLog.list",globCons::mBeamType[beamType].c_str());
  
  std::cout << "inputfile = " << inputfile.c_str() << std::endl;
  std::ifstream file_runList ( inputfile.c_str() );
  if ( !file_runList.is_open() )
  {
    std::cout << "Abort. Fail to read in run list: " << inputfile << std::endl;
    return -1;
  }

  int temp_runId = 0;
  int temp_runIndex = 0;
  std::cout << "reading correction factors: " << std::endl;
  while (file_runList >> temp_runId)
  {
    std::cout << "temp_runId: " << temp_runId << ", temp_runIndex: " << temp_runIndex << endl;
    runId[temp_runIndex] = temp_runId;
    runIndex[temp_runIndex] = temp_runIndex;
    if(beamType == 1) runIndex[temp_runIndex] = temp_runIndex + 2500;
    std::cout << "runId: " << runId[temp_runIndex] << ", runIndex: " << runIndex[temp_runIndex] << endl;
    temp_runIndex++;
    std::cout << endl;
  }
  file_runList.close();

  string outputfile = Form("../../Utility/RunIndex/%s/runIndex_%s.txt",globCons::mBeamType[beamType].c_str(),globCons::mBeamType[beamType].c_str());
  std::ofstream file_runIndex;
  file_runIndex.open(outputfile.c_str());
  if (!file_runIndex.is_open()) 
  {
    std::cout << "failed to open " << outputfile << '\n';
    return -1;
  } 
  else 
  {
    // write
    for(int i_run = 0; i_run < numOfRuns; ++i_run)
    {
      if(runId[i_run] > 0)
      {
	// std::cout << runId[i_run] << "    " << runIndex[i_run] << std::endl;
	file_runIndex << runId[i_run] << "    " << runIndex[i_run] << std::endl;
      }
    }
  }
  file_runIndex.close();

  return 1;
}
