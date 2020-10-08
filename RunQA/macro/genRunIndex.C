#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <iterator>
#include "../StRoot/StRunQAMaker/StRunQACons.h"

int genRunIndex(int energy = 0)
{
  const int numOfRuns = 4000;
  int runId[numOfRuns];
  int runIndex[numOfRuns];
  for(int i_run = 0; i_run < numOfRuns; ++i_run)
  {
    runId[i_run] = -999;
    runIndex[i_run] = -999;
  }

  string inputfile = Form("/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/%s/runNumber.list",runQA::mBeamEnergy[energy].c_str());
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
    std::cout << "runId: " << runId[temp_runIndex] << ", runIndex: " << runIndex[temp_runIndex] << endl;
    temp_runIndex++;
    std::cout << endl;
  }
  file_runList.close();

  std::string outputfile = Form("/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/RunQA/StRoot/StRunQAUtility/RunIndex/runIndex_%s.txt",runQA::mBeamEnergy[energy].c_str());
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
