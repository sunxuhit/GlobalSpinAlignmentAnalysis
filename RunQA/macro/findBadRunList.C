#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLegend.h"

#include "../StRoot/StRunQAMaker/StRunQACons.h"

using namespace std;

std::map<int,int> map_runIndex; // runIndex, runId
bool read_in_runIndex(int energy);
int findRunId(int runIndex);

int findBadRunList(int energy = 0)
{
  map_runIndex.clear();
  bool isOpen_runIndex = read_in_runIndex(energy);
  if(isOpen_runIndex) std::cout << "Run Index read in!" << std::endl;

  // read in bad runIndex list
  std::string inputfile = Form("../StRoot/StRunQAUtility/RunIndex/badRunIndex_%s.txt",runQA::mBeamEnergy[energy].c_str());
  std::cout << "inputfile = " << inputfile.c_str() << std::endl;
  std::ifstream file_badRunIndex ( inputfile.c_str() );
  if ( !file_badRunIndex.is_open() )
  {
    std::cout << "Abort. Fail to read in bad Run Index file: " << inputfile << std::endl;
    return -1;
  }

  int badRunIndex = 0;
  std::vector<int> badRunId;
  badRunId.clear();
  std::cout << "reading bad Bun Index: " << std::endl;
  while (file_badRunIndex >> badRunIndex)
  {
    std::cout << "badRunIndex = " << badRunIndex << std::endl;
    int runId = findRunId(badRunIndex);
    badRunId.push_back(runId);
  }
  file_badRunIndex.close();

  // generate bad runId list
  std::string outputfile = Form("../StRoot/StRunQAUtility/RunIndex/badRunList_%s.txt",runQA::mBeamEnergy[energy].c_str());
  std::ofstream file_badRunList;
  file_badRunList.open(outputfile.c_str());
  if (!file_badRunList.is_open()) 
  {
    std::cout << "failed to open " << outputfile.c_str() << '\n';
    return -1;
  } 
  else 
  {
    // write
    for(int i_runId = 0; i_runId < badRunId.size(); ++i_runId)
    {
      cout << "i_runId = " << i_runId << ", badRunId = " << badRunId[i_runId] << endl;
      file_badRunList << badRunId[i_runId] << std::endl;
    }
  }
  file_badRunList.close();

  return 1;
}

bool read_in_runIndex(int energy = 0)
{
  // read in runId vs. runIndex
  std::string inputfile = Form("../StRoot/StRunQAUtility/RunIndex/runIndex_%s.txt",runQA::mBeamEnergy[energy].c_str());
  std::cout << "inputfile = " << inputfile.c_str() << std::endl;
  std::ifstream file_runIndex ( inputfile.c_str() );
  if ( !file_runIndex.is_open() )
  {
    std::cout << "Abort. Fail to read in run Index file: " << inputfile << std::endl;
    return -1;
  }

  int temp_runId = 0, temp_runIndex = 0;
  std::cout << "reading run Index: " << std::endl;
  while (file_runIndex >> temp_runId >> temp_runIndex)
  {
    // std::cout << "runId = " << temp_runId << ", runIndex = " << temp_runIndex << std::endl;
    map_runIndex[temp_runIndex] = temp_runId;
  }
  file_runIndex.close();
  std::cout << "Run Index read in!" << std::endl;

  return true;
}

int findRunId(int runIndex)
{
  // print map_runIndex content:
  // for (std::map<int,int>::iterator it=map_runIndex.begin(); it!=map_runIndex.end(); ++it)
  // {
  //   std::cout << it->first << " => " << it->second << '\n';
  // }

  std::map<int,int>::iterator it_runIndex = map_runIndex.find(runIndex);
  if(it_runIndex == map_runIndex.end())
  {
    // std::cout << "StRunQAUtility -> could not find in full run list! & send signal to kill the run!" << std::endl;
    return -999;
  }
  else
  {
    std::cout << "runIndex: " << it_runIndex->first << " => runId: " << it_runIndex->second << std::endl;
    return it_runIndex->second;
  }

  return -999;
}
