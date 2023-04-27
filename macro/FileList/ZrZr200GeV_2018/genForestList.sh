#!/bin/bash
date

#. ./genForestList.sh

if [ $# -eq 0 ]
then
  BeamType=ZrZr200GeV_2018
  # JobId=820136BE69B6CDF4F20D30A77C669DDB #generate forest list for this Job
  Task=PhiMesonMaker
  Mode=SE
  # Mode=ME

  OutPutDir="/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/${BeamType}"
  cd $OutPutDir

  ForestList="$OutPutDir/forestRecoPhi${Mode}prod_${BeamType}.list"
  rm $ForestList

  FileDir="/star/data01/pwg/sunxuhit/$BeamType/SpinAlignment/$Task/Forest"
  ls -d $FileDir/file_RecoPhi${Mode}_${BeamType}*.root > $ForestList
  # ForestList="$OutPutDir/Forest_RecoPhi${Mode}_${BeamType}_${JobId}.list"
  # ls -d $FileDir/file_RecoPhi${Mode}_${BeamType}_${JobId}*.root > $ForestList

else
  echo "Wrong number of parameters"
fi
