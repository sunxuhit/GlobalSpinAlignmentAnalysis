#!/bin/bash
date

#. ./cleanOutPut.sh

if [ $# -eq 0 ]
then
  Energy=200GeV_2014
  JobId=9E5703EB6FAE0E39F93C889E8039552F #generate faild list for this Job
  Task=EventPlaneMaker
  Mode=GainCorr
  # Task=RunQA
  # Mode=RunQA

  FileDirectory="/star/u/sunxuhit/AuAu$Energy/SpinAlignment/$Mode"

  OutPutDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/$Energy"
  FailedList="$OutPutDir/condor_failed_${Task}_${JobId}.list"
  echo "delete following files from "
  echo $FailedList
  for item in `cat $FailedList`
  do
    echo deleting $FileDirectory/$item
    ls $FileDirectory/$item
    rm $FileDirectory/$item
  done
  rm $FailedList

else
  echo "Wrong number of parameters"
fi
