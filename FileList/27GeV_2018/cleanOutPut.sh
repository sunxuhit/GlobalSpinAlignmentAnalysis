#!/bin/bash
date

#. ./cleanOutPut.sh

if [ $# -eq 0 ]
then
  Energy=27GeV_2018
  SM=SE
  Mode=QA
  JobId=5EF1F4FCA6A3BD1C7F92CD9DDF087847 #generate faild list for this Job

  FileDirectory="/star/u/sunxuhit/AuAu$Energy/SpinAlignment/$Mode"

  OutPutDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/$Energy"
  FailedList="$OutPutDir/condor_failed_${Mode}_${SM}_${JobId}.list"
  echo "delete following files from "
  echo $FailedList
  for item in `cat $FailedList`
  do
    echo deleting $FileDirectory/$item
    ls $FileDirectory/$item
    rm $FileDirectory/$item
  done
  # rm $FailedList

else
  echo "Wrong number of parameters"
fi
