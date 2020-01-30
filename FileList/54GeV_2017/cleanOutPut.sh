#!/bin/bash
date

#. ./cleanOutPut.sh

if [ $# -eq 0 ]
then
  Energy=54GeV_2017
  SM=SE
  Mode=QA
  JobId=09DDC6912DF49AF532AB39B70796C8FA #generate faild list for this Job

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
