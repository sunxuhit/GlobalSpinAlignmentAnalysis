#!/bin/bash
date

#. ./cleanOutPut.sh

if [ $# -eq 0 ]
then
  Energy=27GeV_2018
  SM=SE
  Mode=QA

  FileDirectory="/star/u/sunxuhit/AuAu$Energy/SpinAlignment/$Mode"

  OutPutDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/$Energy"
  FailedList="$OutPutDir/condor_failed.list"
  echo "delete following files from "
  echo $FailedList
  for item in `cat $FailedList`
  do
    ls $FileDirectory/$item
    echo deleting $FileDirectory/$item
    # rm $item
  done
  # rm $FailedList

else
  echo "Wrong number of parameters"
fi
