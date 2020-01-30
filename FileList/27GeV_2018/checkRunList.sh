#!/bin/bash
date

#. ./VecMesonTree.sh

if [ $# -eq 0 ]
then
  Energy=27GeV_2018
  SM=SE
  Mode=QA
  # JobId=86378FED7FD230281B0DBF5BFDC55F55 #generate faild list for this Job

  OutPutDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/$Energy"
  SubmitDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/TreeProduction/submit/$Energy/JOBS/list"

  CompletedList="$OutPutDir/pico_xrootd_completed_${Mode}_${SM}.list"
  rm $CompletedList
  touch $CompletedList

  TempList="$OutPutDir/pico_xrootd_temp.list"
  rm $TempList
  touch $TempList
  for item in `cat $OutPutDir/condor_completed_${Mode}_${SM}_*.list`
  do
    cat $SubmitDir/$item >> $TempList
  done
  # sort -t '/' -k 16 $TempList > $CompletedList
  sort $TempList > $CompletedList
  rm $TempList
  sed -i "s/[[:blank:]]*$//g" $CompletedList # remove space at the end of each line

  touch $TempList
  # sort -t '/' -k 16 pico_xrootd_production.list > $TempList
  sort pico_xrootd_production.list > $TempList

  FailedList="$OutPutDir/pico_xrootd_failed_${Mode}_${SM}.list"
  rm $FailedList
  touch $FailedList
  comm -13 $CompletedList $TempList > $FailedList
  rm $TempList

else
  echo "Wrong number of parameters"
fi
