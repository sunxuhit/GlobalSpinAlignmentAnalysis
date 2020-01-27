#!/bin/bash
date

#. ./VecMesonTree.sh

if [ $# -eq 0 ]
then
  Energy=27GeV_2018
  SM=SE
  Mode=QA

  OutPutDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/$Energy"
  SubmitDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/TreeProduction/submit/$Energy/JOBS/list"

  LogDirectory="/star/u/sunxuhit/AuAu$Energy/Log"

  CompletedList="$OutPutDir/condor_completed.list"
  rm $CompletedList
  touch $CompletedList
  cd $LogDirectory
  grep -l "exiting normally" *${SM}*.out | sort > $CompletedList
  cd $OutPutDir
  sed -i 's/QA_SE_/sched/g' $CompletedList
  sed -i 's/out/list/g' $CompletedList

  SubmittedList="$OutPutDir/condor_submitted.list"
  rm $SubmittedList
  touch $SubmittedList
  cd $SubmitDir
  ls -d sched*.list | sort > $SubmittedList
  cd $OutPutDir

  FailedList="$OutPutDir/condor_failed.list"
  rm $FailedList
  touch $FailedList
  comm -13 $CompletedList $SubmittedList > $FailedList
  rm $SubmittedList
  rm $CompletedList

  ResubmitList="$OutPutDir/pico_xrootd_resubmit.list"
  rm $ResubmitList
  touch $ResubmitList
  TempList="$OutPutDir/pico_xrootd_temp.list"
  touch $TempList
  for item in `cat $FailedList`
  do
    cat $SubmitDir/$item >> $TempList
  done
  sort $TempList | uniq > $ResubmitList
  rm $TempList

  # generate failed ROOT output list
  sed -i "s/sched/file_"$Energy"_"$Mode"_/g" $FailedList
  sed -i "s/list/root/g" $FailedList


  # TempList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/$Energy/Temp_"$Energy".list"
  # cut -d '/' -f 7 $InPutList | sort | uniq > $TempList

else
  echo "Wrong number of parameters"
fi
