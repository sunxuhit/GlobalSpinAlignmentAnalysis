#!/bin/bash
date

#. ./VecMesonTree.sh

if [ $# -eq 0 ]
then
  Energy=54GeV_2017
  SM=SE
  Mode=QA
  JobId=8E7207D11F7C30A914199841DFBE6A97 #generate faild list for this Job

  OutPutDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/$Energy"
  SubmitDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/TreeProduction/submit/$Energy/JOBS/list"

  LogDirectory="/star/u/sunxuhit/AuAu$Energy/Log"

  CompletedList="$OutPutDir/condor_completed_${Mode}_${SM}_${JobId}.list"
  rm $CompletedList
  touch $CompletedList
  cd $LogDirectory
  # grep -l "exiting normally" *${JOBS}*.out | sort > $CompletedList
  grep -l "Work done" *${JobId}*.out | sort > $CompletedList
  cd $OutPutDir
  sed -i 's/QA_SE_/sched/g' $CompletedList
  sed -i 's/out/list/g' $CompletedList

  SubmittedList="$OutPutDir/condor_submitted_${Mode}_${SM}_${JobId}.list"
  rm $SubmittedList
  touch $SubmittedList
  cd $SubmitDir
  ls -d *${JobId}*.list | sort > $SubmittedList
  cd $OutPutDir

  FailedList="$OutPutDir/condor_failed_${Mode}_${SM}_${JobId}.list"
  rm $FailedList
  touch $FailedList
  comm -13 $CompletedList $SubmittedList > $FailedList
  # rm $SubmittedList
  # rm $CompletedList

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

else
  echo "Wrong number of parameters"
fi
