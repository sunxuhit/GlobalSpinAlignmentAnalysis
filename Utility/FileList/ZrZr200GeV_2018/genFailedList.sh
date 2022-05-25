#!/bin/bash
date

#. ./genFailedList.sh

if [ $# -eq 0 ]
then
  BeamType=ZrZr200GeV_2018
  JobId=33CD03E2749A129BF5688C2682C98EB5 #generate faild list for this Job
  Task=RunQA
  Mode=RunQA

  LogDirectory="/star/u/sunxuhit/$BeamType/SpinAlignment/${Mode}/Log"

  OutPutDir="/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/${BeamType}"
  cd $OutPutDir

  cd $LogDirectory
  CompletedLog="$OutPutDir/condor_completedLog_${Task}_${JobId}.list" # get completed list from run log
  rm $CompletedLog
  touch $CompletedLog
  grep -l "Work done" *${JobId}*.log | sort > $CompletedLog
  sed -i 's/^/sched/g' $CompletedLog
  sed -i 's/log/list/g' $CompletedLog

  CompletedOut="$OutPutDir/condor_completedOut_${Task}_${JobId}.list" # get completed list from stdout
  rm $CompletedOut
  touch $CompletedOut
  grep -l "exiting normally" *${JOBS}*.out | sort > $CompletedOut
  sed -i 's/^/sched/g' $CompletedOut
  sed -i 's/out/list/g' $CompletedOut

  CompletedList="$OutPutDir/condor_completed_${Task}_${JobId}.list" # common completed list from run log & stdout
  rm $CompletedList
  touch $CompletedList
  comm -12 $CompletedLog $CompletedOut | sort > $CompletedList

  rm $CompletedLog
  rm $CompletedOut

  SubmitDir="/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/${Task}/${BeamType}/JOBS/list"
  echo $SubmitDir
  SubmittedList="$OutPutDir/condor_submitted_${Task}_${JobId}.list"
  rm $SubmittedList
  touch $SubmittedList
  cd $SubmitDir
  ls -d *${JobId}*.list | sort > $SubmittedList
  cd $OutPutDir

  FailedList="$OutPutDir/condor_failed_${Task}_${JobId}.list"
  rm $FailedList
  touch $FailedList
  comm -13 $CompletedList $SubmittedList > $FailedList

  rm $SubmittedList
  rm $CompletedList

  ResubmitList="$OutPutDir/pico_xrootd_resubmit_${Task}_${JobId}.list"
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

  #generate failed ROOT output list
  # sed -i "s/sched/file_"$BeamType"_"$Mode"_/g" $FailedList
  # sed -i "s/list/root/g" $FailedList

else
  echo "Wrong number of parameters"
fi
