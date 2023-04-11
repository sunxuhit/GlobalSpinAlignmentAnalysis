#!/bin/bash
date

#. ./genFailedListLocal.sh

if [ $# -eq 0 ]
then
  BeamType=RuRu200GeV_2018
  JobId=0FD8D0C8A5ACD2E8B1D6DECF6753E365 #generate faild list for this Job
  Task=EventPlaneMaker
  # Mode=EventPlaneMaker

  LogDirectory="/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/${Task}/${BeamType}/${JobId}"

  OutPutDir="/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/${BeamType}"

  cd $LogDirectory
  CompletedLog="$OutPutDir/condor_completedLog_${Task}_${JobId}.list" # get completed list from run log
  rm $CompletedLog
  touch $CompletedLog
  grep -l "Work done" *${JobId}*.log | sort > $CompletedLog
  sed -i "" 's/^/sched/g' $CompletedLog
  sed -i "" 's/log/list/g' $CompletedLog

  CompletedOut="$OutPutDir/condor_completedOut_${Task}_${JobId}.list" # get completed list from stdout
  rm $CompletedOut
  touch $CompletedOut
  grep -l "exiting normally" *${JOBS}*.out | sort > $CompletedOut
  sed -i "" 's/^/sched/g' $CompletedOut
  sed -i "" 's/out/list/g' $CompletedOut

  CompletedList="$OutPutDir/condor_completed_${Task}_${JobId}.list" # common completed list from run log & stdout
  rm $CompletedList
  touch $CompletedList
  comm -12 $CompletedLog $CompletedOut | sort > $CompletedList

  rm $CompletedLog
  rm $CompletedOut

  SubmitDir="/Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/${Task}/${BeamType}/${JobId}/list"
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

else
  echo "Wrong number of parameters"
fi
