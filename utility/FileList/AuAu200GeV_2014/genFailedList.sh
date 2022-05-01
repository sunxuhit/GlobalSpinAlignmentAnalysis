#!/bin/bash
date

#. ./genFailedList.sh

if [ $# -eq 0 ]
then
  Energy=200GeV_2014
  Luminosity=low
  JobId=9E5703EB6FAE0E39F93C889E8039552F #generate faild list for this Job
  Task=EventPlaneMaker
  Mode=GainCorr
  # Task=RunQA
  # Mode=RunQA

  OutPutDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/${Energy}"
  SubmitDir="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/${Task}/submit/${Energy}_${Luminosity}/JOBS/list"
  echo $SubmitDir

  LogDirectory="/star/u/sunxuhit/AuAu$Energy/Log/${Mode}"

  CompletedList="$OutPutDir/condor_completed_${Task}_${JobId}.list"
  rm $CompletedList
  touch $CompletedList
  cd $LogDirectory
  # grep -l "exiting normally" *${JOBS}*.out | sort > $CompletedList
  grep -l "Work done" *${JobId}*.log | sort > $CompletedList
  cd $OutPutDir
  sed -i 's/^/sched/g' $CompletedList
  sed -i 's/log/list/g' $CompletedList

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
