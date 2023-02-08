#!/bin/bash
date

#. ./cleanOutPut.sh

if [ $# -eq 0 ]
then
  BeamType=ZrZr200GeV_2018
  JobId=33CD03E2749A129BF5688C2682C98EB5 #generate faild list for this Job
  Task=RunQA

  FileDirectory="/star/u/sunxuhit/$BeamType/SpinAlignment/$Task/Data"
  OutPutDir="/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/${BeamType}"
  FailedList="$OutPutDir/condor_failed_${Task}_${JobId}.list" # failed condor list

  FailedRootFiles="$OutPutDir/condor_rootFailed_${Task}_${JobId}.list"
  cp $FailedList $FailedRootFiles
  sed -i "s/sched/file_"$BeamType"_"$Task"_/g" $FailedRootFiles
  sed -i "s/list/root/g" $FailedRootFiles

  cd $FileDirectory
  ProducedRootFiles="$OutPutDir/condor_rootProduced_${Task}_${JobId}.list" # all produced ROOT files
  rm $ProducedRootFiles
  touch $ProducedRootFiles
  ls -d *${JobId}*.root | sort > $ProducedRootFiles

  RootFilesToDelete="$OutPutDir/condor_rootToDelete_${Task}_${JobId}.list" # all produced ROOT files
  rm $RootFilesToDelete
  touch $RootFilesToDelete
  comm -12 $FailedRootFiles $ProducedRootFiles | sort > $RootFilesToDelete

  rm $FailedRootFiles
  rm $ProducedRootFiles

  echo "delete following files from "
  echo $RootFilesToDelete
  for item in `cat $RootFilesToDelete`
  do
    echo deleting $FileDirectory/$item
    ls $FileDirectory/$item
    rm $FileDirectory/$item
  done
  rm $RootFilesToDelete

else
  echo "Wrong number of parameters"
fi