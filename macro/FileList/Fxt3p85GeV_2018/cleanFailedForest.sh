#!/bin/bash
date

#. ./cleanFailedForest.sh

if [ $# -eq 0 ]
then
  BeamType=Fxt3p85GeV_2018
  JobId=CF3F4859C750A7F521F48D5C095513E3 #generate faild list for this Job
  Task=PhiMesonMaker
  # Mode=RecoPhiSE
  Mode=RecoPhiME

  FileDirectory="/star/u/sunxuhit/$BeamType/SpinAlignment/$Task/Forest"
  OutPutDir="/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/${BeamType}"
  FailedList="$OutPutDir/condor_failed_${Task}_${JobId}.list" # failed condor list

  FailedRootFiles="$OutPutDir/condor_rootFailed_${Task}_${JobId}.list"
  cp $FailedList $FailedRootFiles
  sed -i "s/sched/file_"$Mode"_"$BeamType"_/g" $FailedRootFiles
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
