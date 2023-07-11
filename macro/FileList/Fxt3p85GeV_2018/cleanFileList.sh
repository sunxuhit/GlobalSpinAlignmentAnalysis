#!/bin/bash
date

#. ./cleanFileList.sh

if [ $# -eq 0 ]
then
  BeamType=Fxt3p85GeV_2018
  # Task=RunQA
  # Task=EventPlaneMaker
  # Task=PhiMesonMaker
  Task=PhiMesonAnalyzer

  # OutPutDir="/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/${BeamType}"
  OutPutDir="../../../Utility/FileList/${BeamType}"

  echo "delete following files:"
  ls -d $OutPutDir/condor_failed_${Task}_*.list
  ls -d $OutPutDir/pico_xrootd_resubmit_${Task}_*.list

  rm $OutPutDir/condor_failed_${Task}_*.list
  rm $OutPutDir/pico_xrootd_resubmit_${Task}_*.list

else
  echo "Wrong number of parameters"
fi
