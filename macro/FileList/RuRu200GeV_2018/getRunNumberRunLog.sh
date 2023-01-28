#!/bin/bash
date

#. ./getRunNumberRunLog.sh

if [ $# -eq 0 ]
then
  OutPutDir="../../../Utility/FileList/RuRu200GeV_2018"
  rm $OutPutDir/runNumberRunLog.list

  cat ../../../Utility/FileList/Isobar200GeV_2018/allIsobarRunNumber.list | grep "Ru" > $OutPutDir/runNumberRunLogTemp.list

  awk -F' ' '{print $1}' $OutPutDir/runNumberRunLogTemp.list | sort | uniq > $OutPutDir/runNumberRunLog.list

  rm $OutPutDir/runNumberRunLogTemp.list

else
  echo "Wrong number of parameters"
fi
