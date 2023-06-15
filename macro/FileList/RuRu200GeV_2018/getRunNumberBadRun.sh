#!/bin/bash
date

#. ./.getRunNumberBadRun.sh

if [ $# -eq 0 ]
then
  OutPutDir="../../../Utility/FileList/RuRu200GeV_2018"
  OutPutList="$OutPutDir/badRunRuRu200GeV_2018.list"
  rm $OutPutList
  touch $OutPutList

  BadRunIdList="$OutPutDir/badRunStRefMultCorr.list"
  BadRunQaList="$OutPutDir/badRunStRunQAMaker.list"

  BadRunTempList="$OutPutDir/badRunRuRu200GeV_2018_temp.list"
  rm $BadRunTempList
  touch $BadRunTempList
  cat $BadRunIdList >> $BadRunTempList # get bad run Id from StRefMultCorr
  cat $BadRunQaList >> $BadRunTempList # get bad run Id from StRunQAMaker

  awk -F' ' '{print $1}' $BadRunTempList | sort | uniq > $OutPutList

  rm $BadRunTempList

else
  echo "Wrong number of parameters"
fi


