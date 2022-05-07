#!/bin/bash
date

#. ./removeBadRun.sh

if [ $# -eq 0 ]
then

  OutPutList="./pico_xrootd_production.list"
  rm $OutPutList
  touch $OutPutList

  FullPicoList="./pico_xrootd_full.list"

  BadRunIdList="./badRunStRefMultCorr.list"
  BadPicoList="./badPico.list"

  BadRunTempList="./pico_xrootd_badRun_temp.list"
  rm $BadRunTempList
  touch $BadRunTempList
  for item in `cat $BadRunIdList` # get picos with bad run Id
  do
    cat $FullPicoList | grep $item >> $BadRunTempList
  done
  cat $BadPicoList >> $BadRunTempList # get bad picos identfied through the QA test

  BadRunPicoList="./pico_xrootd_badRun.list"
  rm $BadRunPicoList
  touch $BadRunPicoList
  sort -t '/' -k 16 $BadRunTempList | uniq > $BadRunPicoList
  rm $BadRunTempList

  grep -Fvxf $BadRunPicoList $FullPicoList > $OutPutList
fi

