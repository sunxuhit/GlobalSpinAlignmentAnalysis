#!/bin/bash
date

#. ./removeBadRun.sh

if [ $# -eq 0 ]
then
  OutPutDir="../../../Utility/FileList/Fxt3p85GeV_2018"
  OutPutList="$OutPutDir/pico_xrootd_production.list"
  rm $OutPutList
  touch $OutPutList

  FullPicoList="$OutPutDir/pico_xrootd_full.list"

  BadRunIdList="$OutPutDir/badRunStRefMultCorr.list"
  BadRunQaList="$OutPutDir/badRunStRunQAMaker.list"

  BadRunTempList="$OutPutDir/pico_xrootd_badRun_temp.list"
  rm $BadRunTempList
  touch $BadRunTempList
  for item in `cat $BadRunIdList` # get picos with bad run Id from StRefMultCorr
  do
    cat $FullPicoList | grep $item >> $BadRunTempList
  done
  for item in `cat $BadRunQaList` # get picos with bad run identfied from StRunQAMaker
  do
    cat $FullPicoList | grep $item >> $BadRunTempList
  done
  # cat $BadRunQaList >> $BadRunTempList # get picos with bad run identfied through the QA test

  BadRunPicoList="$OutPutDir/pico_xrootd_badRun.list"
  rm $BadRunPicoList
  touch $BadRunPicoList
  sort -t '/' -k 16 $BadRunTempList | uniq > $BadRunPicoList
  rm $BadRunTempList

  grep -Fvxf $BadRunPicoList $FullPicoList > $OutPutList

else
  echo "Wrong number of parameters"
fi

