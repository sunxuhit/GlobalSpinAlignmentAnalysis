#!/bin/bash
date

#. ./removeBadRun.sh

if [ $# -eq 0 ]
then
  Energy=200GeV_2014

  OutPutList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/${Energy}/pico_xrootd_production.list"
  rm $OutPutList
  touch $OutPutList

  FullPicoList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/${Energy}/pico_xrootd_full.list"

  BadRunIdList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/${Energy}/badRunStRefMultCorr.txt"
  BadPicoList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/${Energy}/badPico.txt"

  BadRunTempList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/${Energy}/pico_xrootd_badRun_temp.list"
  rm $BadRunTempList
  touch $BadRunTempList
  for item in `cat $BadRunIdList` # get picos with bad run Id
  do
    cat $FullPicoList | grep $item >> $BadRunTempList
  done
  cat $BadPicoList >> $BadRunTempList # get bad picos identfied through the QA test

  BadRunPicoList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/${Energy}/pico_xrootd_badRun.list"
  rm $BadRunPicoList
  touch $BadRunPicoList
  sort -t '/' -k 16 $BadRunTempList | uniq > $BadRunPicoList
  rm $BadRunTempList

  grep -Fvxf $BadRunPicoList $FullPicoList > $OutPutList
fi

