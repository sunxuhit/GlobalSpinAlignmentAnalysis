#!/bin/bash
date

#. ./removeBadRun.sh

if [ $# -eq 0 ]
then
  OutPutList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/27GeV_2018/pico_xrootd_production.list"
  rm $OutPutList
  touch $OutPutList

  FullPicoList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/27GeV_2018/pico_xrootd_full.list"

  BadRunIdList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/27GeV_2018/badRunStRefMultCorr.txt"

  BadRunTempList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/27GeV_2018/pico_xrootd_badRun_temp.list"
  rm $BadRunTempList
  touch $BadRunTempList
  for item in `cat $BadRunIdList`
  do
    cat $FullPicoList | grep $item >> $BadRunTempList
  done

  BadRunPicoList="/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/27GeV_2018/pico_xrootd_badRun.list"
  rm $BadRunPicoList
  touch $BadRunPicoList
  sort -t '/' -k 16 $BadRunTempList | uniq > $BadRunPicoList
  rm $BadRunTempList

  grep -Fvxf $BadRunPicoList $FullPicoList > $OutPutList
fi

