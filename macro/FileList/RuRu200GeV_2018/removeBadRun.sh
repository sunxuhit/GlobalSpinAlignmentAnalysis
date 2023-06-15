#!/bin/bash
date

#. ./removeBadRun.sh

if [ $# -eq 0 ]
then
  OutPutDir="../../../Utility/FileList/RuRu200GeV_2018"
  OutPutList="$OutPutDir/pico_xrootd_production.list"
  rm $OutPutList
  touch $OutPutList

  FullPicoList="$OutPutDir/pico_xrootd_full.list"

  BadRunList="$OutPutDir/runNumberBadRun.list"

  BadRunTempList="$OutPutDir/pico_xrootd_badRun_temp.list"
  rm $BadRunTempList
  touch $BadRunTempList
  for item in `cat $BadRunList` # get picos with bad run identfied from both StRefMultCorr & StRunQAMaker
  do
    cat $FullPicoList | grep $item >> $BadRunTempList
  done

  BadRunPicoList="$OutPutDir/pico_xrootd_badRun.list"
  rm $BadRunPicoList
  touch $BadRunPicoList
  awk -F/ '{print $NF, $0}' $BadRunTempList | sort | uniq | cut -f2- -d ' ' > $BadRunPicoList
  rm $BadRunTempList

  grep -Fvxf $BadRunPicoList $FullPicoList > $OutPutList

else
  echo "Wrong number of parameters"
fi

