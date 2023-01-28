#!/bin/bash
date

#. ./makeXrootd.sh

if [ $# -eq 0 ]
then
  OutPutDir="../../../Utility/FileList/ZrZr200GeV_2018"
  rm $OutPutDir/pico_xrootd_full.list
  sed -e 's#^#root://xrdstar.rcf.bnl.gov:1095/#' $OutPutDir/pico_sorted.list > $OutPutDir/pico_xrootd_full.list

else
  echo "Wrong number of parameters"
fi
