#!/bin/bash
date

#. ./getFileList.sh

if [ $# -eq 0 ]
then
  OutPutDir="../../../Utility/FileList/RuRu200GeV_2018"
  rm $OutPutDir/pico.list
  rm $OutPutDir/pico_sorted.list
  rm $OutPutDir/runNumberPicoDst.list

  get_file_list.pl -keys 'path,filename' -cond 'production=P20ic,library=SL20c,trgsetupname=production_isobar_2018,collision=RuRu200,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico.list

  awk -F/ '{print $NF, $0}' $OutPutDir/pico.list | sort | uniq | cut -f2- -d ' ' > $OutPutDir/pico_sorted.list
  awk -F/ '{print $(NF-1), $0}' $OutPutDir/pico.list | cut -d ' ' -f 1 | sort | uniq > $OutPutDir/runNumberPicoDst.list

else
  echo "Wrong number of parameters"
fi
