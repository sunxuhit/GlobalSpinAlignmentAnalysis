#!/bin/bash
date

#. ./sortBadRunIndex.sh

if [ $# -eq 0 ]
then
  BeamType=ZrZr200GeV_2018
  # BeamType=RuRu200GeV_2018
  InPutList="./${BeamType}/badRunIndexUnSorted_${BeamType}.txt"
  OutPutList="./${BeamType}/badRunIndex_${BeamType}.txt"

  rm $OutPutList
  touch $OutPutList

  cat $InPutList | sort -n | uniq > $OutPutList

else
  echo "Wrong number of parameters"
fi
