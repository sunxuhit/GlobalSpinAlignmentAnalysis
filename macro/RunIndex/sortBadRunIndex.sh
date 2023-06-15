#!/bin/bash
date

#. ./sortBadRunIndex.sh

if [ $# -eq 0 ]
then
  # BeamType=ZrZr200GeV_2018
  # BeamType=RuRu200GeV_2018
  BeamType=Fxt3p85GeV_2018
  InPutList="../../Utility/RunIndex/${BeamType}/badRunIndexUnSorted_${BeamType}.txt"
  OutPutList="../../Utility/RunIndex/${BeamType}/badRunIndex_${BeamType}.txt"

  rm $OutPutList
  touch $OutPutList

  cat $InPutList | sort -n | uniq > $OutPutList

else
  echo "Wrong number of parameters"
fi
