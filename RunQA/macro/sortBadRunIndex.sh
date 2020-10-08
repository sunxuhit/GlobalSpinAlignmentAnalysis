#!/bin/bash
date

#. ./sortBadRunIndex.sh

if [ $# -eq 0 ]
then
  Energy=200GeV_2014
  InPutList="../StRoot/StRunQAUtility/RunIndex/badRunIndexUnSorted_${Energy}.txt"
  OutPutList="../StRoot/StRunQAUtility/RunIndex/badRunIndex_${Energy}.txt"

  rm $OutPutList
  touch $OutPutList

  cat $InPutList | sort -n | uniq > $OutPutList

else
  echo "Wrong number of parameters"
fi
