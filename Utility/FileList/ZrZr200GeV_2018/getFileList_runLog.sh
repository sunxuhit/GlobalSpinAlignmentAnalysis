#!/bin/bash
rm runNumberFromRunLogTemp.list
rm runNumberFromRunLog.list

cat ../Isobar200GeV_2018/allIsobarRun.list | grep "Zr" > runNumberFromRunLogTemp.list

awk -F' ' '{print $1}' runNumberFromRunLogTemp.list | sort | uniq > runNumberFromRunLog.list
