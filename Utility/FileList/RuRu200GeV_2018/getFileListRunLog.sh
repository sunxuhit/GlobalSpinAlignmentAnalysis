#!/bin/bash
rm runNumberRunLogTemp.list
rm runNumberRunLog.list

cat ../Isobar200GeV_2018/allIsobarRun.list | grep "Ru" > runNumberRunLogTemp.list

awk -F' ' '{print $1}' runNumberRunLogTemp.list | sort | uniq > runNumberRunLog.list
