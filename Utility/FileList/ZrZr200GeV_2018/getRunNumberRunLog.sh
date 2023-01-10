#!/bin/bash
rm runNumberRunLog.list

cat ../Isobar200GeV_2018/allIsobarRun.list | grep "Zr" > runNumberRunLogTemp.list

awk -F' ' '{print $1}' runNumberRunLogTemp.list | sort | uniq > runNumberRunLog.list

rm runNumberRunLogTemp.list
