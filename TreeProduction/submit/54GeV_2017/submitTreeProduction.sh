#!/bin/sh

codePath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/TreeProduction

##########Energy Selection##########
# energy=0  # 200GeV
# library=SL18h
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/200GeV_2014
# outPath=/star/data01/pwg/sunxuhit/AuAu200GeV_2014
 
energy=1  # 54.0GeV
library=SL18c
listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/54GeV_2017
outPath=/star/data01/pwg/sunxuhit/AuAu54GeV_2017

# energy=2  # 27GeV
# library=SL19b
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/27GeV_2018
# outPath=/star/data01/pwg/sunxuhit/AuAu27GeV_2018
##########Energy Selection##########

##########Mode Selection##########
mode=0
outDir=QA

# mode=1
# outDir=ReCenterParameter

# mode=2
# outDir=ShiftParameter

# mode=3
# outDir=Resolution

# mode=4
# outDir=Phi/Forest
##########Mode Selection##########

##########Mixed Event Selection##########
flag_ME=0 # 0 for SE | 1 for ME
SM=SE
# flag_ME=1 # 0
# SM=ME
##########Mixed Event Selection##########


mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

##########Test Production##########
star-submit-template -template testProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template TreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
# star-submit-template -template resubmitTreeProductionTemp.xml -entities mode=$mode,energy=$energy,flag_ME=$flag_ME,SM=$SM,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
