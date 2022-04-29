#!/bin/sh

codePath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/RunQA

##########Energy Selection##########
energy=0  # 200GeV
library=SL20a
listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/200GeV_2014
outPath=/star/data01/pwg/sunxuhit/AuAu200GeV_2014
 
# energy=1  # 54.0GeV
# library=SL18c
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/54GeV_2017
# outPath=/star/data01/pwg/sunxuhit/AuAu54GeV_2017

# energy=2  # 27GeV
# library=SL19b
# listPath=/star/u/sunxuhit/WorkSpace/VecMesonSpinAlignment_BESII/FileList/27GeV_2018
# outPath=/star/data01/pwg/sunxuhit/AuAu27GeV_2018
##########Energy Selection##########

##########Mode Selection##########
mode=0
outDir=RunQA
##########Mode Selection##########

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

##########Test Production##########
# star-submit-template -template testRunQA_prod.xml -entities mode=$mode,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
star-submit-template -template RunQA_prod.xml -entities mode=$mode,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
# star-submit-template -template resubmitRunQATemp.xml -entities mode=$mode,energy=$energy,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
