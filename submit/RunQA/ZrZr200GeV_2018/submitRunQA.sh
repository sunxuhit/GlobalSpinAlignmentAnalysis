#!/bin/sh

codePath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis

##########Beam Type Selection##########
beamType=0  # ZrZr200GeV_2018
library=SL20c
listPath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList
outPath=/star/data01/pwg/sunxuhit/ZrZr200GeV_2018/SpinAlignment
##########Beam Type Selection##########

##########Mode Selection##########
mode=0
outDir=RunQA
##########Mode Selection##########

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

##########Test Production##########
star-submit-template -template testRunQATemp.xml -entities mode=$mode,beamType=$beamType,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Test Production##########

##########Full Production##########
# star-submit-template -template RunQATemp.xml -entities mode=$mode,beamType=$beamType,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Full Production##########

##########Re-Submit##########
# star-submit-template -template resubmitRunQATemp.xml -entities mode=$mode,beamType=$beamType,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath,outDir=$outDir
##########Re-Submit##########
