#!/bin/sh

codePath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis

##########Beam Type Selection##########
mode=0 # 0: phi meson TTree
beamType=0  # ZrZr200GeV_2018
flagME=0 # 0: Same Event | 1: Mixed Event
library=SL20c
listPath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/ZrZr200GeV_2018
outPath=/star/data01/pwg/sunxuhit/ZrZr200GeV_2018/SpinAlignment/PhiMesonMaker
##########Beam Type Selection##########

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

##########Test Production##########
star-submit-template -template testPhiMesonTemp.xml -entities mode=$mode,beamType=$beamType,flagME=$flagME,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Test Production##########

##########Full Production##########
# star-submit-template -template PhiMesonTemp.xml -entities mode=$mode,beamType=$beamType,flagME=$flagME,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Full Production##########

##########Re-Submit##########
# star-submit-template -template resubmitPhiMesonTemp.xml -entities mode=$mode,beamType=$beamType,flagME=$flagME,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Re-Submit##########
