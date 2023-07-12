#!/bin/sh

codePath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis

##########Beam Type Selection##########
beamType=2  # Fxt3p85GeV_2018
mode=0 # 0: phi meson QA
flagME=0 # 0: Same Event
strME=SE
# flagME=1 # 1: Mixed Event
# strME=ME
startEvt=0
stopEvt=100000024
library=SL20d
listPath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/Fxt3p85GeV_2018
outPath=/star/data01/pwg/sunxuhit/Fxt3p85GeV_2018/SpinAlignment/PhiMesonAnalyzer
##########Beam Type Selection##########

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

##########Test Production##########
star-submit-template -template testPhiMesonTemp.xml -entities beamType=$beamType,mode=$mode,flagME=$flagME,strME=$strME,startEvt=$startEvt,stopEvt=$stopEvt,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Test Production##########

##########Full Production##########
# star-submit-template -template PhiMesonTemp.xml -entities beamType=$beamType,mode=$mode,flagME=$flagME,strME=$strME,startEvt=$startEvt,stopEvt=$stopEvt,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Full Production##########

##########Re-Submit##########
# star-submit-template -template resubmitPhiMesonTemp.xml -entities beamType=$beamType,mode=$mode,flagME=$flagME,strME=$strME,startEvt=$startEvt,stopEvt=$stopEvt,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Re-Submit##########
