#!/bin/sh

codePath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis

##########Beam Type Selection##########
beamType=1  # RuRu200GeV_2018
library=SL20c
listPath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/RuRu200GeV_2018
outPath=/star/data01/pwg/sunxuhit/RuRu200GeV_2018/SpinAlignment/RunQA
##########Beam Type Selection##########

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

##########Test Production##########
star-submit-template -template testRunQATemp.xml -entities beamType=$beamType,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Test Production##########

##########Full Production##########
# star-submit-template -template RunQATemp.xml -entities beamType=$beamType,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Full Production##########

##########Re-Submit##########
# star-submit-template -template resubmitRunQATemp.xml -entities beamType=$beamType,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Re-Submit##########
