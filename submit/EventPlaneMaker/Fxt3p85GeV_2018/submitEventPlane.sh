#!/bin/sh

codePath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis

##########Beam Type Selection##########
mode=5 # 0: gain & phi wgt | 1: sub EP recenter | 2: sub EP shift | 4: resolution | 5: flow
beamType=2  # Fxt3p85GeV_2018
library=SL20d
listPath=/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/FileList/Fxt3p85GeV_2018
outPath=/star/data01/pwg/sunxuhit/Fxt3p85GeV_2018/SpinAlignment/EventPlaneMaker
##########Beam Type Selection##########

mkdir -p JOBS/report
mkdir -p JOBS/csh
mkdir -p JOBS/list

##########Test Production##########
star-submit-template -template testEventPlaneTemp.xml -entities mode=$mode,beamType=$beamType,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Test Production##########

##########Full Production##########
# star-submit-template -template EventPlaneTemp.xml -entities mode=$mode,beamType=$beamType,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Full Production##########

##########Re-Submit##########
# star-submit-template -template resubmitEventPlaneTemp.xml -entities mode=$mode,beamType=$beamType,library=$library,codePath=$codePath,outPath=$outPath,listPath=$listPath
##########Re-Submit##########
