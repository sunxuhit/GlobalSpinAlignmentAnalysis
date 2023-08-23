#!/bin/bash
date

#. ./syncEventPlaneUtil.sh

if [ $# -eq 0 ]
then
  # echo "synyc Correction files of ZrZr200GeV_2018"
  # rsync -avz /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/EventPlaneMaker/ZrZr200GeV_2018/ rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/EventPlaneMaker/ZrZr200GeV_2018/

  # echo "synyc Correction files of RuRu200GeV_2018"
  # rsync -avz /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/EventPlaneMaker/RuRu200GeV_2018/ rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/EventPlaneMaker/RuRu200GeV_2018/

  echo "synyc Correction files of Fxt3p85GeV_2018"
  rsync -avz /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/EventPlaneMaker/Fxt3p85GeV_2018/ rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/EventPlaneMaker/Fxt3p85GeV_2018/

  # echo "synyc Correction files"
  # rsync -avz /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/EventPlaneMaker/ rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/Utility/EventPlaneMaker/

else
  echo "Wrong number of parameters"
fi
