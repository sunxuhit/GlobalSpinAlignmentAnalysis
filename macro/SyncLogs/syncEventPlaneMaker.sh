#!/bin/bash
date

#. ./syncEventPlaneMaker.sh

if [ $# -eq 0 ]
then
  echo "synyc logs & Jobs of ZrZr200GeV_2018"
  rsync -avz rcf:/star/u/sunxuhit/ZrZr200GeV_2018/SpinAlignment/EventPlaneMaker/Log/archivedFiles/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/EventPlaneMaker/ZrZr200GeV_2018
  rsync -avz rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/EventPlaneMaker/ZrZr200GeV_2018/archivedJobs/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/EventPlaneMaker/ZrZr200GeV_2018

  echo "synyc logs & Jobs of RuRu200GeV_2018"
  rsync -avz rcf:/star/u/sunxuhit/RuRu200GeV_2018/SpinAlignment/EventPlaneMaker/Log/archivedFiles/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/EventPlaneMaker/RuRu200GeV_2018
  rsync -avz rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/EventPlaneMaker/RuRu200GeV_2018/archivedJobs/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/EventPlaneMaker/RuRu200GeV_2018

  echo "synyc logs & Jobs of Fxt3p85GeV_2018"
  rsync -avz rcf:/star/u/sunxuhit/Fxt3p85GeV_2018/SpinAlignment/EventPlaneMaker/Log/archivedFiles/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/EventPlaneMaker/Fxt3p85GeV_2018
  rsync -avz rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/EventPlaneMaker/Fxt3p85GeV_2018/archivedJobs/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/EventPlaneMaker/Fxt3p85GeV_2018

else
  echo "Wrong number of parameters"
fi
