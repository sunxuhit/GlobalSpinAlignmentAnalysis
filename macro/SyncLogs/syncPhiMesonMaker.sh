#!/bin/bash
date

#. ./syncPhiMesonMaker.sh

if [ $# -eq 0 ]
then
  echo "synyc logs & Jobs of ZrZr200GeV_2018"
  rsync -avz rcf:/star/u/sunxuhit/ZrZr200GeV_2018/SpinAlignment/PhiMesonMaker/Log/archivedFiles/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/PhiMesonMaker/ZrZr200GeV_2018
  rsync -avz rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/PhiMesonMaker/ZrZr200GeV_2018/archivedJobs/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/PhiMesonMaker/ZrZr200GeV_2018

  echo "synyc logs & Jobs of RuRu200GeV_2018"
  rsync -avz rcf:/star/u/sunxuhit/RuRu200GeV_2018/SpinAlignment/PhiMesonMaker/Log/archivedFiles/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/PhiMesonMaker/RuRu200GeV_2018
  rsync -avz rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/PhiMesonMaker/RuRu200GeV_2018/archivedJobs/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/PhiMesonMaker/RuRu200GeV_2018

  echo "synyc logs & Jobs of Fxt3p85GeV_2018"
  rsync -avz rcf:/star/u/sunxuhit/Fxt3p85GeV_2018/SpinAlignment/PhiMesonMaker/Log/archivedFiles/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/PhiMesonMaker/Fxt3p85GeV_2018
  rsync -avz rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/PhiMesonMaker/Fxt3p85GeV_2018/archivedJobs/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/PhiMesonMaker/Fxt3p85GeV_2018

else
  echo "Wrong number of parameters"
fi
