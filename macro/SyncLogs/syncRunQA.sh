#!/bin/bash
date

#. ./syncRunQA.sh

if [ $# -eq 0 ]
then
  echo "synyc logs & Jobs of ZrZr200GeV_2018"
  rsync -avz rcf:/star/u/sunxuhit/ZrZr200GeV_2018/SpinAlignment/RunQA/Log/archivedFiles/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/RunQA/ZrZr200GeV_2018
  rsync -avz rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/RunQA/ZrZr200GeV_2018/archivedJobs/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/RunQA/ZrZr200GeV_2018

  echo "synyc logs & Jobs of RuRu200GeV_2018"
  rsync -avz rcf:/star/u/sunxuhit/RuRu200GeV_2018/SpinAlignment/RunQA/Log/archivedFiles/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/RunQA/RuRu200GeV_2018
  rsync -avz rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/RunQA/RuRu200GeV_2018/archivedJobs/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/RunQA/RuRu200GeV_2018

  echo "synyc logs & Jobs of Fxt3p85GeV_2018"
  rsync -avz rcf:/star/u/sunxuhit/Fxt3p85GeV_2018/SpinAlignment/RunQA/Log/archivedFiles/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/log/RunQA/Fxt3p85GeV_2018
  rsync -avz rcf:/star/u/sunxuhit/WorkSpace/SpinAlignment/GlobalSpinAlignmentAnalysis/submit/RunQA/Fxt3p85GeV_2018/archivedJobs/ /Users/xusun/WorkSpace/STAR/SpinAlignment/GlobalSpinAlignmentAnalysis/jobs/RunQA/Fxt3p85GeV_2018

else
  echo "Wrong number of parameters"
fi
