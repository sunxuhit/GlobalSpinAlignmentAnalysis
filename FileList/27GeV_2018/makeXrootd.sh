#!/bin/bash

rm pico_xrootd_full.list
sed -e 's#^#root://xrdstar.rcf.bnl.gov:1095/#' pico_sorted.list > pico_xrootd_full.list
