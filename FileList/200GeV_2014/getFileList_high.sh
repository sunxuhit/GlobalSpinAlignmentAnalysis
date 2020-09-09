#!/bin/bash
rm pico_high.list
rm pico_high_sorted.list
rm runNumber_high.list

get_file_list.pl -keys 'path,filename' -cond 'production=P15ic,library=SL18f,trgsetupname=AuAu_200_production_high_2014,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico_high.list

awk -F/ '{print $NF, $0}' pico_high.list | sort | uniq | cut -f2- -d ' ' > pico_high_sorted.list
awk -F/ '{print $(NF-1), $0}' pico_high.list | cut -d ' ' -f 1 | sort | uniq > runNumber_high.list
