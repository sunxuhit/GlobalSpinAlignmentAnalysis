#!/bin/bash
rm pico_low.list
rm pico_low_sorted.list
rm runNumber_low.list

get_file_list.pl -keys 'path,filename' -cond 'production=P16id,library=SL18f,trgsetupname=AuAu_200_production_low_2014,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico_low.list

awk -F/ '{print $NF, $0}' pico_low.list | sort | uniq | cut -f2- -d ' ' > pico_low_sorted.list
awk -F/ '{print $(NF-1), $0}' pico_low.list | cut -d ' ' -f 1 | sort | uniq > runNumber_low.list
