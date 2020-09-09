#!/bin/bash
rm pico_mid.list
rm pico_mid_sorted.list
rm runNumber_mid.list

get_file_list.pl -keys 'path,filename' -cond 'production=P16id,library=SL18f,trgsetupname=AuAu_200_production_mid_2014,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico_mid.list

awk -F/ '{print $NF, $0}' pico_mid.list | sort | uniq | cut -f2- -d ' ' > pico_mid_sorted.list
awk -F/ '{print $(NF-1), $0}' pico_mid.list | cut -d ' ' -f 1 | sort | uniq > runNumber_mid.list
