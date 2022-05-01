#!/bin/bash
rm pico.list
rm pico_sorted.list
rm runNumber.list

get_file_list.pl -keys 'path,filename' -cond 'production=P18ic,trgsetupname=AuAu54_production_2017,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico.list

awk -F/ '{print $NF, $0}' pico.list | sort | uniq | cut -f2- -d ' ' > pico_sorted.list
awk -F/ '{print $(NF-1), $0}' pico.list | cut -d ' ' -f 1 | sort | uniq > runNumber.list
