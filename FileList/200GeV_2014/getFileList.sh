#!/bin/bash
rm pico.list
rm pico_sorted.list
rm runNumber.list

get_file_list.pl -keys 'path,filename' -cond 'production=P18ih,library=SL20a,trgsetupname=AuAu_200_production_2014,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico.list
get_file_list.pl -keys 'path,filename' -cond 'production=P16id,library=SL18f,trgsetupname=AuAu_200_production_low_2014||AuAu_200_production_mid_2014,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct >> pico.list
get_file_list.pl -keys 'path,filename' -cond 'production=P15ic,library=SL18f,trgsetupname=AuAu_200_production_high_2014,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct >> pico.list

awk -F/ '{print $NF, $0}' pico.list | sort | uniq | cut -f2- -d ' ' > pico_sorted.list
awk -F/ '{print $(NF-1), $0}' pico.list | cut -d ' ' -f 1 | sort | uniq > runNumber.list
