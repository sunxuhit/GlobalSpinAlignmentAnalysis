#!/bin/bash
rm pico_prod.list
rm pico_prod_sorted.list
rm runNumber_prod.list

get_file_list.pl -keys 'path,filename' -cond 'production=P18ih,library=SL20a,trgsetupname=AuAu_200_production_2014,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico_prod.list

awk -F/ '{print $NF, $0}' pico_prod.list | sort | uniq | cut -f2- -d ' ' > pico_prod_sorted.list
awk -F/ '{print $(NF-1), $0}' pico_prod.list | cut -d ' ' -f 1 | sort | uniq > runNumber_prod.list
