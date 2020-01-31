#!/bin/bash
get_file_list.pl -keys 'path,filename' -cond 'production=P18ic,trgsetupname=AuAu54_production_2017,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico.list
sort -t '/' -k 13 pico.list | uniq > pico_sorted.list
cut -d '/' -f 12 pico_sorted.list | sort | uniq > runNumber.list
