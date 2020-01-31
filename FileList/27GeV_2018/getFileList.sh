#!/bin/bash

get_file_list.pl -keys 'path,filename' -cond 'production=P19ib,trgsetupname=27GeV_production_2018,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico.list
sort -t '/' -k 13 pico.list | uniq > pico_sorted.list
cut -d '/' -f 12 pico_sorted.list | sort | uniq > runNumber.list
