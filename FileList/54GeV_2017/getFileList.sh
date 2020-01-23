#!/bin/bash
get_file_list.pl -keys 'path,filename' -cond 'production=P18ic,trgsetupname=AuAu54_production_2017,filetype=daq_reco_picoDst,filename~st_physics,storage!=hpss' -limit 0 -delim / -distinct > pico_hpss.list
sort -t '/' -k 10 pico_hpss.list > pico_hpss_sorted.list
cut -d '/' -f 10 pico_hpss_sorted.list | sort | uniq >>newrunNumber_54GeV.list
