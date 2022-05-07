#!/bin/bash
rm pico_merged.list
rm pico_merged_sorted.list
rm runNumber_merged.list

echo "merge pico.list"
cat pico_prod.list >> pico_merged.list
cat pico_low.list >> pico_merged.list
cat pico_mid.list >> pico_merged.list
cat pico_high.list >> pico_merged.list

echo "merge pico_sorted.list"
cat pico_prod_sorted.list >> pico_merged_sorted.list
cat pico_low_sorted.list >> pico_merged_sorted.list
cat pico_mid_sorted.list >> pico_merged_sorted.list
cat pico_high_sorted.list >> pico_merged_sorted.list

echo "merge runNumber.list"
cat runNumber_*.list | sort | uniq >> runNumber_merged.list
