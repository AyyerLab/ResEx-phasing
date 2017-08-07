#!/bin/sh

fsc_file='FSC_run_sequence.txt'
fsc_value='0.400'
awk -v val="$fsc_value" '$1 == val' fsc-0-0.dat

