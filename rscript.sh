#!/bin/bash

for i in {119..120}
#for i in 116
do	
	printf -v j "%.2d" $i
	echo Running number $j...
	./recon 400 100
	
	cp data/recon.raw results/output_${j}.raw
	cp data/frecon.raw results/foutput_${j}.raw
	cp PHASING.log logs/iter_${j}.log
	cp prtf.dat logs/prtf_${j}.dat
	
	printf "output_%.2d.raw\t(400,100)\t0.65\t4.5 voxel support\n" $i >> results/KEY
	echo
done
