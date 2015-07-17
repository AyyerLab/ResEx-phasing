#!/bin/bash

i=1
while [ $i -lt 100 ]
do
	printf -v j "%.2d" $i
	if [ -f results/output_${j}.raw ]
	then
		((i++))
	else
		break
	fi
done

echo Cleaning to number $j...
#echo Give results message
#read comment
cp data/recon.raw results/output_${j}.raw
cp data/frecon.raw results/foutput_${j}.raw
cp PHASING.log logs/iter_${j}.log
cp prtf.dat logs/prtf_${j}.dat
#printf "output_%.2d.raw\t(600,400)\t0.6\t\tSame as \`29', but tight support\n" $i
#printf "output_%.2d.raw\t(600,400)\t0.6\t\tSame as \`29', but tight support\n" $i >> results/KEY
echo

