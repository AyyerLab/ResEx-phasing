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
cp PHASING.log results/log_${j}.dat
cp prtf.dat results/prtf_${j}.dat
cp src/config.ini results/conf_${j}.ini
#printf "output_%.2d.raw\t(600,400)\t0.6\t\tSame as \`29', but tight support\n" $i
#printf "output_%.2d.raw\t(600,400)\t0.6\t\tSame as \`29', but tight support\n" $i >> results/KEY
echo Generating map
supp_fname=`grep support_fname src/config.ini|awk '{print $3}'`
./utils/gen_map results/output_${j}.raw 301 300. $supp_fname
cd data/maps
rm -f output_${j}.mtz
phenix.map_to_structure_factors output_${j}.map.ccp4 d_min=2.3 output_file_name=output_${j}.mtz
cd ..
cd ..
pwd
