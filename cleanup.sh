#!/bin/zsh

i=1
while [ $i -lt 200 ]
do
	j=`printf "%.2d" $i`
	if [ -f results/output_${j}.raw ]
	then
		((i++))
	else
		break
	fi
done

set -e
voxres=800.
if [ $# -gt 0 ]
then
	if [ $1 = '-v' ]
	then
		set -x
	else
		voxres=$1
	fi
fi


prefix=`grep output_prefix config.ini|awk -F"[ =]" '{print $2}'`
size=`grep size config.ini|awk -F"[ =]" '{print $2}'`
center=$((size/2))
voxsize=$((voxres / 2. / center))

echo Cleaning to number $j with data from $prefix
cp ${prefix}-recon.raw results/output_${j}.raw
cp ${prefix}-frecon.raw results/foutput_${j}.raw
cp ${prefix}-log.dat results/log_${j}.dat
cp ${prefix}-prtf.dat results/prtf_${j}.dat
cp config.ini results/conf_${j}.ini
echo Generating map
supp_fname=`grep support_fname config.ini|awk -F"[ =]" '{print $2}'`
echo support = $supp_fname
echo ./utils/gen_map results/output_${j}.raw $size $voxsize $voxsize $voxsize $supp_fname
./utils/gen_map results/output_${j}.raw $size $voxsize $voxsize $voxsize $supp_fname
cd data/maps
rm -f output_${j}.mtz
phenix.map_to_structure_factors output_${j}.map.ccp4 box=True output_file_name=output_${j}.mtz
cd -
