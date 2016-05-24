#!/bin/zsh

if [ $# -lt 3 ]
then
	echo Format: $0 map size res_at_edge
	exit
fi

skip=0
if [ $# -gt 3 ]
then
	skip=$4
fi

size=$2
res_at_edge=$3
rad=$(( size / 2 ))
res=$(( 1. * res_at_edge * rad ))

map_name=`basename $1`
mapnoext="${map_name%.*}"
log_name=results/${mapnoext}.log

sx=`grep Stretch $log_name|awk -F'[=,()]' '{print $3}'`
sy=`grep Stretch $log_name|awk -F'[=,()]' '{print $4}'`
sz=`grep Stretch $log_name|awk -F'[=,()]' '{print $5}'`
padmodel=`grep "Saving padded" $log_name|awk '{print $5}'`
padsize=`grep "volume size" $log_name|awk '{print $5}'`
padnoext="${padmodel%.*}"
strmodel=data/${mapnoext}-str.cpx
lowresmodel_name=data/${mapnoext}-recon.raw
supp_name=data/${mapnoext}-3.supp
supprecon_name=data/${mapnoext}-srecon.raw
supp_radius=3
supp_thresh=10

echo skip = $skip
echo Log file: $log_name

if [ $skip -eq 0 ]
then
	echo ./utils/read_map $1 $res | tee $log_name
	./utils/read_map $1 $res &>> $log_name

	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/gen_fdens $padmodel $padsize | tee -a $log_name
	./utils/gen_fdens $padmodel $padsize &>> $log_name

	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/fstretch ${padnoext}-fdens.cpx $padsize $size $sx $sy $sz $strmodel | tee -a $log_name
	./utils/fstretch ${padnoext}-fdens.cpx $padsize $size $sx $sy $sz $strmodel &>> $log_name

	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/gen_dens $strmodel $size $lowresmodel_name | tee -a $log_name
	./utils/gen_dens $strmodel $size $lowresmodel_name &>> $log_name
fi

if [ $skip -lt 2 ]
then
	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/create_support $lowresmodel_name $size $supp_radius $supp_thresh $supp_name | tee -a $log_name
	./utils/create_support $lowresmodel_name $size $supp_radius $supp_thresh $supp_name &>> $log_name

	echo -------------------------------------------------------------------------------- >> $log_name
	echo Constraining lowresmodel by support | tee -a $log_name
	python <<-EOF
	import numpy as np
	m = np.fromfile('$lowresmodel_name', '=f4')
	s = np.fromfile('$supp_name', '=u1')
	m *= s
	m.tofile('$supprecon_name')
	EOF

	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/gen_fdens $supprecon_name $size data/${mapnoext}.cpx | tee -a $log_name
	./utils/gen_fdens $supprecon_name $size data/${mapnoext}.cpx &>> $log_name

	echo -------------------------------------------------------------------------------- >> $log_name
	supp_name=data/${mapnoext}-$supp_radius.supp
	echo ./utils/create_support $lowresmodel_name $size $supp_radius $supp_thresh $supp_name | tee -a $log_name
	./utils/create_support $lowresmodel_name $size $supp_radius $supp_thresh $supp_name &>> $log_name
fi
