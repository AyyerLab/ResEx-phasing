#!/usr/bin/env bash
####!/usr/bin/env zsh

if [ $# -lt 3 ]
then
	echo Format: $0 map size res_at_edge
	echo Optional: skip supp_radius supp_thresh point_group
	exit
fi

size=$2
res_at_edge=$3

# If zsh 
#rad=$(( size / 2 ))
#res=$(( 1. * res_at_edge * rad ))

# If bash
export size
rad=`echo|awk -v size=$size '{print int(size/2)}'`
export rad res_at_edge
res=`echo |awk '{print ENVIRON["res_at_edge"] * ENVIRON["rad"]}'`

supp_radius=3
supp_thresh=0.1
skip=0
point_group=222

# Parse optional arguments
if [ $# -gt 3 ]
then
	skip=$4
fi

if [ $# -gt 4 ]
then
	supp_radius=$5
	supp_thresh=$6
	echo supp_radius = $supp_radius supp_thresh = $supp_thresh
fi

if [ $# -gt 6 ]
then
	point_group=$7
	echo point_group = $point_group
fi

if [ $skip -eq 0 ]
then
	echo Doing entire pipeline
elif [ $skip -lt 2 ]
then
	echo Only doing support wrapping
fi

# Calculate filenames
map_name=`basename $1`
mapnoext="${map_name%.*}"
log_name=results/${mapnoext}.log
echo Log file: $log_name
strmodel=data/convert/${mapnoext}-str.cpx
lowresmodel_name=data/convert/${mapnoext}-recon.raw
supp_name=data/convert/${mapnoext}.supp
supprecon_name=data/convert/${mapnoext}-srecon.raw
cpx_name=data/convert/${mapnoext}.cpx

# Read map and generate 3D model
if [ $skip -eq 0 ]
then
	echo $1 > $log_name
	echo ./utils/read_map $1 $res | tee $log_name
	./utils/read_map $1 $res >> $log_name 2>&1
	
	padmodel=`grep "Saving padded" $log_name|awk '{print $5}'`
	padsize=`grep "volume size" $log_name|awk '{print $5}'`
	padnoext="${padmodel%.*}"
	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/gen_fdens $padmodel $padsize ${padnoext}-fdens.cpx 1| tee -a $log_name
	./utils/gen_fdens $padmodel $padsize ${padnoext}-fdens.cpx 1 >> $log_name 2>&1
	
	sx=`grep Stretch $log_name|awk -F'[=,()]' '{print $3}'`
	sy=`grep Stretch $log_name|awk -F'[=,()]' '{print $4}'`
	sz=`grep Stretch $log_name|awk -F'[=,()]' '{print $5}'`
	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/fstretch ${padnoext}-fdens.cpx $padsize $size $sx $sy $sz $strmodel | tee -a $log_name
	./utils/fstretch ${padnoext}-fdens.cpx $padsize $size $sx $sy $sz $strmodel >> $log_name 2>&1
	
	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/gen_dens $strmodel $size $lowresmodel_name | tee -a $log_name
	./utils/gen_dens $strmodel $size $lowresmodel_name >> $log_name 2>&1
	
	echo -------------------------------------------------------------------------------- >> $log_name
	echo Removing temporary files
	rm -v $padmodel ${padnoext}-fdens.cpx ${padnoext}-fdens-sym.raw ${strmodel} ${strmodel%.*}-sym.raw >> $log_name
fi

# Calculate support and constrain model by support
if [ $skip -lt 2 ]
then
	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/create_support $lowresmodel_name $size $supp_radius $supp_thresh $supp_name | tee -a $log_name
	./utils/create_support $lowresmodel_name $size $supp_radius $supp_thresh $supp_name >> $log_name 2>&1
	
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
	echo ./utils/gen_fdens $supprecon_name $size $cpx_name $point_group | tee -a $log_name
	./utils/gen_fdens $supprecon_name $size $cpx_name $point_group >> $log_name 2>&1
	
	echo -------------------------------------------------------------------------------- >> $log_name
	echo ./utils/create_support $supprecon_name $size $supp_radius $supp_thresh $supp_name | tee -a $log_name
	./utils/create_support $supprecon_name $size $supp_radius $supp_thresh $supp_name >> $log_name 2>&1
fi

