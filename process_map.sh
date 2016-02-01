#!/bin/zsh

if [ $# -lt 2 ]
then
	echo Format: $0 map size res_at_edge
	exit
fi

size=$2
res_at_edge=$3
rad=$(( size / 2 ))
res=$(( 1. * res_at_edge * rad ))

map_name=`basename $1`
mapnoext="${map_name%.*}"
log_name=${mapnoext}.log

echo $log_name
echo ./utils/read_map $1 $res
./utils/read_map $1 $res &>$log_name

sx=`grep Stretch $log_name|awk -F'[=,()]' '{print $3}'`
sy=`grep Stretch $log_name|awk -F'[=,()]' '{print $4}'`
sz=`grep Stretch $log_name|awk -F'[=,()]' '{print $5}'`
padmodel=`grep "Saving padded" $log_name|awk '{print $5}'`
padsize=`grep "volume size" $log_name|awk '{print $5}'`
padnoext="${padmodel%.*}"
strmodel=data/${mapnoext}-str.cpx
lowresmodel_name=data/${mapnoext}-recon.raw
supp_name=data/${mapnoext}-3.supp

echo ./utils/gen_fdens $padmodel $padsize
#echo ================================================================================ >> $log_name
./utils/gen_fdens $padmodel $padsize &>> $log_name

echo ./utils/fstretch ${padnoext}-fdens.cpx $padsize $size $sx $sy $sz $strmodel
#echo ================================================================================ >> $log_name
./utils/fstretch ${padnoext}-fdens.cpx $padsize $size $sx $sy $sz $strmodel &>> $log_name

echo ./utils/gen_dens $strmodel $size $lowresmodel_name
#echo ================================================================================ >> $log_name
./utils/gen_dens $strmodel $size $lowresmodel_name &>> $log_name

echo ./utils/create_support $lowresmodel_name $size 3. 5 $supp_name
#echo ================================================================================ >> $log_name
./utils/create_support $lowresmodel_name $size 3. 5 $supp_name &>> $log_name

