#!/bin/zsh

# This script takes the output from a PHENIX Phaser or Refinement run
# and produces a CCP4 map of a single asymmetric unit using cut_out_density
# Parameters:
# 	folder: Path to PHENIX output directory
# 	folder_type: The string 'phaser' or 'refine'
# 	res_cutoff: High resolution cutoff in Angstroms
# 	tag: Prefix of output file
# Output:
# 	<tag>_2mFo-DFc.ccp4 file in the folder from which the script is launched

if [ $# -lt 4 ]
then
	echo Format: $0 folder folder_type res_cutoff tag
	exit
fi

cd $1
folder_type=$2
res_cutoff=$3
if [ $folder_type = 'phaser' ]
	phenix.cut_out_density high_resolution=3.0 cutout_type=model \
		cutout_model_radius=${res_cutoff} *phaser*.pdb *phaser*.mtz
elif [ $folder_type = 'refine' ]
	phenix.cut_out_density high_resolution=2.5 cutout_type=model \
		cutout_model_radius=${res_cutoff} *refine*.pdb `ls *refine*.mtz|grep -v data`
else
	echo Unknown folder_type, $2
	exit
fi

phenix.fft prefix=$4 scale=volume cutout.mtz
cp ${4}*.ccp4 $OLDPWD
cd -
ls -l ${4}*.ccp4
