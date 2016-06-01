#!/bin/zsh

if [ $# -lt 2 ]
then
	echo Format: $0 folder tag
	exit
fi

cd $1
phenix.cut_out_density high_resolution=3.0 cutout_type=model cutout_model_radius=2.0 *phaser*.pdb *phaser*.mtz
phenix.fft prefix=$2 cutout.mtz
cp *.ccp4 ..
cd ..
