#!/bin/sh

#SBATCH -p cfel
#SBATCH -t 18:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

#SBATCH -o .%j.out
#SBATCH -e .%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=moshir.harsh@desy.de

#SBATCH -J ResEx_sim

cfile='4et8_config-7.ini'
fsc_file='data/recon_liq_0.8-100-sym-outer150-bgfit-randstartDM/FSC.txt'
fsc_value='0.400'

for i in `seq 1 20`;
do  
	ls > /dev/null
	ls .. > /dev/null
	modulecmd bash load mpi/openmpi-x86_64
	export OMP_NUM_THREADS=`nproc`

	dirname="data/recon_liq_0.8-100-sym-outer150-bgfit-randstartDM/$i/"
	mkdir $dirname

	sed -i "5s/.*/num_threads = $OMP_NUM_THREADS/g" $cfile
	./recon -c $cfile

	#calculate FSC
	./utils/calc_fsc ./data/4et8_full_2mFo-DFc-srecon.raw ./data/recon_temp7/4et8_full-pd.raw 301 2

	#use Python Script to create plots
	python plot_fsc_prtf.py

	#write fsc at fsc_value to a file against the number of iterations
	awk -v val="$fsc_value" -v OFS='\t' -v iter_index="$i" '$1 == val {print iter_index, $0}' fsc-0-0.dat >> $fsc_file

	#Cleanup and Organize Files
	mv fsc_plot.pdf prtf_plot.pdf $dirname
	cp data/recon_temp7/4et8_full-log.dat $cfile data/recon_temp7/4et8_full-last.raw fsc-0-0.dat data/recon_temp7/4et8_full-prtf.dat $dirname
done

