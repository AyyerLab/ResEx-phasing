#!/bin/sh

if [ $# -lt 1 ]
then
	echo Need folder name/number
	exit
fi
echo Saving to data/results_par_search_4et8_full_2mFo-DFc-sym_run_sequence/${1}/

#SBATCH -p cfel
#SBATCH -t 6:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

#SBATCH -o .%j.out
#SBATCH -e .%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=moshir.harsh@desy.de

#SBATCH -J ResEx_sim
#cfile_iter='4et8_config_run_sequence_iter.ini'
c_file_DM='4et8_config_run_sequence_DM.ini'
c_file_ER='4et8_config_run_sequence_ER.ini'
fsc_file_DM='FSC_run_sequence_DM.txt'
fsc_file_ER='FSC_run_sequence_ER.txt'
fsc_value='0.400'
dirname_DM="data/results_par_search_4et8_full_2mFo-DFc-sym_run_sequence_DM/${1}/"
dirname_ER="data/results_par_search_4et8_full_2mFo-DFc-sym_run_sequence_ER/${1}/"

ls > /dev/null
ls .. > /dev/null
modulecmd bash load mpi/openmpi-x86_64
export OMP_NUM_THREADS=`nproc`

mkdir $dirname_DM
mkdir $dirname_ER

#sed -i "5s/.*/num_threads = $OMP_NUM_THREADS/g" $cfile_iter
#./recon -c $cfile_iter

sed -i "5s/.*/num_threads = $OMP_NUM_THREADS/g" $cfile_DM
./recon -c $cfile_DM

#calculate FSC
./utils/calc_fsc ./data/4et8_full_2mFo-DFc-srecon.raw ./data/recon_avg/4et8_full-pd.raw 301 2

#use Python Script to create plots
python plot_fsc_prtf.py

#write fsc at fsc_value to a file against the number of iterations
awk -v val="$fsc_value" '$1 == val' fsc-0-0.dat >> $fsc_file_DM

#Cleanup and Organize Files
mv fsc_plot.svg prtf_plot.svg $dirname_DM
cp data/recon_avg/4et8_full-log.dat 4et8_config.ini data/recon_avg/4et8_full-last.raw fsc-0-0.dat $dirname_DM

#Same steps for ER
sed -i "5s/.*/num_threads = $OMP_NUM_THREADS/g" $cfile_ER
./recon -c $cfile_ER

#calculate FSC
./utils/calc_fsc ./data/4et8_full_2mFo-DFc-srecon.raw ./data/recon_avg/4et8_full-pd.raw 301 2
python plot_fsc_prtf.py

#write fsc at fsc_value to a file against the number of iterations
awk -v val="$fsc_value" '$1 == val' fsc-0-0.dat >> $fsc_file_ER

#Cleanup and Organize Files
mv fsc_plot.svg prtf_plot.svg $dirname_ER
cp data/recon/4et8_full-log.dat 4et8_config.ini data/recon/4et8_full-last.raw fsc-0-0.dat $dirname_ER

