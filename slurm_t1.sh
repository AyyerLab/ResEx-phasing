#!/bin/sh

if [ $# -lt 1 ]
then
	echo Need folder name/number
	exit
fi
echo Saving to data/results_par_search_4et8_full_2mFo-DFc-sym_run_sequence/${1}/

cfile_iter='4et8_config_run_sequence_iter.ini'
cfile_DM='4et8_config_run_sequence_DM.ini'
cfile_ER='4et8_config_run_sequence_ER.ini'
fsc_file_DM='data/recon_liq-mod-0.8-100-outer150-bgfit-manyDM-ERDM/FSC-DM.txt'
fsc_file_ER='data/recon_liq-mod-0.8-100-outer150-bgfit-manyDM-ERDM/FSC-ER.txt'
fsc_value='0.400'
dirname_DM="data/recon_liq-mod-0.8-100-outer150-bgfit-manyDM-ERDM/DM/$1/"
dirname_ER="data/recon_liq-mod-0.8-100-outer150-bgfit-manyDM-ERDM/ER/$1/"

ls > /dev/null
ls .. > /dev/null
modulecmd bash load mpi/openmpi-x86_64
export OMP_NUM_THREADS=`nproc`

mkdir $dirname_DM
mkdir $dirname_ER

sed -i "5s/.*/num_threads = $OMP_NUM_THREADS/g" $cfile_iter
./recon -c $cfile_iter

sed -i "5s/.*/num_threads = $OMP_NUM_THREADS/g" $cfile_DM
./recon -c $cfile_DM

#calculate FSC
./utils/calc_fsc ./data/4et8_full_2mFo-DFc-srecon.raw ./data/recon_avg_2/4et8_full-pd.raw 301 2

#use Python Script to create plots
python plot_fsc_prtf.py

#write fsc at fsc_value to a file against the number of iterations
awk -v val="$fsc_value" -v OFS='\t' -v iter_index="$1" '$1 == val {print iter_index, $0}' fsc-0-0.dat >> $fsc_file_DM

#Cleanup and Organize Files
mv fsc_plot.pdf prtf_plot.pdf $dirname_DM
cp data/recon_avg_2/4et8_full-log.dat $cfile_iter $cfile_DM data/recon_avg_2/4et8_full-last.raw fsc-0-0.dat data/recon_avg_2/4et8_full-prtf.dat $dirname_DM

#Same steps for ER
sed -i "5s/.*/num_threads = $OMP_NUM_THREADS/g" $cfile_ER
./recon -c $cfile_ER

#calculate FSC
./utils/calc_fsc ./data/4et8_full_2mFo-DFc-srecon.raw ./data/recon_avg_2/4et8_full-pd.raw 301 2

#use python to plot
python plot_fsc_prtf.py

#write fsc at fsc_value to a file against the number of iterations
awk -v val="$fsc_value" -v OFS='\t' -v iter_index="$1" '$1 == val {print iter_index, $0}' fsc-0-0.dat >> $fsc_file_ER

#Cleanup and Organize Files
mv fsc_plot.pdf prtf_plot.pdf $dirname_ER
cp data/recon_avg_2/4et8_full-log.dat $cfile_iter $cfile_ER data/recon_avg_2/4et8_full-last.raw fsc-0-0.dat data/recon_avg_2/4et8_full-prtf.dat $dirname_ER
