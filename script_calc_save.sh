#!/bin/bash

# Script to run the Bragg_Count with the initial config file, plot the FSC and the PRTF and save FSC, PRTF, log and the config file

# Author: Moshir Harsh, July 26, 2017

# config_file = './4et8_config.ini'

#gedit ./4et8_config.ini

#echo "Enter the Folder name"

#read folder_name

./recon -c 4et8_config.ini

#Calculate FSC

./utils/calc_fsc ./data/4et8_full_2mFo-DFc-srecon.raw ./data/recon/4et8_full-pd.raw 301 2

#Use Python Script to create plots
python plot_fsc_prtf.py

#$base = "/data/data_different_parameters/" 
#$path = $base$folder_name

#echo $path

mkdir data/results_par_search_4et8_full_2mFo-DFc-P1/17
mv fsc_plot.svg prtf_plot.svg data/results_par_search_4et8_full_2mFo-DFc-P1/17
cp data/recon/4et8_full-log.dat 4et8_config.ini data/results_par_search_4et8_full_2mFo-DFc-P1/17

#mkdir folder_name
#mv fsc_plot.svg prtf_plot.svg folder_name
#cp data/recon/4et8_full-log.dat 4et8_config.ini folder_name




