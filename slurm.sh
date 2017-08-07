#!/bin/sh

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
cfile='4et8_config.ini'

ls > /dev/null
ls .. > /dev/null
modulecmd bash load mpi/openmpi-x86_64
export OMP_NUM_THREADS=`nproc`

sed -i "5s/.*/num_threads = $OMP_NUM_THREADS/g" $cfile

./recon -c $cfile
