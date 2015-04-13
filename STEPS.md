# Reconstruction pipeline
* * *

For the phased-Bragg resolution-extension project, the following are the steps 
needed to perform the reconstruction:

* Generate hkl grid of complex Fourier amplitudes
* Generate input high-q merged intensities
* Generate support file
* Modify config file
* Run reconstruction
* View results

## Generation of hkl values:
This is the procedure to get the complex hkl grid which will be used at low 
resolutions.

### utils/fstretch
Stretch the true-scale simulations to take into account the integer Bragg 
spacing. The stretch factors are (2.2045, 1.9475, 1.913).
	./utils/fstretch ../partial/data/out3d_1280c.cpx 1280 501 2.2045 1.9475 1.913

This produces `../partial/data/out3d_1280c_str_501.cpx` which is a 501^3 volume.

### utils/bragg_gen
Generate Bragg hkl values by sampling the stretched complex intensities. The 
reciprocal lattice constants are (6, 4, 3)
	./utils/bragg_gen ../partial/data/out3d_1280c_str_501.cpx 501 6 4 3

This generates `data/hkl_167...cpx`

## Generation of merged data
We start with a 501^3 volume of real numbers. The main pre-processing step is 
to set low-q values to -1. and cube-corner values to 0. 

### utils/zero_outer

## Generate support file

## Modify config file

## Run reconstruction
After the config file has been processed, the reconstruction is run using
	./recon <iter> <ave_iter>
where `<iter>` is the number of iterations and `<ave_iter>` is the iteration 
number after which it starts a running average and calculates the PRTF. The 
log file is PHASING.log

## View results
The reconstruction produces data/recon.raw and its intensities data/frecon.raw.
Every 10 iterations, it also write the intermediate state data/interm.raw.

All of these can be viewed using the Python script view.py. The usage is:
	./view.py data/recon.raw
	./view.py data/frecon.raw 1
The extra flag is added if one wants to look at the full volume and not just 
the central region.
