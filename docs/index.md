---
layout: default
---

# ResEx Phasing
### Resolution Extension using Continuous Diffraction

ResEx is used for **RES**olution **EX**tension of a protein crystal structure using continuous diffraction data at higher resolution, where Bragg peaks are no longer visible. 

Here you will find usage instructions as well as other information about the program.

### Installation
There is no real installation process to speak of. Just clone the repository and compile the C code:
```
$ git clone https://github.com/kartikayyer/bragg-cont.git
$ cd bragg-cont
$ make
```
However, there are a few dependencies:
 * GCC (Other C compilers will also work, with minor modifications to the Makefile)
 * GNU Scientific Library (GSL) (for special functions and random number generation)
 * FFTW3 (for Fourier transforms. Single precision with threading)
 * python 2.7 (for the GUIs, preferably from Anaconda)
    * numpy
    * matplotlib

### Usage
A quick start page for a basic reconstruction from simulated data can be found [here]("{{ sim-tutorial.md | relative_url }}").

The inputs to the program are a 3D volume of oversampled intensities and a CCP4/MRC map containing a low-resolution model obtained from the Bragg data. There are helper scripts to extract this map from the output of a PHENIX run (CCP4 support coming soon).

