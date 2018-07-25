---
layout: default
---

# Tutorial
## Simulated Data

Here is a step-by-step guide to phasing some simulated data. This should be the
first thing you try before attempting to work with experimental data and will
provide a jumping-off point for exploring more options which may be required
with tricky real data.

### Process map
We will use the provided 4et8_sim.ccp4 map shipped with the code. This map of a
single lysozyme molecule was made from data collected at LCLS (Boutet et al.
Science 2011). Since this is a simulation, we will generate the intensity file
also from the map. The steps to do this are:
```
$ ./process_map.sh data/maps/4et8_sim.ccp4 301 2.0 0 3.0 2.0 222
$ cp data/convert/4et8_sim-sym.raw data/merges/4et8_intens.raw
```
The `process_map.sh` script parses the CCP4 map and produces a complex-valued
volume representing the Fourier transform of the molecule and a 3D support mask,
among a few other files.

### Set up system
### Run reconstruction
### Examine results
