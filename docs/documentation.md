---
layout: default
---

<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS_CHTML"></script>

# Documentation [Under construction]

Welcome to the detailed documentation page. Here you will find algorithmic details of what the reconstruction program does, descriptions of all the options you can pass in the configuration file, as well as a list of the various utility scripts.

We recommend the user goes through the [basic tutorial](sim-tutorial.html) first in order to get a feel for 'what' the program does. This page discusses more the 'how' of the program.

## Table of Contents

 * [Algorithmic details](#algorithm-details)
   * [Basic Projections](#basic-projections)
   * [Two data sets](#two-data-sets)
   * [Histogram matching](#histogram-matching)
   * [Dynamic support update (solvent flattening)](#dynamic-support)
   * [Spherically symmetric background fitting](#background-fitting)
   * [Composing the projections](#composing-the-projections)
 * [Configuration options](#configuration-options)
   * [Parameters](#parameters)
   * [Files](#files)
   * [Algorithm](#algorithm)
 * [Utilities](#utilities)

## Algorithm Details

The reconstruction is done by a family of phase retrieval methods known as iterative projection algorithms. We start with a random initial condition and incrementally improve it by applying some mapping to it. The "projection" part means that each mapping is composed of projection operations. A projection to a set with a given input returns the member of the set closest to the input, where 'closest' is defined under some metric (usually Euclidean).

In conventional phase retrieval performed in coherent diffractive imaging (CDI), two projections are one which makes the model consistent with the data and the other makes it consistent with a support constraint which says that the density must be non-zero only in a small region. These two projection operations are then composed together in special ways to improve convergence properties.

In our experiment, we have conventional Bragg data at low resolution and twinned continuous diffraction beyond that. To start with, the two data sets are treated independently, with the Bragg diffraction phased using conventional methods like molecular replacement, and the diffuse data is untwinned and phased to produce a higher resolution model.

### Basic Projections
### Two data sets
### Histogram Matching
### Dynamic Support
### Background fitting
### Composing the projections
The reconstruction algorithm works by composing the projections together in certain specific ways in order to update the current iterate. They are expressed in terms of a mapping by which the next iterate $$x_{n+1}$$ is calculated from the current iterate $$x_n$$ and its two projections $$P_D(x_n)$$ and $$P_F(x_n)$$. The following update rules are implemented:

 * `ER` (Error reduction or alternating projections): The simplest update rule, which takes you to the nearest local minima in the error.
   * <span>$$x_{n+1} = P_D[P_F(x_n)]$$</span>
 * `HIO` (Hybrid Input-Output): Relaxation of error reduction to only partially enforce the support constraint. Here it is expressed in the projection mapping form and extended to more general direct-space projections.
   * <span>$$x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right) P_F(x_n) - \frac{x_n}{\beta}\right] - P_F(x_n)\right\}$$</span>
 * `DM` (Difference map): Update rule avoids local minima and efficiently searches solution space. Sometimes has trouble with jumping around too much with noisy data.
   * <span>$$x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right)P_F(x_n) - \frac{x_n}{\beta}\right] - P_F\left[\left(1-\frac{1}{\beta}\right)P_D(x_n) + \frac{x_n}{\beta}\right]\right\}$$</span>
 * `mod-DM` (Modified difference map): Modification of difference map to stay closer to data projection of current iterate.
   * <span>$$x_n' = \beta x_n + (1-\beta) P_F(x_n)$$</span>
   * <span>$$x_{n+1} = x_n' + P_F\left[2 P_D(x_n') - x_n'\right] - P_D(x_n')$$</span>
 * `RAAR` (Relaxed averaged alternating reflections): Another modification of difference map to be more noise-tolerant, originally for the case where the direct space constraint is support plus positivity.
   * <span>$$x_{n+1} = \beta \left\{x_n + P_D\left[2 P_F(x_n)\right] + P_D\left[-x_n\right] - P_F(x_n)\right\} + (1-\beta) P_F(x_n)$$</span>
   * If $$P_D$$ can be assumed to be linear, then the update rule simplifies to
     * <span>$$x_{n+1} = \beta \left\{x_n + P_D\left[2 P_F(x_n) - x_n\right] - P_F(x_n)\right\} + (1-\beta) P_F(x_n)$$</span>

All the update rules except `ER` are equivalent if $$\beta = 1$$. 

## Configuration options

The implementation of the various algorithmic details are implemented via the  configuration file. Here all the possible options are listed and described. The file is split into three sections, `[parameters]`, `[files]` and `[algorithm]`. The first two set up the system and the third describes what the program does.

### [parameters]

| `size` | Edge length of the voxel grid. Currently, only cubic grids are supported |
| `bragg_qmax` | fraction of edge resolution to which the low-Q Bragg data should be used. If the edge resolution is 2 &#8491; and the Bragg data is trustworthy up to 5 &#8491;, the value will be 2/5 = 0.4 |
| `scale_factor` | Relative scaling between the Bragg and diffuse Fourier magnitudes. Conventionally calculated using the `calc_scale` utility. |
| `num_threads` | Number of threads to use in the multi-threaded FFT and other parallelizable sub-routines. |
| `point_group` | Point group of the diffuse intensities. Currently supported options are '1', '222' and '4' |

### [files]
The following tokens are all paths to various files. Parentheses indicate optional values.

| `intens_fname` | Diffuse intensity volume. Raw `size`<sup>3</sup> grid of float32 values. |
| `bragg_fname` | Complex valued Fourier amplitudes of low-resolution model. Raw `size`<sup>3</sup> grid of complex64 values. note that only values in a central sphere of radius `floor(bragg_qmax * size / 2)` need to be meaningful. |
| `support_fname` | Support mask defining region of real-space where the density can be non-zero. Raw `size`<sup>3</sup> grid of uint8 values. |
| `output_prefix` | Prefix to output files. This is usually a path plus the initial part of a filename. For example, the log file will be at `<output_prefix>-log.dat` |
| (`input_fname`) | Optional path to initial guess for electron densities. If you want to continue a previous reconstruction pass in the `<output_prefix>-last.raw` file. Raw `size`<sup>3</sup> grid of float32 values. |
| (`inputbg_fname`) | Optional path to initial guess for background intensities. Only used if the `bg_fitting` option is turned on. Raw `size`<sup>3</sup> grid of float32 values. |

### [algorithm]
As before parentheses indicate optional values. Flags can take the value `0` for off/false and `1` for on/true. All flags are off by default.

| `algorithm` | Space-separated list of alternating numbers and algorithm codes. Code options are `ER`, `DM`, `HIO`, `RAAR`, `mod-DM`. See [description](#composing-the-projections) above for details of what they do. For example `100 DM 50 ER 100 DM` will do 100 iterations of difference map followed by 50 iterations of error reduction, followed by 100 more iterations of difference map. |
| `avg_algorithm` | Same format as in the `algorithm` option, but a running average of the projections is calculated over these iterations. |
| `beta` | Relaxation parameter for all the algorithms above except `ER`. Note that for `beta = 1`, all algorithms except `ER` are equivalent. |
| (`bg_fitting`) | FLAG for whether to retrieve a spherically symmetric background in addition to the real-space model. See [description](#background-fitting) above for details. |
| (`local_variation`) | FLAG for whether to dynamically update the support using a local variation threshold *a la* solvent flattening. If `0`, the support stays constant. |
| (`positivity`) | FLAG for whether to assume the electron densities in the support are non-negative. Usually not applicable since the zero value in protein crystallography is set to be the solvent level. |
| (`normalize_prtf`) | FLAG for whether to sharpen the output model by the reciprocal of the PRTF. The averaging process by which the PRTF is calculated reduces the fourier magnitudes. If this option is turned on, the output reconstruction is high-pass filtered to cancel out this effect. This makes the I vs q dependence match the diffuse intensities, which may be more suitable for refinement. However, it also enhances noise, so the right resolution cutoff must be chosen. |
| (`histogram`) | FLAG for whether to apply a histogram projection. Requires `hist_fname` to be set. See [description](#histogram-matching) above for details. |
| (`hist_fname`) | Path to two-column text file specifying the histogram to be matched against. The first line of the file is the number of rows. Following that, the first column has the value of the density in the center of the bin and the second column has the fraction of voxels with that value. |
| (`blurring`) | FLAG for whether to rotationally blur the calculated Fourier intensities before comparing with the data. Needs either `quat_fname` or both `num_div` and `sigma_deg`.
| (`quat_fname`) | Path to quaternion file in the [*Dragonfly*](http://github.com/duaneloh/Dragonfly/wiki) format. |
| (`num_div`) | Refinement level of quasi-uniform SO(3) rotation sampling. Number of samples of the whole rotation group is `10*(5*num_div^2 + num_div)`. Needs `sigma_deg` to also be specified. |
| (`sigma_deg`) | Width in degrees of Gaussian distribution of angles centered on the identity quaternion. Combined with `num_div` this generates a quaternion list, which can be used for rotational blurring |

## Utilities

In addition to the main reconstruction code, quite a number of helper utilities
are provided. They can be used to prepare or modify the input data or analyze
the reconstruction output. Additionally, there are some convenience scripts
which employ one or more utilities.
