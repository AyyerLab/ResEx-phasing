---
layout: default
---

# Documentation

Welcome to the detailed documentation page. Here you will find algorithmic
details of what the reconstruction program does, descriptions of all the options
you can pass in the configuration file, as well as a list of the various utility
scripts.

We recommend the user goes through the [basic tutorial](sim-tutorial.md) first
in order to get a feel for 'what' the program does. This page discusses more the
'how' of the program.

## Table of Contents

 * [Algorithmic details](#algorithm-details)
   * [Basic Projections](#basic-projections)
   * [Two data sets](#two-data-sets)
   * [Histogram matching](#histogram-matching)
   * [Dynamic support update (solvent flattening)](#dynamic-support)
   * [Spherically symmetric background fitting](#background-fitting)
   * [Composing the projections](#composing-projections)
 * [Configuration options](#config-options)
   * [Parameters](#parameters)
   * [Files](#files)
   * [Algorithm](#algorithm)
     * [Iteration types](#iteration-types)
	 * [Extra flags](#extra-flags)
 * [Utilities](#utilities)

## Algorithm Details

The reconstruction is done by a family of phase retrieval methods known as 
iterative projection algorithms. We start with a random initial condition and 
incrementally improve it by applying some mapping to it. The "projection" part
means that each mapping is composed of projection operations. A projection to a 
set with a given input returns the member of the set closest to the input, where
'closest' is defined under some metric (usually Euclidean).

In conventional phase retrieval performed in coherent diffractive imaging (CDI),
two projections are one which makes the model consistent with the data and the
other makes it consistent with a support constraint which says that the density
must be non-zero only in a small region. These two projection operations are
then composed together in special ways to improve convergence properties.

In our experiment, we have conventional Bragg data at low resolution and twinned
continuous diffraction beyond that. To start with, the two data sets are treated
independently, with the Bragg diffraction phased using conventional methods like
molecular replacement, and the diffuse data is untwinned and phased to produce a
higher resolution model.

### Basic Projections
### Two data sets
### Histogram Matching
### Dynamic Support
### Background fitting
### Composing the projections

## Configuration options

The implementation of the various algorithmic details are implemented via the 
configuration file. Here all the possible options are listed and described. The
file is split into three sections, `[parameters]`, `[files]` and `[algorithm]'.
The first two set up the system and the third describes what the program does.

### [parameters]
### [files]
### [algorithm]
#### Iteration types
#### Extra flags

## Utilities

In addition to the main reconstruction code, quite a number of helper utilities
are provided. They can be used to prepare or modify the input data or analyze
the reconstruction output. Additionally, there are some convenience scripts
which employ one or more utilities.
