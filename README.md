# 1Dslab_multiphysic_analytical_benchmark

Validation of Greisheimer and Kooreman Analytical Benchmark from PHYSOR 2022

This repository aims to use Cardinal to benchmark against a Greisheimer and Kooreman paper, "Analytical Benchmark Solution for 1D Neutron Transport Coupled with Themral Conduction and Material Expansion." (2022).

The a PDF of the paper is included in this repository.

# There are various subdirectories for the cases as well as some other important aspects
* The subdirectories with numbers correspond to each case with the directory number of x-elements in the mesh. The neutronics and MOOSE meshes have the same dimensions.
* The `test_tracks` directory is for verifying the particle tracks are obeying the S2 patch as expected (i.e. only move in the +/- x direction).
* The `test` directory was used to verify an OpenMC error related to this GitHub [Issue](https://github.com/openmc-dev/openmc/issues/2296) that was experienced when setting up this simulation.
* The `inactve_study` directory has the capability to examine Shannon entropy and make decisions about the convergence of the fission source and number of required inactive batches.
* The `visualization` directory is used for post processing results into plots for the paper