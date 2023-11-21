# 1Dslab_multiphysic_analytical_benchmark

Validation of Greisheimer and Kooreman Analytical Benchmark from PHYSOR 2022

This repository aims to use Cardinal to benchmark against a Greisheimer and Kooreman paper, "Analytical Benchmark Solution for 1D Neutron Transport Coupled with Themral Conduction and Material Expansion." (2022).

A PDF of this paper is included in this repository. Additionally the executable that was used for this simulation is located at `cardinal-opt`, though it isn't guaranteed to work in general, and in fact might not work on some operating systems/different versions of some packages. Best to make your own `cardinal-opt` that takes advantage of my [patched version of OpenMC on my fork](https://github.com/lewisgross1296/openmc/tree/s2patch). (branch named s2patch) This can be accomplished using the Cardinal environment variable `OPENMC_DIR`.

## There are various subdirectories for the cases as well as some other important aspects of modeling this problem
* All of the simulations are located in `simulations`. There are three subdirectories
    - `more_active_batches` which contains the most recent versions of the simulations. It has mesh sizes from 5 elements to 1000 elements
        * An impoprtant subdirectory of `more_active_batches` is named `restarts`. In it, each case reruns a transport simulation with the converged temperature distribution from the multiphysics simulation with more particles to get better statistics. The `restarts` directory has a subdirectory named `reviz` for visualizing the outpupt of these simulations.
    - `no_ss_no_mesh_template_flux` is an older directory that was used while some of the simulation settings were being determined. It has mesh sizes between 50 and 1000 elements.
    - `visualization` which houses some visualizations for the `more_active_batches` simulation and the `no_ss_no_mesh_template_flux_relaxation_on` subdirectories
* The `test_tracks` directory is for verifying the particle tracks are obeying the S2 patch as expected (i.e. only move in the +/- x direction).
* The `inactve_study` directory has the capability to examine Shannon entropy to help make decisions about the convergence of the fission source and number of required inactive batches.

Feel free to reach out with questions or to open an issue for any problems you notice!