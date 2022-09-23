# 1Dslab_multiphysic_analytical_benchmark
Validation of Greisheimer and Kooreman Analytical Benchmark from PHYSOR 2022

This repository aims to use Cardinal to benchmark against a Greisheimer and Kooreman paper, "Analytical Benchmark Solution for 1D Neutron Transport Coupled with Themral Conduction and Material Expansion." (2022).

The a PDF of the paper is included in this repository.

# In order to run the simulation
## create the OpenMC XML files for Cardinal
`python make_openmc_model.py`
## run cardinl simulation using (optional) mpi paralellism and openmp parallelism
`mpiexec -np 4 ~/cardinal/cardinal-opt -i openmc.i --n-threads=10`