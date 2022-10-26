# 1Dslab_multiphysic_analytical_benchmark
Validation of Greisheimer and Kooreman Analytical Benchmark from PHYSOR 2022

This repository aims to use Cardinal to benchmark against a Greisheimer and Kooreman paper, "Analytical Benchmark Solution for 1D Neutron Transport Coupled with Themral Conduction and Material Expansion." (2022).

The a PDF of the paper is included in this repository.

The subdirectories with numbers correspond to the number of x elements in each
simulation. The test_tracks directory is for verifying the S2 patch is as expected.

# In order to run the simulation
## create the OpenMC XML files for Cardinal
`python make_openmc_model.py`
## Generate the mesh
`~/cardinal/cardinal-opt -i mesh.i --mesh-only`
## run cardinal simulation using (optional) mpi paralellism and openmp parallelism from either of the following commands
### entire simulation
`mpiexec -np 4 ~/cardinal/cardinal-opt -i solid.i --n-threads=10`
### just openmc subapp
`mpiexec -np 4 ~/cardinal/cardinal-opt -i openmc.i --n-threads=10`