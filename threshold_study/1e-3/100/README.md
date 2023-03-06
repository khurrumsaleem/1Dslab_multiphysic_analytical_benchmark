# This directory is for mesh size of 100 x elements
## create the OpenMC XML files for Cardinal
`python make_openmc_model_100.py`
## Generate the mesh
`~/cardinal/cardinal-opt -i mesh_100.i --mesh-only`
## run cardinal simulation using (optional) mpi paralellism and openmp parallelism from either of the following commands
### entire simulation
`mpiexec -np 4 ~/cardinal/cardinal-opt -i openmc_100_1e-3.i --n-threads=10`
### just solid subapp
`mpiexec -np 4 ~/cardinal/cardinal-opt -i solid_100.i --n-threads=10`