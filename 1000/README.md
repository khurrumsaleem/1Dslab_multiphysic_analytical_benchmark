# This directory is for mesh size of 1000 x elements
## create the OpenMC XML files for Cardinal
`python make_openmc_model_1000.py`
## Generate the mesh
`~/cardinal/cardinal-opt -i mesh_1000.i --mesh-only`
## run cardinal simulation using (optional) mpi paralellism and openmp parallelism from either of the following commands
### entire simulation
`mpiexec -np 4 ~/cardinal/cardinal-opt -i solid_1000.i --n-threads=10`
### just openmc subapp
`mpiexec -np 4 ~/cardinal/cardinal-opt -i openmc_1000.i --n-threads=10`