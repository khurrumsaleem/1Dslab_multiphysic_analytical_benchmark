# This directory is for mesh size of 50 x elements
## create the OpenMC XML files for Cardinal
`python make_openmc_model_50.py`
## Generate the mesh
`~/cardinal/cardinal-opt -i mesh_50.i --mesh-only`
## run cardinal simulation using (optional) mpi paralellism and openmp parallelism from either of the following commands
### entire simulation
`mpiexec -np 2 ~/cardinal/cardinal-opt -i openmc_50.i --n-threads=18`
### just solid subapp
`mpiexec -np 2 ~/cardinal/cardinal-opt -i solid_50.i --n-threads=18`