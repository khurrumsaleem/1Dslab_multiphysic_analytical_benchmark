# This directory is for mesh size of 50 x elements
## create the OpenMC XML files for Cardinal
`python make_openmc_model_50.py`
## Generate the mesh
`~/cardinal/cardinal-opt -i mesh_50.i --mesh-only`
## run cardinal simulation using (optional) mpi paralellism and openmp parallelism from either of the following commands
### entire simulation
`mpiexec -np 4 ~/cardinal/cardinal-opt -i openmc_50_1e-4.i --n-threads=10`
### just solid subapp
`mpiexec -np 4 ~/cardinal/cardinal-opt -i solid_50.i --n-threads=10`

### On obtaining the flux results
Run the script via `python plot_results.py`. Note now that triggers have been turned on, the statepoint file used to plot the flux will need to be manually entered (unless I put in a command line argument for it sometime later)