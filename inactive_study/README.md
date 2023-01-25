# This directory hosts a Shannon Entropy study of the system to determine the number of inactive batches necessary for converging the fission distribution and eigenvalue keff so that it is free of contamination from the initial guess
In the command below, `-i` is the argument to python OpenMC model (w/o the `.py`) and `-input` is the argument to provide the Cardinal input file to run
`python inactive_study.py -i shannon_study_1000_elem -input openmc_1000.i -n-threads=16`
