# This directory hosts a Shannon Entropy study of the system to determine the number of inactive batches necessary for converging the fission distribution and eigenvalue keff so that it is free of contamination from the initial guess
In the command below, `-i` is the argument to python OpenMC model (w/o the `.py`) and `-input` is the argument to provide the Cardinal input file to run
`python inactive_study.py -i shannon_study_1000_elem -input openmc_1000.i -n-threads=16`
To use stationary detection, go to `/inactive_study/inactive_study` and run `python find_stationary.py` providing command line arguments for the method to detect steady state. There are three options
- all: use all batches for assessing convergence `python find_stationary.py --method=all`
- half: use the last half of batches `python find_stationary.py --method=half`
- window: specify a window length via `--window_length` and asssess convergence over that many past entries `python find_stationary.py --method=window --window_length=25`