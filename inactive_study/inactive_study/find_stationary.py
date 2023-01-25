import numpy as np
import openmc
import sys
import statistics as stat

# eigenvalue batching
n_inactive = 250
n_active = 50

n = 3
sp_file = 'statepoint_' + str(n) + '_layers.h5'

entropy = []
k = []
k_generation_avg = []
std_devs = []

max_entropy = sys.float_info.min
min_entropy = sys.float_info.max
max_k = sys.float_info.min
min_k = sys.float_info.max
averaging_batches = 10

with openmc.StatePoint(sp_file) as sp:
    entropy.append(sp.entropy)
    k.append(sp.k_generation)
    max_entropy = max(np.max(sp.entropy[:n_inactive][:]), max_entropy)
    min_entropy = min(np.min(sp.entropy[:n_inactive][:]), min_entropy)
    max_k = max(np.max(sp.k_generation[:n_inactive][:]), max_k)
    min_k = min(np.min(sp.k_generation[:n_inactive][:]), min_k)

    averaging_k = sp.k_generation[(n_inactive - averaging_batches):n_inactive]
    k_generation_avg.append(sum(averaging_k) / len(averaging_k))

# detecting stationarity. Accordinig to Brown, the below is done in MCNP5:
# find the first cycle when Hsrc is within one standard deviation of its average for the last half of cycles
# TODO decide whether to use this, all points, or a windowed average
window_length = 10
# detect when entropy of current batch is within window (mean +/- std of window)
for j in range(window_length,len(entropy[0])):
    window = entropy[0][(j-window_length):j]
    window_mean = np.average(window)
    window_dev = np.std(window)
    window_low = window_mean - window_dev
    window_high = window_mean + window_dev
    # print("window bounds: ", window_low , " to ", window_high )
    if(window_high - entropy[0][j] > 0 and entropy[0][j]-window_low > 0):
        # if entropy[j] is less than window_high and greater than
        # window_low, then it is within the window and we've succeeded
        print("Iteration ", j, " produced a value in the window. It is recommended "
              "to do at least this many inactive cycles.")
        break
    else:
        continue

# compute on the fly the average of the eigenvalue (last half batches) +/- the std dev of the last half cycles
# and see if the current k is within this interval
