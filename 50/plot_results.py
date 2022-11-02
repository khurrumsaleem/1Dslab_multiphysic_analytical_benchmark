import openmc
import numpy as np
import matplotlib.pyplot as plt

# Load statepoint file and get the number of binss
sp = openmc.StatePoint('statepoint.500.h5')
tally = sp.get_tally()
bins = tally.shape[0] # this works because the mesh is N x 1 x 1, so the total number of bins is N
# get scores for each tally type, compute realtive errors
# flux
flux = tally.get_slice(scores=['flux'])
# Determine relative error
flux_relative_error = np.zeros_like(flux.std_dev)
nonzero = flux.mean > 0
flux_relative_error[nonzero] = flux.std_dev[nonzero] / flux.mean[nonzero]
# distribution of relative errors
plt.hist(flux_relative_error[nonzero], bins=50)
plt.title('Relative Error Distribution for Flux')
plt.xlabel('Relative Error')
plt.ylabel('Frequency')
plt.savefig('relative_errors_50.png')
plt.clf()
# kappa fisison
kappa_fission = tally.get_slice(scores=['kappa-fission'])

# plot flux vs position
L = 10 # slab length
# plot using histogram
flux.mean.shape = (bins)
flux.std_dev.shape= (bins)
space_bins = np.linspace(-L,L,bins)
# mean,bin_edges = np.histogram(flux.mean,bins=bins)
plt.step(space_bins,flux.mean)
plt.title('Flux in Each Mesh Element')
plt.xlabel('Position')
plt.ylabel('Flux')
plt.savefig('flux_dist_50.png')