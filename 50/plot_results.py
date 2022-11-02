import openmc
import numpy as np
import matplotlib.pyplot as plt

# Load statepoint file and get the number of bins
sp = openmc.StatePoint('statepoint.500.h5')
tally = sp.get_tally()
num_voxels = tally.shape[0] # this works because the mesh is N x 1 x 1, so the total number of bins is N

# get scores for each flux
flux = tally.get_slice(scores=['flux'])

L = 106.47 # slab length from paper
# plot using hist
flux.mean.shape = (num_voxels)
flux_min = min(flux.mean)
flux_max = max(flux.mean)
flux.std_dev.shape= (num_voxels)
bin_boundaries = np.linspace(-L,L,num_voxels+1)
# give one value in each bin paired with the bin boundaries.
# then scale the weight of each bin by the flux
plt.hist(bin_boundaries[:-1],bin_boundaries,weights=flux.mean,histtype='step') # plot flux vs position
# TODO add error bars for flux bins
bin_average_values = (bin_boundaries[1:] + bin_boundaries[:-1])/2 # get bin centers
plt.title('Flux in Each Mesh Element')
plt.xlabel('Position [cm]')
plt.ylabel('Flux [n/cm^2-s]')
plt.xlim([-(L+2),L+2])
plt.xticks(bin_boundaries[::5])
plt.ylim([flux_min-0.01,flux_max+0.01])
plt.savefig('flux_dist_50.png')
plt.clf()
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

# TODO compute error from analytical solutioin
