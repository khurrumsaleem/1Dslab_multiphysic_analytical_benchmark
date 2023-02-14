import numpy as np
import matplotlib.pyplot as plt

# demo block
L = 100
num_voxels = 4
fluxes = np.array([1,2,3,2]) # true values
bin_boundaries = np.linspace(-L,L,num_voxels+1)
# plot
plt.hist(bin_boundaries[:-1],bin_boundaries,weights=fluxes,bottom=0,histtype='step')
plt.title('Flux in Each Voxel')
plt.xlabel('x Position')
plt.ylabel('Flux')
plt.savefig('mwe.png')
plt.clf()