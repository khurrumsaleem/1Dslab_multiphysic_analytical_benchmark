import openmc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# list mesh element sizes that will be used to iterate through all cases
n_elems = [50,100,250,500,1000]

# mesh size vs error
# collect error norm data
norms = []
for n in n_elems:
    filename = 'openmc_' + str(n) + '_out.csv'
    df = pd.read_csv(filename)
    last_val  = df.iloc[len(df)-1,-1]
    norms.append(last_val)

# plot element size vs computed error
plt.plot(n_elems,norms,'-ro',label=r'$\epsilon=\frac{||T_{a}-T_{x}||_{2}}{||T_{a}||_{2}}$')
plt.xticks([0,50,100,250,500,1000])
plt.yticks(np.linspace(0.0010,0.0017,6))
# plt.yticks(np.linspace(0.0010,0.0017,5))
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('Number of X Mesh Elements')
plt.ylabel('Error Norm $\epsilon$')
plt.title('Error from Analytical Solution for Varying Mesh Elements',y=1.05)
plt.legend(bbox_to_anchor=[0.9,0.985])
plt.grid()
plt.savefig('error_study.png')
plt.clf()

# plotting flux distributions TODO determine if we even want this..
# populate list of statepoint file names and extract flux from each statepoint
sp_files = []
flux_dists = []
tallied_powers = []
voxel_volumes = []
source_strengths = []
for n in n_elems:
    sp_files.append('statepoint_' + str(n) + '.h5')

# grab the flux and tallied power for each mesh and append them to the earlier lists
for f in sp_files:
    sp = openmc.StatePoint(f)
    flux_dists.append(sp.get_tally().get_slice(scores=['flux']))
    tallied_powers.append(sp.get_tally(id=2).get_slice(scores=['kappa-fission']))
    # tallied power append
    # compute voxel and source strength, convert flux units via flux = c*flux/V_voxel, apppend or modify flux_dists

