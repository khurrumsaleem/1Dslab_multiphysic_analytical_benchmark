import openmc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# list mesh element sizes that will be used to iterate through all cases
n_elems = [50,100,250,500,1000]

# mesh size vs error
# collect error norm data
log_norms = []
for n in n_elems:
    filename = 'openmc_' + str(n) + '_out.csv'
    df = pd.read_csv(filename)
    last_val  = df.iloc[len(df)-1,-1]
    log_norms.append(np.log(last_val))

log_N = np.log(n_elems)
# plot element size vs computed error
plt.plot(log_N,log_norms,'-ro',label=r'$\log(\epsilon)=log(\frac{||T_{a}-T_{x}||_{2}}{||T_{a}||_{2}})$')
plt.yticks(np.log(np.linspace(0.0010,0.0017,6)))
plt.yticks(np.linspace(-6.35,-6.85,6))
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('Log of number of X elements $\log(N)$')
plt.ylabel('Log of Error Norm $log(\epsilon)$')
plt.title('Error from Analytical Solution for Varying Mesh Elements',y=1.05)
plt.legend(bbox_to_anchor=[0.985,0.985])
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

