import openmc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# list mesh element sizes that will be used to iterate through all cases
n_elems = [50,100,250,500,1000]
L = 106.47 # equilibrium length from paper
P = 1.0e22 # eV/s

# mesh size vs error
# collect error norm data
# thresholds = ['1e-2','1e-3','1e-4','1e-5','1e-6']
thresholds = ['1e-2']
log_norms = []
for t in thresholds:
    for n in n_elems:
        filename = 'openmc_' + str(n) + '_' + t + '_out.csv'
        df_norms = pd.read_csv(filename)
        last_val  = df_norms.iloc[len(df_norms)-1,-1]
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

# lists to contain data for each mesh size
sp_files = []
flux_mesh_tallies = []
tallied_global_system_powers = []
kappa_fission_mesh_tallies = []
voxel_volumes = []
source_strengths = []

# generate list of statepoints
for n in n_elems:
    sp_files.append('statepoint_' + str(n) + '.h5')
    voxel_volumes.append(L/n)

# populate flux, system tallied power, and local kappa fission rates with
# openmc.Tally objects from each simulation
for f in sp_files:
    sp = openmc.StatePoint(f)
    flux_mesh_tallies.append(sp.get_tally().get_slice(scores=['flux']))
    tallied_global_system_powers.append(sp.get_tally(id=2).get_slice(scores=['kappa-fission']))
    kappa_fission_mesh_tallies.append(sp.get_tally(id=3))

# TODO what is going on with tally 1 vs tally 3 in tallies.out

# compute source strengths and convert flux to proper n/cm^2-s units
for power,voxel_volume in zip(tallied_global_system_powers,voxel_volumes):
    source_strengths.append(P/(voxel_volume*float(power.mean[0])))

# # convert to correct units for flux and kappa fission TODO fix syntax
# for flux,power,ss in zip(flux_mesh_tallies,kappa_fission_mesh_tallies,source_strengths):
#     flux = flux * ss
#     power = power * ss

# grab temperature for 50 case
df_temp_50 = pd.read_csv("openmc_50_out_temp_0004.csv")
temp_50 = df_temp_50.loc[:,"temp"]

# print ratio of temp to flux and see if f is the same or not (i.e. T=f*phi)
# for i in range (len(temp_50)):
#     print(temp_50[i])