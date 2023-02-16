import openmc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# list mesh element sizes that will be used to iterate through all cases
n_elems = [50,100,250,500,1000]
thresholds = ['1e-2'] # TODO run simulation for each threshold and add the necessary csv files to this directory, then run w all thresholds
# thresholds = ['1e-2','1e-3','1e-4','1e-5','1e-6']
L = 106.47 # equilibrium length from paper
P = 1.0e22 # eV/s
phi0 = 2.5e14 # n/cm^2-s
q = 1e8 # ev/s
lam = 0.5*(1+np.sqrt(1+(16*q*q*phi0*phi0)/(P*P)))

# dictionary of errors for each triger threshold
error_data = {}
# collect error norm data organized by threshold
for t in thresholds:
    error_data[t] = {}
    for n in n_elems:
        filename = f"openmc_{n}_{t}_out.csv"
        df_norms = pd.read_csv(filename)
        last_val  = df_norms.iloc[len(df_norms)-1,-1]
        error_data[t][n] =  np.log(last_val)

log_N = np.log(n_elems)
error_df = pd.DataFrame.from_dict(error_data)
# plot the error for each threshold against the number of mesh elements
for t in thresholds:
    plt.plot(log_N,error_df.loc[:,t],'-o',label=t)
plt.yticks(np.log(np.linspace(0.0010,0.0017,6))) # TODO may need to adjust once they're all together
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
raw_flux_mesh_tallies = [] # units from openmc
flux_mesh_tallies = [] # proper flux units n/cm^2-s
tallied_global_system_powers = []
raw_kappa_fission_mesh_tallies = [] # units from openmc
kappa_fission_mesh_tallies = [] # proper power units ev/cm^3-s
voxel_volumes = []
source_strengths = []

# generate list of statepoints
for n in n_elems:
    sp_files.append('statepoint_' + str(n) + '.h5')
    voxel_volumes.append(L/n)

# populate flux, system tallied power, and local kappa fission rates with
# openmc.Tally objects from each simulation (openmc units)
for f in sp_files:
    sp = openmc.StatePoint(f)
    raw_flux_mesh_tallies.append(sp.get_tally().get_slice(scores=['flux']))
    tallied_global_system_powers.append(sp.get_tally(id=2).get_slice(scores=['kappa-fission']))
    raw_kappa_fission_mesh_tallies.append(sp.get_tally(id=3))

# compute source strengths and convert flux to proper n/cm^2-s unit
# P [=] ev/s , voxel_volume [=] cm^3, power.mean[0], ev/sp
for power,voxel_volume in zip(tallied_global_system_powers,voxel_volumes):
    source_strengths.append(P/(voxel_volume*float(power.mean[0]))) # units of sp/cm^3-s

# convert to correct units for flux and kappa fission
for flux,power,ss in zip(raw_flux_mesh_tallies,raw_kappa_fission_mesh_tallies,source_strengths):
    flux_mesh_tallies.append(flux * ss)
    kappa_fission_mesh_tallies.append(power * ss)

# PLOTTING T/phi for the 50 and 10000 case
ratio_50 = []
ratio_1000 = []
# generate x coords
xx_50 = np.linspace(-L/2,L/2,n_elems[0])
xx_1000 = np.linspace(-L/2,L/2,n_elems[-1])
# populate data frames
df_temp_50 = pd.read_csv("openmc_50_out_temp_0004.csv")
temp_50 = df_temp_50.loc[:,"temp"]
df_temp_1000 = pd.read_csv("openmc_1000_out_temp_0003.csv")
temp_1000 = df_temp_1000.loc[:,"temp"]
# populate ratio vars
for i in range (len(temp_50)):
    ratio_50.append(float(temp_50[i]/flux_mesh_tallies[0].mean[i]))
for i in range (len(temp_1000)):
    ratio_1000.append(float(temp_1000[i]/flux_mesh_tallies[-1].mean[i]))
plt.plot(xx_50,ratio_50,'-bo',label="50")
plt.plot(xx_1000,ratio_1000,'-ro',label="1000")
plt.xlabel('X coordniate')
plt.ylabel('$T(x)/\phi(x)$')
plt.title('temp to flux ratio for each element')
plt.legend()
plt.savefig('ratio_T_to_flux.png')
plt.clf()

num_to_analy = []
# flux numerical to analytical ratio
phi = phi0*np.sqrt(1- ((lam - 1)*P*P*np.multiply(xx_50,xx_50))/(L*L*q*q*phi0*phi0))
for i in range(50):
    num_to_analy.append(float(flux_mesh_tallies[0].mean[i] / phi[i]))
plt.plot(xx_50,num_to_analy,'-go',)
plt.xlabel('X coordniate')
plt.ylabel('numerical to analytical ratio $\phi(x)$')
plt.title('Ratio numerical to analytical flux')
plt.savefig('num_to_analytical_flux.png')
plt.clf()

# fission_source = []
# # flux error ratio
# for i in range(50):
#     # print(phi[i])
#     fission_source.append(float(kappa_fission_mesh_tallies[0].mean[i]))
# # plt.savefig('heat_source.png')