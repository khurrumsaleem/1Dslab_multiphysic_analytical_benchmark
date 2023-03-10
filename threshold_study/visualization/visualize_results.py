import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pandas as pd

# list mesh element sizes that will be used to iterate through all cases
n_elems = [50,100,250,500,1000]
log_N = np.log(n_elems)
L = 106.47 # equilibrium length from paper
P = 1.0e22 # eV/s
phi0 = 2.5e14 # n/cm^2-s
q = 1e8 # ev/s
lam = 0.5*(1+np.sqrt(1+(16*q*q*phi0*phi0)/(P*P)))

# dictionary of errors for each triger threshold
error_data = []
# collect error norm data organized by threshold
for n in n_elems:
    filename = f"openmc_{n}_1e-2_out.csv"
    df_norms = pd.read_csv(filename)
    last_val  = df_norms.iloc[len(df_norms)-1,-1]
    error_data.append(np.log(last_val))

# plot the error for each threshold against the number of mesh elements
plt.plot(log_N,error_data,'-ro',label=r'$\epsilon=\frac{||T_{a} - T_{x} ||_{2}}{|| T_{a} ||_{2}}$')
plt.yticks(np.linspace(-6.5,-10,6))
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel('Log of number of X elements $\log(N)$') 
plt.ylabel('Log of Error Norm $log(\epsilon)$')
plt.title('Error from Analytical Solution for Varying Mesh Elements \n Steady State Tolerance 1e-5',y=1)
plt.legend(bbox_to_anchor=[0.985,0.985],fontsize=14)
plt.grid()
plt.savefig('error_study.png')
plt.clf()

# lists to contain data for each mesh size
sp_files = {}
raw_flux_mesh_tallies = [] # units from openmc
flux_mesh_tallies = [] # proper flux units n/cm^2-s
tallied_global_system_powers = []
raw_kappa_fission_mesh_tallies = [] # units from openmc
kappa_fission_mesh_tallies = [] # proper power units ev/cm^3-s
voxel_volumes = []
source_strengths = []
eigs = [] # eigenvalues to fix normalization of flux

# generate list of statepoints
for n in n_elems:
    sp_files.append('statepoint_' + str(n) + '.h5')
    voxel_volumes.append(L/n)

# populate flux, system tallied power, and local kappa fission rates with
# openmc.Tally objects from each simulation (openmc units)
for f in sp_files:
    sp = openmc.StatePoint(f)
    raw_flux_mesh_tallies.append(sp.get_tally(id=1).get_slice(scores=['flux']))
    tallied_global_system_powers.append(sp.get_tally(id=3).get_slice(scores=['kappa-fission']))
    raw_kappa_fission_mesh_tallies.append(sp.get_tally(id=1).get_slice(scores=['kappa-fission']))
    eigs.append(sp.keff.n)


nu_fission_rate = float(openmc.StatePoint(sp_files[0]).get_tally(id=2).get_slice(scores=['nu-fission']).mean[0])
fission_rate = float(openmc.StatePoint(sp_files[0]).get_tally(id=2).get_slice(scores=['fission']).mean[0])
# print(tallied_global_system_powers[0].mean[0]/fission_rate) # fission rate computed is very close to provided value
# print(nu_fission_rate/fission_rate) # nu is as expected

# compute source strengths and convert flux to proper n/cm^2-s unit
# P [=] ev/s , voxel_volume [=] cm^3, power.mean[0], ev/sp
for power,voxel_volume,keff in zip(tallied_global_system_powers,voxel_volumes,eigs):
    source_strengths.append(P*nu_fission_rate/(keff*voxel_volume*float(power.mean[0]))) # units of sp/cm^3-s

# convert to correct units for flux and kappa fission
for flux,power,ss in zip(raw_flux_mesh_tallies,raw_kappa_fission_mesh_tallies,source_strengths):
    flux_mesh_tallies.append(flux * ss)
    kappa_fission_mesh_tallies.append(power * ss)

# PLOTTING T/phi for the 50 and 10000 case
# dictionaries to house data for plotting
temps = {} # key is number of elements, value is data frame read in from csv
ratios = {} # key is number of elements, value is numpary array of ratio of temp to flux for each element
xx = {} # key is number of elements, value is numpy array of mesh center points
# generate x coords
for n in n_elems:
    xx[n] = np.linspace(-L/2,L/2,n)
    temps[n] = pd.read_csv(f"openmc_{n}_temp.csv")
# populate dictionaries

for n in n_elems:
    plt.plot(xx[n],temps[n].loc[:,"temp"],'-ro')
    plt.xlabel("X coordinate [cm]")
    plt.ylabel("Temperature [K]")
    plt.title(f"Temperature for {n} Mesh Elements")
    plt.savefig(f'temp_{n}.png')
    plt.clf()

# populate ratio vars
# for i in range (len(temp_50)):
#     ratio_50.append(float(temp_50[i]/flux_mesh_tallies[0].mean[i]))
# plt.plot(xx_50,ratio_50,'-bo',label="50 x-elem")
# plt.xlabel('X coordniate')
# plt.ylabel('$T(x)/\phi(x)$')
# plt.title('temp to flux ratio for each element')
# plt.legend()
# plt.savefig('ratio_T_to_flux.png')
# plt.clf()

# num_to_analy = []
# # flux numerical to analytical ratio
# phi = phi0*np.sqrt(1- ((lam - 1)*P*P*np.multiply(xx_50,xx_50))/(L*L*q*q*phi0*phi0))
# for i in range(50):
#     num_to_analy.append(float(flux_mesh_tallies[0].mean[i] / phi[i]))
# plt.plot(xx_50,num_to_analy,'-go',)
# plt.xlabel('X coordniate')
# plt.ylabel('numerical $\phi(x)$ to analytical $\phi(x)$ ratio')
# plt.title('Ratio numerical to analytical flux')
# plt.savefig('num_to_analytical_flux.png')
# plt.clf()


# plt.plot(xx_50,flux_mesh_tallies[0].mean.flat,'-ko')
# plt.title("numericall flux vs position")
# plt.xlabel('X coordniate')
# plt.ylabel('Numerical $\phi(x)$')
# plt.savefig("numerical_flux.png")
# plt.clf()

# kappa_fission_means = []
# kappa_fission_std_dev = []
# # flux error ratio
# for i in range(50):
#     kappa_fission_means.append(float(kappa_fission_mesh_tallies[0].mean[i]))
#     kappa_fission_std_dev.append(float(kappa_fission_mesh_tallies[0].std_dev[i]))
# plt.plot(xx_50,kappa_fission_means,'-ko')
# plt.errorbar(xx_50,kappa_fission_means,yerr=kappa_fission_std_dev,marker = '|',fmt='none',elinewidth=1,capsize=3,capthick=1)
# plt.xlabel('X coordniate')
# plt.ylabel('heat source [ev/cm^3-s]')
# plt.title('Tallied Kappa Fission (Units Converted Using SS)')
# plt.savefig('heat_source.png')
# plt.clf()


# iterations to steady state plot
m = [50,100,250,500,1000]
its = [31,46,72,101,143]

fig, ax = plt.subplots()
plt.plot(m,its,'-bo')
ax.loglog()
ax.set_xticks(m)
ax.set_yticks([25,50,75,100,125,150])
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.minorticks_off()
ax.grid(which='major')
plt.xlabel("Mesh Elements")
plt.ylabel("Iterations to Steady State")
plt.title("Comparing the Iterations Required for MOOSE \n to Detect Steady State (Tolerance = 1e-5)")
plt.savefig("its_to_ss.png")