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

error_data = [] # list to house the error for each mesh size
for n in n_elems:
    filename = f"openmc_{n}_1e-2_out.csv"
    df_norms = pd.read_csv(filename)
    last_val  = df_norms.iloc[len(df_norms)-1,-1]
    error_data.append(np.log(last_val))

# plot the error for each threshold against the number of mesh elements
plt.plot(log_N,error_data,'-ro',label=r'$\epsilon=\frac{||T_{a} - T_{x} ||_{2}}{|| T_{a} ||_{2}}$')
plt.yticks(np.linspace(-6.5,-10,6))
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel("Log of number of X elements $\log(N)$")
plt.ylabel("Log of Error Norm $log(\epsilon)$")
plt.title("Error from Analytical Solution for Varying Mesh Elements \n Steady State Tolerance 1e-5",y=1)
plt.legend(bbox_to_anchor=[0.985,0.985],fontsize=14)
plt.grid()
plt.savefig("error_study.png")
plt.clf()

# linear polynomial fit to get slope of line of best fit
p = np.polyfit(log_N,error_data,1)
# print("polyfit:" , p)

# dictionaries of statepoint filenames and voxel volumes
# key in each case is n (the number of mesh elements).
# convention applies to rest of dictionaries unless otherwise noted
sp_filenames = {}
voxel_volumes = {}

# generate filenames and compute voxel volumes
for n in n_elems:
    sp_filenames[n] = f'statepoint_{n}.h5'
    voxel_volumes[n] = L/n # implied division by 1 x 1 (y by z) cross sectional are

# dictionaries of raw tally data for each mesh size
raw_flux_mesh_tallies = {} # units of n-cm/sp
tallied_global_system_powers = {} # units of eV-cm/sp
raw_kappa_fission_mesh_tallies = {} # units of eV-cm/sp
nu_fission_rates = {} # units of particles/sp - INTEGRAL QUANTITY, store foat
fission_rates = {} # units of fissions/sp - INTEGRAL QUANTITY, store foat
eigs = {} # k effective - one data point per mesh size, store float

# populate flux, system tallied power, and local kappa fission rates with
# openmc.Tally objects from each simulation (openmc units)
for n in n_elems:
    sp = openmc.StatePoint(sp_filenames[n])
    raw_flux_mesh_tallies[n] = sp.get_tally(id=1).get_slice(scores=['flux'])
    tallied_global_system_powers[n] = sp.get_tally(id=3).get_slice(scores=['kappa-fission'])
    raw_kappa_fission_mesh_tallies[n] = sp.get_tally(id=1).get_slice(scores=['kappa-fission'])
    nu_fission_rates[n] = float(openmc.StatePoint(sp_filenames[n]).get_tally(id=2).get_slice(scores=['nu-fission']).mean[0])
    fission_rates[n] = float(openmc.StatePoint(sp_filenames[n]).get_tally(id=2).get_slice(scores=['fission']).mean[0])
    eigs[n] = [float(sp.keff.n),float(sp.keff.s)]

# print(eigs)

source_strengths = {}
# compute source strengths and convert flux to proper n/cm^2-s unit
# P [=] ev/s, nu_fission_rates [=] n/sp, voxel_volume [=] cm^3, eigs [=] n/sp, voxel_volumes [=] cm^3, tallied_global_system_powers.mean[0], ev/sp
for n in n_elems:
    source_strengths[n] = P*nu_fission_rates[n]/(eigs[n][0]*voxel_volumes[n]*float(tallied_global_system_powers[n].mean[0])) # units of sp/cm^3-s

flux_mesh_tallies = {} # proper flux units n/cm^2-s
kappa_fission_mesh_tallies = {} # proper power units ev/cm^3-s
# convert to correct units for flux and kappa fission
for n in n_elems:
    flux_mesh_tallies[n] = raw_flux_mesh_tallies[n] * source_strengths[n]
    kappa_fission_mesh_tallies[n] = raw_kappa_fission_mesh_tallies[n] * source_strengths[n]

# dictionaries to house data for plotting
xx = {} # key is number of elements, value is numpy array of mesh center points
temps = {} # key is number of elements, value is data frame read in from csv

for n in n_elems:
    xx[n] = np.linspace(-L/2,L/2,n)
    temps[n] = pd.read_csv(f"openmc_{n}_temp.csv").loc[:,"temp"]

for n in n_elems:
    plt.plot(xx[n],temps[n],'-ro')
    plt.xlabel("X coordinate [cm]")
    plt.ylabel("Temperature [K]")
    plt.title(f"Temperature for {n} Mesh Elements")
    plt.savefig(f"temp_{n}.png")
    plt.clf()

# dictioinaries for flux results and comparisons to temp and analytical solution key is number of elements
ratios_T_to_phi = {} # list of ratios of temp to flux for each element
analytical_phi = {}
ratios_flux_num_to_analy = {}
# compute T/phi for each case
for n in n_elems:
    ratios_T_to_phi[n] = [] # empty list to store ratios T to phi
    ratios_flux_num_to_analy[n] = [] # empty list to store ratios numerical to analytical flux
    analytical_phi[n] = phi0*np.sqrt(1- ((lam - 1)*P*P*np.multiply(xx[n],xx[n]))/(L*L*q*q*phi0*phi0))
    for i in range(n):
        ratios_flux_num_to_analy[n].append(float(flux_mesh_tallies[n].mean[i]/analytical_phi[n][i]))
        ratios_T_to_phi[n].append(float(temps[n][i]/flux_mesh_tallies[n].mean[i]))
    # numerical to analytical ratio
    plt.plot(xx[n],ratios_flux_num_to_analy[n],'-go',label=f"{n} x-elem")
    plt.xlabel("X Coordinate [cm]")
    plt.ylabel(r"Numerical $\phi(x)$ to Analytical $\phi(x)$ Ratio")
    plt.title(f"Ratio Numerical to Analytical Flux {n} Mesh Elements")
    plt.legend()
    plt.savefig(f"num_to_analytical_flux_{n}.png")
    plt.clf()
    # T to phi ratio
    plt.plot(xx[n],ratios_T_to_phi[n],'-bo',label=f"{n} x-elem")
    plt.xlabel("X Coordinate [cm]")
    plt.ylabel(r"Ratio $\frac{T(x)}{\phi(x)}$ ")
    plt.title(f"Temp to Flux Ratio for {n} Mesh Elements")
    plt.legend()
    plt.savefig(f"ratio_T_to_flux_{n}.png")
    plt.clf()



# plot all error ratios on one plot
for n in n_elems:
    plt.plot(xx[n],ratios_flux_num_to_analy[n],label=f"{n} x-elem")
plt.xlabel("X Coordinate [cm]")
plt.ylabel(r"Numerical $\phi(x)$ to Analytical $\phi(x)$ Ratio")
plt.title(f"Ratio Numerical to Analytical Flux All Meshes")
plt.legend()
plt.savefig("flux_num_to_analy_ratios.png")
plt.clf()

# compute error norm
flux_means = {}
flux_std_dev = {}
flux_deviations = {}
flux_dev_norms = {}
analytical_norms = {}
error_norm = {}
# flux error ratio
for n in n_elems:
    flux_means[n] = []
    flux_std_dev[n] = []
    flux_deviations[n] = []
    # compute mean, stdev, and analytical - numerical
    for i in range(n):
        flux_means[n].append(float(flux_mesh_tallies[n].mean[i]))
        flux_std_dev[n].append(float(flux_mesh_tallies[n].std_dev[i]))
        flux_deviations[n].append(analytical_phi[n][i]-float(flux_mesh_tallies[n].mean[i]))
    # take the norm of the deviation and analytical solution
    flux_dev_norms[n] = np.linalg.norm(flux_deviations[n],ord=2)
    analytical_norms[n] = np.linalg.norm(analytical_phi[n],ord=2)
    # ratio is error norm
    error_norm[n] = float(flux_dev_norms[n]/analytical_norms[n])
    # plt.plot(xx[n],flux_means[n],'-ko',label=f"{n} x-elem")
    # plt.errorbar(xx[n],flux_means[n],yerr=flux_std_dev[n],marker = '|',fmt='none',elinewidth=1,capsize=3,capthick=1)
    # plt.xlabel('X Coordniate [cm]')
    # plt.ylabel("Flux [n/cm^2-s]")
    # plt.title("Tallied Flux (Units Converted Using SS)")
    # plt.grid()
    # plt.savefig(f"flux_{n}.png")
    plt.clf()

for n in n_elems:
    plt.plot([n],[error_norm[n]],'-o',label=f"{n}")
plt.xlabel("number of elements")
plt.ylabel("error norm")
plt.xticks([n for n in n_elems])
plt.grid()
plt.savefig("flux_error_norms.png")
plt.clf()


kappa_fission_means = {}
kappa_fission_std_dev = {}
# flux error ratio
for n in n_elems:
    kappa_fission_means[n] = []
    kappa_fission_std_dev[n] = []
    for i in range(n):
        kappa_fission_means[n].append(float(kappa_fission_mesh_tallies[n].mean[i]))
        kappa_fission_std_dev[n].append(float(kappa_fission_mesh_tallies[n].std_dev[i]))
    plt.plot(xx[n],kappa_fission_means[n],'-ko',label=f"{n} x-elem")
    plt.xlabel('X Coordniate [cm]')
    plt.ylabel("Heat Source [ev/cm^3-s]")
    plt.title("Tallied Kappa Fission (Units Converted Using SS)")
    plt.savefig(f"heat_source_{n}.png")
    plt.clf()


# iterations to steady state plot
m = [50,100,250,500,1000]
its = [33,48,72,102,144]

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