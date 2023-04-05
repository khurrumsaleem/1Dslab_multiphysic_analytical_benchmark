import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import pandas as pd

# list mesh element sizes that will be used to iterate through all cases
# n_elems = [50,100,250,500,1000]
n_elems = [50,100,250,500,1000]
log_N = np.log(n_elems)
L = 106.47 # equilibrium length from paper
P = 1.0e22 # eV/s
T0 = 293 # K
phi0 = 2.5e14 # n/cm^2-s
q = 1e8 # ev/s
k0= 1.25e19 # eV/(s-cm-K^2) k(T) = k0 T(x)
lam = 0.5*(1+np.sqrt(1+(16*q*q*phi0*phi0)/(P*P)))
Sig_t0 = np.sqrt(P/((lam-1)*k0*L))/(T0)

log_temp_error_norms = [] # list to house the error for each mesh size
for n in n_elems:
    filename = f"openmc_{n}_out.csv"
    df_norms = pd.read_csv(filename)
    last_val  = df_norms.iloc[len(df_norms)-1,-1]
    log_temp_error_norms.append(np.log(last_val))

# plot the error for each threshold against the number of mesh elements
plt.plot(log_N,log_temp_error_norms,'-ro',label=r'$\epsilon_{T}=\frac{||T_{a} - T_{x} ||_{2}}{|| T_{a} ||_{2}}$')
plt.yticks(np.linspace(-6,-10,6))
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel(r"Log of number of X elements $\log(N)$",fontsize=16)
plt.ylabel(r"Log of Error Norm $log(\epsilon_{T})$",fontsize=16)
# plt.title("Temperature Error Norm vs Mesh Size \n 200 Picard Iterations",y=1)
plt.legend(bbox_to_anchor=[0.985,0.985],fontsize=14)
plt.grid()
plt.savefig("temp_error_norms.png",bbox_inches='tight')
plt.clf()

# linear polynomial fit to get slope of line of best fit temp
pt = np.polyfit(log_N,log_temp_error_norms,1)
print("polyfit:" , pt)

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

print(eigs)

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
cardinal_fluxes = {}
for n in n_elems:
    xx[n] = np.linspace(-L/2,L/2,n)
    temps[n] = pd.read_csv(f"openmc_{n}_temp.csv").loc[:,"temp"]
    cardinal_fluxes[n] = pd.read_csv(f"openmc_{n}_flux.csv").loc[:,"flux"]
    if(n==50):
        plt.plot(xx[n],temps[n],'-ro')
    else:
        plt.plot(xx[n],temps[n])
    plt.xticks([-60,-40,-20,0,20,40,60])
    plt.xlabel("X coordinate [cm]",fontsize=16)
    plt.ylabel("Temperature [K]",fontsize=16)
    # plt.title(f"Temperature for {n} Mesh Elements")
    plt.grid()
    plt.savefig(f"temp_{n}.png",bbox_inches='tight')
    plt.clf()

# dictioinaries for flux results and comparisons to temp and analytical solution key is number of elements
ratios_T_to_phi = {} # list of ratios of temp to flux for each element
analytical_phi = {}
analytical_T = {}
ratios_flux_num_to_analy = {}
ratios_temp_num_to_analy = {}
# compute T/phi for each case
for n in n_elems:
    ratios_T_to_phi[n] = [] # empty list to store ratios T to phi
    ratios_temp_num_to_analy[n] = []
    ratios_flux_num_to_analy[n] = [] # empty list to store ratios numerical to analytical flux
    analytical_phi[n] = phi0*np.sqrt(1- ((lam - 1)*P*P*np.multiply(xx[n],xx[n]))/(L*L*q*q*phi0*phi0))
    analytical_T[n] = Sig_t0*T0*np.sqrt( (q*q*L*L*phi0*phi0)/(P*P) - (lam-1)*np.multiply(xx[n],xx[n]))
    for i in range(n):
        # ratios_flux_num_to_analy[n].append(float(flux_mesh_tallies[n].mean[i]/analytical_phi[n][i])) # openmc flux tally way
        ratios_flux_num_to_analy[n].append(float(cardinal_fluxes[n][i]/analytical_phi[n][i]))
        # ratios_T_to_phi[n].append(float(temps[n][i]/flux_mesh_tallies[n].mean[i]))
        ratios_temp_num_to_analy[n].append(float(temps[n][i]/analytical_T[n][i]))
    # numerical to analytical ratio
    plt.plot(xx[n],ratios_flux_num_to_analy[n],'-go',label=f"{n} x-elem")
    plt.xticks([-60,-40,-20,0,20,40,60])
    plt.xlabel("X Coordinate [cm]")
    plt.ylabel(r"Numerical $\phi(x)$ to Analytical $\phi(x)$ Ratio")
    plt.title(f"Ratio Numerical to Analytical Flux {n} Mesh Elements")
    plt.legend()
    plt.grid()
    plt.savefig(f"num_to_analytical_flux_{n}.png")
    plt.clf()

# plot all temp error ratios on one plot
for n in n_elems:
    plt.plot(xx[n],ratios_temp_num_to_analy[n],label=f"{n} x-elem")
plt.xticks([-60,-40,-20,0,20,40,60])
plt.yticks([0.9995,1.000,1.001,1.002,1.003,1.004,1.0045])
plt.xlabel("X Coordinate [cm]",fontsize=16)
plt.ylabel(r"Numerical $T(x)$ to Analytical $T(x)$ Ratio",fontsize=16)
# plt.title(f"Ratio Numerical to Analytical Temp All Meshes. 200 Picard Iterations")
plt.legend(ncol=2)
plt.grid()
plt.savefig("temp_num_to_analy_ratios.png",bbox_inches='tight')
plt.clf()

# compute error norm
flux_means = {}
flux_std_dev = {}
flux_deviations = {}
flux_dev_norms = {}
analytical_norms = {}
flux_error_norms = {}
sigma_r = {}
# flux error ratio
for n in n_elems:
    flux_means[n] = []
    flux_std_dev[n] = []
    flux_deviations[n] = []
    sigma_r[n] = []
    # compute mean, stdev, and analytical - numerical
    for i in range(n):
        # openmc flux ways
        # flux_means[n].append(float(flux_mesh_tallies[n].mean[i]))
        # flux_std_dev[n].append(float(flux_mesh_tallies[n].std_dev[i]))
        # flux_deviations[n].append(analytical_phi[n][i]-float(flux_mesh_tallies[n].mean[i]))
        flux_deviations[n].append(analytical_phi[n][i]-float(cardinal_fluxes[n][i]))
    # compute uncertainty in ratio
    # sigma_r^2 = flux_std_dev^2 / analytical_phi^2
    # for i in range(n):
    #     sigma_r[n].append(float( flux_std_dev[n][i]/analytical_phi[n][i]) )
    # take the norm of the deviation and analytical solution
    flux_dev_norms[n] = np.linalg.norm(flux_deviations[n],ord=2)
    analytical_norms[n] = np.linalg.norm(analytical_phi[n],ord=2)
    # ratio is error norm
    flux_error_norms[n] = float(flux_dev_norms[n]/analytical_norms[n])
    # plot flux with error bars
    # plt.plot(xx[n],flux_means[n],'-ko',label=f"{n} x-elem")
    # plt.errorbar(xx[n],flux_means[n],yerr=flux_std_dev[n],marker = '|',fmt='none',elinewidth=1,capsize=3,capthick=1)
    # plt.xlabel('X Coordniate [cm]')
    # plt.ylabel("Flux [n/cm^2-s]")
    # plt.title("Tallied Flux (Units Converted Using SS)")
    # plt.grid()
    # plt.savefig(f"flux_{n}.png")
    # plt.clf()

# plot all flux error ratios on one plot
for n in n_elems:
    plt.plot(xx[n],ratios_flux_num_to_analy[n],label=f"{n} x-elem")
plt.xticks([-60,-40,-20,0,20,40,60])
plt.yticks([0.9995,1.000,1.001,1.002,1.003,1.004,1.0045])
plt.xlabel("X Coordinate [cm]",fontsize=16)
plt.ylabel(r"Numerical $\phi(x)$ to Analytical $\phi(x)$ Ratio",fontsize=16)
# plt.title("Ratio Numerical to Analytical Flux All Meshes. 200 Picard Iterations")
plt.legend(ncol=2)
plt.grid()
plt.savefig("flux_num_to_analy_ratios.png",bbox_inches='tight')
plt.clf()

# # plot individual C/E with error bars
# for n in n_elems:
#     plt.plot(xx[n],ratios_flux_num_to_analy[n],label=f"{n} x-elem")
#     plt.errorbar(xx[n],ratios_flux_num_to_analy[n],yerr=sigma_r[n],marker = '|',fmt='none',elinewidth=0.25,capsize=3,capthick=1)
#     plt.xlabel("X Coordinate [cm]")
#     plt.ylabel(r"Numerical $\phi(x)$ to Analytical $\phi(x)$ Ratio")
#     plt.axis([-54, 54, 0.998, 1.006])
#     plt.title(f"Computed to Expected Ratio with Error Bars. \n {n} Elements. 200 Picard Iterations")
#     plt.grid()
#     plt.savefig(f"flux_num_to_analy_ratio_w_error_bars_{n}_elem.png")
#     plt.clf()

# # compute and plot flux relative errors all on one plot
# rel_err = {}
# for n in n_elems:
#     rel_err[n] = []
#     for i in range(n):
#         rel_err[n].append(np.divide(flux_std_dev[n][i],flux_means[n][i]))
#     plt.plot(xx[n],rel_err[n],label=f"{n}")

# plt.xlabel("Position [cm]")
# plt.ylabel("Relative Error")
# plt.title("Relative Errors for Each Mesh Size (200 Picard)")
# plt.grid()
# plt.legend()
# plt.savefig("relative_errors.png",bbox_inches='tight')
# plt.clf()


log_flux_error_norms = [np.log(flux_error_norms[n]) for n in n_elems]
print(log_flux_error_norms)
# linear polynomial fit to get slope of line of best fit
pf = np.polyfit(log_N,log_flux_error_norms,1)
print("polyfit:" , pf)

plt.plot(log_N,log_flux_error_norms,'-o',label=r'$\epsilon_{\phi}=\frac{||\phi_{a} - \phi_{x} ||_{2}}{|| \phi_{a} ||_{2}}$')
plt.xlabel(r"Log of number of X elements $\log(N)$",fontsize=16)
plt.ylabel(r"Log of Error Norm $log(\epsilon_{\phi})$",fontsize=16)
# plt.title("Flux Error Norm vs Mesh Size \n 200 Picard Iterations with Relaxation")
plt.yticks(np.linspace(-6,-10,6))
plt.grid()
plt.legend(bbox_to_anchor=[0.985,0.985],fontsize=14)
plt.savefig("flux_error_norms.png",bbox_inches='tight')
plt.clf()
