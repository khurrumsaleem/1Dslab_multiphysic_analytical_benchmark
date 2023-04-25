import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker as mtick
import pandas as pd

# mesh element sizes that will be used to iterate through all cases
n_elems = np.array([50,100,250,500,1000])
log_N = np.log10(n_elems)
L = 106.47 # equilibrium length from paper
P = 1.0e22 # eV/s
T0 = 293 # K
phi0 = 2.5e14 # n/cm^2-s
q = 1e8 # ev/s
k0= 1.25e19 # eV/(s-cm-K^2) k(T) = k0 T(x)
lam = 0.5*(1+np.sqrt(1+(16*q*q*phi0*phi0)/(P*P)))
Sig_t0 = np.sqrt(P/((lam-1)*k0*L))/(T0)

# read in x coordinates from CSV
xx = {} # key is number of elements, value is numpy array of mesh center points
ones = {} # need a ones vector of every mesh size for the C/E
for n in n_elems:
    xx[n] = pd.read_csv(f"openmc_{n}_temp.csv").loc[:,"x"]
    ones[n] = np.ones(len(xx[n]))

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
    raw_kappa_fission_mesh_tallies[n] = sp.get_tally(id=1).get_slice(scores=['kappa-fission'])
    tallied_global_system_powers[n] = sp.get_tally(id=2).get_slice(scores=['kappa-fission'])
    eigs[n] = [float(sp.keff.n),float(sp.keff.s)]

print(eigs)

source_strengths = {}
# compute source strengths and convert flux to proper n/cm^2-s unit
# P [=] ev/s, nu_fission_rates [=] n/sp, voxel_volume [=] cm^3, eigs [=] n/sp, voxel_volumes [=] cm^3, tallied_global_system_powers.mean[0], ev/sp
for n in n_elems:
    source_strengths[n] = P/(voxel_volumes[n]*float(tallied_global_system_powers[n].mean[0])) # units of sp/cm^3-s


flux_mesh_tallies = {} # proper flux units n/cm^2-s
kappa_fission_mesh_tallies = {} # proper power units ev/cm^3-s
# convert to correct units for flux and kappa fission
for n in n_elems:
    flux_mesh_tallies[n] = raw_flux_mesh_tallies[n] * source_strengths[n]
    kappa_fission_mesh_tallies[n] = raw_kappa_fission_mesh_tallies[n] * source_strengths[n]


# dictioinaries for flux results and comparisons to temp and analytical solution key is number of elements
analytical_phi = {}
ratios_flux_num_to_analy = {}
# compute T/phi for each case
for n in n_elems:
    ratios_flux_num_to_analy[n] = np.zeros(n)
    analytical_phi[n] = phi0*np.sqrt(1- ((lam - 1)*P*P*np.multiply(xx[n],xx[n]))/(L*L*q*q*phi0*phi0))
    for i in range(n):
        ratios_flux_num_to_analy[n][i] = float(flux_mesh_tallies[n].mean[i]/analytical_phi[n][i]) # openmc flux tally way

# compute error norm and uncertainty in ratio
flux_means = {}
flux_std_devs = {}
KF_means = {}
KF_std_devs = {}
flux_deviations = {}
flux_dev_norms = {}
analytical_norms = {}
flux_error_norms = {}
sigma_r = {}
two_sigma_r = {}
# flux error ratio
for n in n_elems:
    flux_means[n] = np.zeros(n)
    flux_std_devs[n] = np.zeros(n)
    flux_deviations[n] = np.zeros(n)
    KF_means[n] = np.zeros(n)
    KF_std_devs[n] = np.zeros(n)
    sigma_r[n] = np.zeros(n)
    # analytical - numerical
    for i in range(n):
        flux_means[n][i] = float(flux_mesh_tallies[n].mean[i])
        flux_deviations[n][i] = analytical_phi[n][i]-float(flux_means[n][i])
        flux_std_devs[n][i] = float(flux_mesh_tallies[n].std_dev[i])
        KF_means[n][i] = float(kappa_fission_mesh_tallies[n].mean[i])
        KF_std_devs[n][i] = float(kappa_fission_mesh_tallies[n].std_dev[i])
    # take the norm of the deviation and analytical solution
    # compute uncertainty in ratio
    # sigma_r^2 = (flux_std_dev/analytical_phi)^2 + ([flux_mean*voxel_volume*kf_std_dev]/[analytical_phi*P])^2
    for i in range(n):
        sigma_r[n][i] =np.sqrt( ( ( flux_std_devs[n][i]*flux_std_devs[n][i] ) / ( analytical_phi[n][i]*analytical_phi[n][i] ) ) +\
        ( flux_means[n][i]*flux_means[n][i]*voxel_volumes[n]*voxel_volumes[n]*KF_std_devs[n][i]*KF_std_devs[n][i] )/( analytical_phi[n][i]*analytical_phi[n][i]*P*P ) )
    flux_dev_norms[n] = np.linalg.norm(flux_deviations[n],ord=2)
    analytical_norms[n] = np.linalg.norm(analytical_phi[n],ord=2)
    # ratio is error norm
    flux_error_norms[n] = float(flux_dev_norms[n]/analytical_norms[n])
    two_sigma_r[n] = 2*sigma_r[n]


# plot numerical flux for 50 mesh elements
plt.plot(xx[50],flux_means[50],'-o')
plt.errorbar(xx[50],flux_means[50],yerr=flux_std_devs[50],marker = '|',fmt='none',elinewidth=1,capsize=4,capthick=2)
plt.xlabel("X Coordinate [cm]",fontsize=16)
plt.xticks([-60,-40,-20,0,20,40,60])
plt.ylabel(r"Flux [n/$cm^{2}$-s]",fontsize=16)
plt.grid()
plt.rcParams['figure.figsize'] = [4, 4]
plt.savefig("flux_50.png", bbox_inches="tight")
plt.clf()


# plot all flux C/E
for n in n_elems:
    plt.plot(xx[n],ratios_flux_num_to_analy[n],label=f"{n} x-elem")
    # plt.errorbar(xx[n],ratios_flux_num_to_analy[n],yerr=two_sigma_r[n],marker = '|',fmt='none',elinewidth=0.25,capsize=3,capthick=1)
plt.plot(xx[n],ones[n],'k',label="exact")
plt.xticks([-60,-40,-20,0,20,40,60])
plt.yticks([1-0.6e-3,1-0.4e-3,1-0.2e-3,1,1+0.2e-3,1+0.4e-3,1+0.6e-3])
plt.xlabel("X Coordinate [cm]",fontsize=16)
plt.ylabel(r"Flux C/E",fontsize=16)
plt.gca().yaxis.tick_right()
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4f'))
plt.legend(ncol=2)
plt.grid()
plt.savefig("flux_num_to_analy_ratios.png",bbox_inches='tight')
plt.clf()

# plot C/E = 1 with error bounds on same plot as C/E
for n in n_elems:
    plt.plot(xx[n],ratios_flux_num_to_analy[n],label=f"{n} x-elem")
    plt.plot(xx[n],ones[n],'k',label="exact")
    plt.plot(xx[n],ones[n]+two_sigma_r[n],'-m')
    plt.plot(xx[n],ones[n]-two_sigma_r[n],'-m')
    plt.gca().fill_between(xx[n],ones[n]+two_sigma_r[n],ones[n]-two_sigma_r[n],color='lightgray')
    plt.xticks([-60,-40,-20,0,20,40,60])
    plt.yticks([1-1e-3,1-0.5e-3,1,1+0.5e-3,1+1e-3])
    plt.xlabel("X Coordinate [cm]",fontsize=16)
    plt.ylabel(r"Flux C/E",fontsize=16)
    plt.gca().yaxis.tick_right()
    plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4f'))
    plt.legend()
    plt.grid()
    plt.savefig(f"{n}_flux_CE_error_bars.png",bbox_inches='tight')
    plt.clf()


log_flux_error_norms = np.array([np.log10(flux_error_norms[n]) for n in n_elems])

# results_T = stats.linregress(log_N,log_temp_error_norms)
# print("slope = ", results_T.slope, " r-squared = ",results_T.rvalue*results_T.rvalue)

results_flux = stats.linregress(log_N,log_flux_error_norms)
print("slope = ", results_flux.slope, " r-squared = ",results_flux.rvalue*results_flux.rvalue)


# plot mesh elements vs flux error norms
plt.plot(log_N,log_flux_error_norms,'-o',label=r'$\epsilon_{\phi}=\frac{||\phi_{a} - \phi_{x} ||_{2}}{|| \phi_{a} ||_{2}}$')
plt.xlabel(r"Log of number of X elements $\log(N)$",fontsize=16)
plt.ylabel(r"Log of Error Norm $log(\epsilon_{\phi})$",fontsize=16)
# plt.title("Flux Error Norm vs Mesh Size \n 200 Picard Iterations with Relaxation")
plt.xticks(log_N)
plt.yticks([-3.6,-3.7,-3.8,-3.9])
plt.grid()
plt.legend(bbox_to_anchor=[0.985,0.215],fontsize=14)
plt.savefig("flux_error_norms.png",bbox_inches='tight')
plt.clf()
