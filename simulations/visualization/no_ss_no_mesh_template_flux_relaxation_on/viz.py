# a tidier visualization script that takes Cardinal outputs
# and plots C/E and error norms for the temperature and flux

import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker as mtick
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

#### COMPUTED TO EXPECTED TEMP AND FLUX ####

# dictionaries to house data
# keys are the number of elements and values are pandas data frame
xx = {} # x-centroids
temps = {} # temperature [K]
cardinal_fluxes = {} # flux from Cardinal [n/cm^2-s]
for n in n_elems:
    xx[n] = pd.read_csv(f"openmc_{n}_temp.csv").loc[:,"x"]
    temps[n] = pd.read_csv(f"openmc_{n}_temp.csv").loc[:,"temp"]
    cardinal_fluxes[n] = pd.read_csv(f"openmc_{n}_flux.csv").loc[:,"flux"]


# plot 50 mesh element temperature
plt.plot(xx[50],temps[50],'-ro')
plt.xticks([-60,-40,-20,0,20,40,60])
plt.xlabel("X coordinate [cm]",fontsize=16)
plt.ylabel("Temperature [K]",fontsize=16)
# plt.title(f"Temperature for {50} Mesh Elements")
plt.grid()
plt.savefig(f"temp_{n}.png",bbox_inches='tight')
plt.clf()

# dictioinaries for C/E plots
# key is number of mesh elements
# analytical quantities
analytical_phi = {}
analytical_T = {}
# C/E using VPP data from Cardinal
ratios_flux_num_to_analy = {}
ratios_temp_num_to_analy = {}
for n in n_elems:
    analytical_phi[n] = phi0*np.sqrt(1- ((lam - 1)*P*P*np.multiply(xx[n],xx[n]))/(L*L*q*q*phi0*phi0))
    analytical_T[n] = Sig_t0*T0*np.sqrt( (q*q*L*L*phi0*phi0)/(P*P) - (lam-1)*np.multiply(xx[n],xx[n]))
    ratios_temp_num_to_analy[n] = []
    ratios_flux_num_to_analy[n] = [] 
    for i in range(n):
        ratios_flux_num_to_analy[n].append(cardinal_fluxes[n][i]/analytical_phi[n][i]-1)
        ratios_temp_num_to_analy[n].append(temps[n][i]/analytical_T[n][i]-1)

# #plot individual C/E
# for n in n_elems:
#     plt.plot(xx[n],ratios_flux_num_to_analy[n],'-go',label=f"{n} x-elem")
#     plt.xticks([-60,-40,-20,0,20,40,60])
#     plt.yticks(np.linspace(np.min(ratios_flux_num_to_analy[n]),np.max(ratios_flux_num_to_analy[n]),5))
#     plt.xlabel("X Coordinate [cm]")
#     plt.ylabel(r"Flux C/E - 1")
#     plt.title(f"Ratio Numerical to Analytical Flux {n} Mesh Elements")
#     plt.legend()
#     plt.grid()
#     plt.savefig(f"num_to_analytical_flux_{n}.png",bbox_inches='tight')
#     plt.clf()

# plot all temp C/E on one plot
for n in n_elems:
    plt.plot(xx[n],ratios_temp_num_to_analy[n],label=f"{n} x-elem")
plt.xticks([-60,-40,-20,0,20,40,60])
dict_min = min(i for v in ratios_temp_num_to_analy.values() for i in v)
dict_max = max(i for v in ratios_temp_num_to_analy.values() for i in v)
plt.yticks(np.linspace(dict_min,dict_max,8))
plt.gca().yaxis.tick_right()
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
plt.xlabel("X Coordinate [cm]",fontsize=16)
plt.ylabel(r"Temperature C/E - 1",fontsize=16)
plt.yticks()
# plt.title(f"Ratio Numerical to Analytical Temp All Meshes. 200 Picard Iterations")
plt.legend(ncol=2)
plt.grid()
plt.savefig("temp_num_to_analy_ratios.png",bbox_inches='tight')
plt.clf()

# plot all flux C/E on one plot
for n in n_elems:
    plt.plot(xx[n],ratios_flux_num_to_analy[n],label=f"{n} x-elem")
plt.xticks([-60,-40,-20,0,20,40,60])
dict_min = min(i for v in ratios_flux_num_to_analy.values() for i in v)
dict_max = max(i for v in ratios_flux_num_to_analy.values() for i in v)
plt.yticks(np.linspace(dict_min,dict_max,8))
plt.gca().yaxis.tick_right()
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
plt.xlabel("X Coordinate [cm]",fontsize=16)
plt.ylabel(r"Flux C/E - 1",fontsize=16)
# plt.title("Ratio Numerical to Analytical Flux All Meshes. 200 Picard Iterations")
plt.legend(ncol=2)
plt.grid()
plt.savefig("flux_num_to_analy_ratios.png",bbox_inches='tight')
plt.clf()


#### ERROR NORMS ####

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

# compute flux error norm
flux_deviations = {}
flux_dev_norms = {}
analytical_norms = {}
flux_error_norms = {}
# flux error ratio
for n in n_elems:
    flux_deviations[n] = []
    # compute mean, stdev, and analytical - numerical
    for i in range(n):
        flux_deviations[n].append(analytical_phi[n][i]-cardinal_fluxes[n][i])
    flux_dev_norms[n] = np.linalg.norm(flux_deviations[n],ord=2)
    analytical_norms[n] = np.linalg.norm(analytical_phi[n],ord=2)
    # ratio is error norm
    flux_error_norms[n] = flux_dev_norms[n]/analytical_norms[n]

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
