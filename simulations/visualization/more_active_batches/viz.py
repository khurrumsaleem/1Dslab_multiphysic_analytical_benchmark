# a tidier visualization script that takes Cardinal outputs
# and plots C/E and error norms for the temperature

import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker as mtick
import pandas as pd
import scipy.stats as stats

# list mesh element sizes that will be used to iterate through all cases
n_elems = np.array([5,10,25,50,100,250,500,1000])
log_N = np.log10(n_elems)
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
# ones = np.ones(len(xx[1000]))
ones = np.ones(len(xx[1000]))

# plot 50 mesh element temperature
plt.plot(xx[50],temps[50],'-ro')
plt.xticks([-60,-40,-20,0,20,40,60])
plt.yticks(np.linspace(310,345,8))
plt.xlabel("X coordinate [cm]",fontsize=16)
plt.ylabel("Temperature [K]",fontsize=16)
# plt.title(f"Temperature for {50} Mesh Elements")
plt.grid()
plt.savefig(f"temp_50.png",bbox_inches='tight')
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
    ratios_temp_num_to_analy[n] = np.zeros(n)
    ratios_flux_num_to_analy[n] = np.zeros(n)
    for i in range(n):
        ratios_flux_num_to_analy[n][i] = cardinal_fluxes[n][i]/analytical_phi[n][i]
        ratios_temp_num_to_analy[n][i] = temps[n][i]/analytical_T[n][i]

# plot analytical flux for 50
plt.plot(xx[50],analytical_phi[50],'-o')
plt.xlabel("X Coordinate [cm]")
plt.ylabel("Flux [n/cm^2-s]")
plt.grid()
plt.savefig("analytical_flux_50.png")
plt.clf()

# plot all temp C/E on one plot
for n in n_elems:
    plt.plot(xx[n],ratios_temp_num_to_analy[n],label=f"{n} x-elem")
plt.plot(xx[1000],ones,'k',label="exact")
plt.xticks([-60,-40,-20,0,20,40,60])
dict_min = min(i for v in ratios_temp_num_to_analy.values() for i in v)
dict_max = max(i for v in ratios_temp_num_to_analy.values() for i in v)
# plt.yticks(np.linspace(dict_min,dict_max,8))
plt.yticks([0.995,0.999,0.9995,1.0,1.005])
plt.xlabel("X Coordinate [cm]",fontsize=16)
plt.ylabel(r"Temperature C/E",fontsize=16)
plt.gca().yaxis.tick_right()
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
plt.yticks()
# plt.title(f"Ratio Numerical to Analytical Temp All Meshes. 200 Picard Iterations")
plt.legend(ncol=2)
plt.grid()
plt.savefig("temp_num_to_analy_ratios.png",bbox_inches='tight')
plt.clf()

# COARSE TEMP C/E
for n in [5,10,25]:
    plt.plot(xx[n],ratios_temp_num_to_analy[n],label=f"{n} x-elem")
plt.plot(xx[1000],ones,'k',label="exact")
plt.xticks([-60,-40,-20,0,20,40,60])
dict_min = min(i for v in ratios_temp_num_to_analy.values() for i in v)
dict_max = max(i for v in ratios_temp_num_to_analy.values() for i in v)
# plt.yticks(np.linspace(dict_min,dict_max,8))
plt.yticks([0.995,0.996,0.997,0.998,0.999,1,1.001])
plt.xlabel("X Coordinate [cm]",fontsize=16)
plt.ylabel(r"Temperature C/E",fontsize=16)
plt.gca().yaxis.tick_right()
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
plt.yticks()
# plt.title(f"Ratio Numerical to Analytical Temp All Meshes. 200 Picard Iterations")
plt.legend(loc='center',ncol=2)
plt.grid()
plt.savefig("coarse_temp_num_to_analy_ratios.png",bbox_inches='tight')
plt.clf()

# FINE TEMP C/E
for n in [50,100,250,500,1000]:
    plt.plot(xx[n],ratios_temp_num_to_analy[n],label=f"{n} x-elem")
plt.plot(xx[1000],ones,'k',label="exact")
plt.xticks([-60,-40,-20,0,20,40,60])
dict_min = min(i for v in ratios_temp_num_to_analy.values() for i in v)
dict_max = max(i for v in ratios_temp_num_to_analy.values() for i in v)
# plt.yticks(np.linspace(dict_min,dict_max,8))
plt.yticks([0.99994,0.99995,0.99996,0.99997,0.99998,0.99999,1,1.00001,1.00002])
plt.xlabel("X Coordinate [cm]",fontsize=16)
plt.ylabel(r"Temperature C/E",fontsize=16)
plt.gca().yaxis.tick_right()
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5f'))
plt.yticks()
# plt.title(f"Ratio Numerical to Analytical Temp All Meshes. 200 Picard Iterations")
plt.legend(ncol=2)
plt.grid()
plt.savefig("fine_temp_num_to_analy_ratios.png",bbox_inches='tight')
plt.clf()

#### TEMP ERROR NORM ####

log_temp_error_norms = [] # list to house the error for each mesh size
for n in n_elems:
    filename = f"openmc_{n}_out.csv"
    df_norms = pd.read_csv(filename)
    last_val  = df_norms.iloc[len(df_norms)-1,-1]
    log_temp_error_norms.append(np.log10(last_val))
print(log_temp_error_norms)
# plot the error for each threshold against the number of mesh elements
plt.plot(log_N,log_temp_error_norms,'-ro',label=r'$\epsilon_{T}=\frac{||T_{a} - T_{x} ||_{2}}{|| T_{a} ||_{2}}$')
plt.xticks(log_N)
plt.yticks([-1.75,-2,-2.25,-2.5,-2.75,-3,-3.25,-3.5,-3.75,-4,-4.25])
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xlabel(r"Log of Number of X Elements",fontsize=16)
plt.ylabel(r"$log(\epsilon_{T})$",fontsize=16)
# plt.title("Temperature Error Norm vs Mesh Size \n 200 Picard Iterations",y=1)
plt.legend(bbox_to_anchor=[0.985,0.985],fontsize=14)
plt.grid()
plt.savefig("temp_error_norms.png",bbox_inches='tight')
plt.clf()

results_T = stats.linregress(log_N,log_temp_error_norms)
print("slope = ", results_T.slope, " r-squared = ",results_T.rvalue*results_T.rvalue)