import openmc
import numpy as np
import matplotlib.pyplot as plt

# Load statepoint file and get the number of bins
sp = openmc.StatePoint('statepoint.500.h5')
tally = sp.get_tally()
num_voxels = tally.shape[0] # this works because the mesh is N x 1 x 1, so the total number of bins is N
# get scores for each flux
flux = tally.get_slice(scores=['flux'])

# paramters for analytical solution
T0 = 293
L = 106.47 # equilibrium length from paper (TODO perhaps use formula)
P = 1.0e22 # eV/s
q = 1e8 # eV
k0= 1.25e19 # eV/(s-cm-K^2) k(T) = k0 T(x)
phi0 = 2.5e14 # 1/s-cm^2 flux at the origin
eV_to_J = 1.602e-19 # J per eV
lam = 0.5*(1+np.sqrt(1+(16*q*q*phi0*phi0)/(P*P)))
# plot flux data using hist
flux.mean.shape = (num_voxels)
flux_min = min(flux.mean)
flux_max = max(flux.mean)
flux.std_dev.shape= (num_voxels)
bin_boundaries = np.linspace(-L,L,num_voxels+1)
# give one value in each bin paired with the bin boundaries.
# then scale the weight of each bin by the flux
plt.hist(bin_boundaries[:-1],bin_boundaries,weights=flux.mean,histtype='step') # plot flux vs position
# error bars for flux bins
bin_average_values = (bin_boundaries[1:] + bin_boundaries[:-1])/2 # get bin centers
plt.errorbar(bin_average_values,flux.mean,yerr=flux.std_dev,marker = '|',fmt='none',elinewidth=1,capsize=3,capthick=1)
plt.title('Flux in Each Mesh Element')
plt.xlabel('Position [cm]')
plt.ylabel('Flux [n/cm^2-s]')
plt.xlim([-(L+2),L+2])
plt.xticks(bin_boundaries[::5])
plt.ylim([flux_min-0.005,flux_max+0.005])
plt.savefig('flux_dist_50.png')
plt.clf()
# Determine relative error
flux_relative_error = np.zeros_like(flux.std_dev)
nonzero = flux.mean > 0
flux_relative_error[nonzero] = flux.std_dev[nonzero] / flux.mean[nonzero]
# distribution of relative errors
plt.hist(flux_relative_error[nonzero], bins=50)
plt.title('Relative Error Distribution for Flux')
plt.xlabel('Relative Error')
plt.ylabel('Frequency')
# TODO compute error from analytical solutioin
xx = np.linspace(-L,L,num_voxels)
analytical_flux = phi0*np.sqrt(np.ones(num_voxels)-((lam-1)*P*P*np.multiply(xx,xx))/(L*L*q*q*phi0*phi0))
# print(analytical_flux)
# print(max(analytical_flux))
plt.savefig('relative_errors_50.png')
plt.clf()
plt.plot(xx,analytical_flux,'-bo')
plt.xlabel('Position [cm]')
plt.ylabel('Flux [n/cm^2-s]')
plt.title('Analytical Flux for 50 Mesh Points')
plt.savefig('analytical_50.png')
#TODO scaling from total kappa fission tally to turn flux tally into correct units