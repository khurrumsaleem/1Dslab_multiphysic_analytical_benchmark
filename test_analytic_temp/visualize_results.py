import openmc
import numpy as np
import matplotlib.pyplot as plt
P = 1e22
L = 106.47
N = 50
phi0 = 2.5e14 # n/cm^2-s
q = 1e8 # ev/s
lam = 0.5*(1+np.sqrt(1+(16*q*q*phi0*phi0)/(P*P)))

#open statepoint
sp = openmc.StatePoint('statepoint.150.h5')

# obtain raw tallies (need source strength normalization)
raw_flux_mesh_tally = sp.get_tally(id=1).get_slice(scores=['flux'])
raw_kappa_fission_tally = sp.get_tally(id=1).get_slice(scores=['kappa-fission'])

# compute source strength
nu_fission_rate = float(sp.get_tally(id=2).get_slice(scores=['nu-fission']).mean[0])
fission_rate = float(sp.get_tally(id=2).get_slice(scores=['fission']).mean[0])
eig = sp.keff.n
power = sp.get_tally(id=2).get_slice(scores=['kappa-fission']).mean[0]
voxel_volume = L/N
source_strength = P*nu_fission_rate/(eig*voxel_volume*float(power))

# convert units for flux and kappa-fissiont allies, reshape to be 1D arrays
flux = raw_flux_mesh_tally * source_strength
kf = raw_kappa_fission_tally * source_strength
flux.mean.shape = (N)
kf.mean.shape = (N)

# plot
xx = np.linspace(-L/2,L/2,N)
plt.plot(xx,flux.mean,'-bo')
plt.xlabel("x coordinate")
plt.ylabel("flux [n/$cm^{2}$-s]")
plt.title("Flux Computed from Standaone OpenMC with the \n Analytical Temperature Set in Each Cell")
plt.savefig("flux_temp_analytical.png")
plt.clf()

plt.plot(xx,kf.mean,'-go')
plt.xlabel("x coordinate")
plt.ylabel("kappa fission rate")
plt.title("Kappa Fission Computed from Standaone OpenMC with the \n Analytical Temperature Set in Each Cell")
plt.savefig("kappa_fission_temp_analytical.png")
plt.clf()

num_to_analy = []
# flux numerical to analytical ratio
phi = phi0*np.sqrt(1- ((lam - 1)*P*P*np.multiply(xx,xx))/(L*L*q*q*phi0*phi0))
for i in range(N):
    num_to_analy.append(float(flux.mean[i])/phi[i])
plt.plot(xx,num_to_analy,'-co',)
plt.xlabel('X coordniate')
plt.ylabel('numerical to analytical ratio $\phi(x)$')
plt.title('Ratio Numerical to Analytical Flux from Standaone OpenMC with the \n Analytical Temperature Set in Each Cell')
plt.savefig('temp_analytical_num_flux_to_analytical_ratio.png')
plt.clf()
