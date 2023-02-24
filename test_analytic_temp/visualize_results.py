import openmc
import numpy as np
import matplotlib.pyplot as plt
P = 1e22
L = 106.47
N = 50
sp = openmc.StatePoint('statepoint.150.h5')
raw_flux_mesh_tally = sp.get_tally(id=1).get_slice(scores=['flux'])
raw_kappa_fission_tally = sp.get_tally(id=1).get_slice(scores=['kappa-fission'])
nu_fission_rate = float(sp.get_tally(id=2).get_slice(scores=['nu-fission']).mean[0])
fission_rate = float(sp.get_tally(id=2).get_slice(scores=['fission']).mean[0])
eig = sp.keff.n
power = sp.get_tally(id=2).get_slice(scores=['kappa-fission']).mean[0]
voxel_volume = L/N
source_strength = P*nu_fission_rate/(eig*voxel_volume*float(power))
flux = raw_flux_mesh_tally * source_strength
kf = raw_kappa_fission_tally * source_strength
flux.mean.shape = (N)
kf.mean.shape = (N)
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
plt.savefig("kappa_fission_temp_analytical.png")
plt.clf()