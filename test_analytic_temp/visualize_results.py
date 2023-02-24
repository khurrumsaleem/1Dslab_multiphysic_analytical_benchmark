import openmc
import numpy as np
import matplotlib.pyplot as plt

sp = openmc.StatePoint('statepoint.150.h5')
raw_flux_mesh_tally = sp.get_tally(id=1).get_slice(scores=['flux'])
nu_fission_rate = float(sp.get_tally(id=2).get_slice(scores=['nu-fission']).mean[0])
fission_rate = float(sp.get_tally(id=2).get_slice(scores=['fission']).mean[0])
eig = sp.keff.n
power = sp.get_tally(id=2).get_slice(scores=['kappa-fission']).mean[0]
voxel_volume = L/N
source_strength = P*nu_fission_rate/(eig*voxel_volume*float(power))