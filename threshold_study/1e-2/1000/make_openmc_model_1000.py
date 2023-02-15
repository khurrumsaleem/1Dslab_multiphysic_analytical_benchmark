import openmc
import openmc.mgxs as mgxs
import numpy as np
import h5py

# problem physical parameters
T0 = 293
Tmax = 900
Tmin = T0 - 5 # lower bound to give data some range
L0 = 100
L = 106.47 # equilibrium length from paper
rho = 1.2 # g/cc
N_A = 6.022e23 # Avagadro's number
A = 180 # mass number for slab material
num_dens = rho*N_A/A
P = 1.0e22 # eV/s
q = 1e8 # eV
k0= 1.25e19 # eV/(s-cm-K^2) k(T) = k0 T(x)
phi0 = 2.5e14 # 1/s-cm^2 flux at the origin
s = 0.45 # Sigma_s/Sigma_t
f = 1.5 # nu Sigma_f/Sigma_t
nu = f/(1-s) # n per fission, this value comes from assuming there is no non-fisssion absorption and using the provided ratios
lam = 0.5*(1+np.sqrt(1+(16*q*q*phi0*phi0)/(P*P))) # eigenvalue solution
Sig_t0 = np.sqrt(P/((lam-1)*k0*L))/(T0) # 1/cm
sig_t0 = Sig_t0/num_dens # cm^2
# number of regions in the problem
N = 1000
infdim = 0.5 # length at which the reflective boundary conditions will be to simulate infiniteness in YZ dimension

# generate one group cross section data
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0.0, 20.0e6]) #  groups in eV

# create temperature data for cross sections
# isotropic angular data
temps = [float(T) for T in range(Tmin,Tmax+1)]
xsdata = openmc.XSdata('slab_xs', energy_groups=groups, temperatures=temps, num_delayed_groups=0)
xsdata.order = 0

# populate XS data for each temperature
for T in range(Tmin,Tmax+1):
    Sig_t = (Sig_t0 * T0) / T
    Sig_s = s*Sig_t
    nu_Sig_f = f*Sig_t
    Sig_f = nu_Sig_f / nu
    Sig_a = Sig_f # no non fission absorptions
    # add values to xsdata object
    xsdata.set_total(np.array([Sig_t]),temperature=T)
    xsdata.set_scatter_matrix(np.array([[[Sig_s]]]),temperature=T)
    xsdata.set_absorption(np.array([Sig_a]),temperature=T)
    xsdata.set_fission(np.array([Sig_f]),temperature=T)
    xsdata.set_nu_fission(np.array([nu_Sig_f]),temperature=T)
    xsdata.set_kappa_fission(np.array([q*Sig_f]),temperature=T)

# export xsdata
one_g_XS_file = openmc.MGXSLibrary(groups) # initialize the library
one_g_XS_file.add_xsdata(xsdata) # add benchmark XS data
one_g_XS_file.export_to_hdf5('one_gxs.h5') # write to disk

# create macroscopic object and export material with temperature XS library generated above
slab = openmc.Material(1, "slab")
slab.set_density('macro',1.)
slab.add_macroscopic('slab_xs')

materials = openmc.Materials([slab])
materials.cross_sections = 'one_gxs.h5'
materials.export_to_xml()

# create N x 1 x 1 regular mesh
mesh = openmc.RegularMesh()
mesh.dimension = (N,1,1) # N slices in the x dimension, 1 in y and z dimensions
mesh.lower_left = (-L/2,-infdim,-infdim)
mesh.upper_right = (L/2,infdim,infdim)
root_cell, cells = mesh.build_cells(['vacuum','vacuum','reflective','reflective','reflective','reflective']) # xmin , x max, ymin, ymax, zmin, zmax. internal surfaces transmission by default
# fill cells with slab
for cell in cells:
    cell.fill = slab
# create universe from filled cells, create geometry from universe and export the geometry
root_universe = openmc.Universe(name='root universe', cells=[root_cell])
geom = openmc.Geometry(root_universe)
geom.export_to_xml()

# create MGXS tallies
mesh_filter = openmc.MeshFilter(mesh)
tally = openmc.Tally()
tally.filters = [mesh_filter]
tally.scores = ['kappa-fission','flux']
mgxs_tallies = openmc.Tallies([tally])
mgxs_tallies.export_to_xml()

# settings
settings = openmc.Settings()
batches = 150
inactive = 50
particles = 25000
settings.energy_mode = 'multi-group'
settings.batches = batches
settings.inactive = inactive
settings.particles = particles
settings.output = {'tallies': True,'summary':False}

# Create a uniform spatial source distribution over fissionable zones  and a bidirectional angular distribution
bounds = [-L, -infdim, -infdim, L, infdim, infdim]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
mu = openmc.stats.Discrete([0.0],[1.0]) # exactly halfway between 1 and -1, i.e. in the xy plane
phi = openmc.stats.Discrete([0,np.pi],[0.5,0.5]) # either directly x or -x
bidirectional_x = openmc.stats.PolarAzimuthal(mu=mu,phi=phi)
settings.source = openmc.Source(space=uniform_dist,angle=bidirectional_x)
settings.temperature = {'default': T0,
                        'method': 'nearest',
                        'range': (Tmin, Tmax)} # good to load all temperatures you could encounter in multiphysics
settings.export_to_xml()
