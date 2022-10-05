import openmc
import openmc.mgxs as mgxs
import numpy as np
import h5py

# problem physical parameters
T0 = 293
Tmax = 900
L0 = 100
L = 106.47 # equilibrium length from paper (TODO perhaps use formula)
rho = 1.2 # g/cc
N_A = 6.022e23 # Avagadro's number
A = 180 # mass number for slab material
num_dens = rho*N_A/A
P = 1.0e22 # eV/s
q = 1e8 # eV
k0= 1.25e19 # eV/(s-cm-K^2)
phi0 = 2.5e14 # 1/s-cm^2 flux at the origin
s = 0.45 # Sigma_s/Sigma_t
f = 1.5 # nu Sigma_f/Sigma_t
nu = 30/11 # neutrons per fission
lam = 0.5*(1+np.sqrt(1+(16*q*q*phi0*phi0)/(P*P))) # eigenvalue solution
Sig_t0 = np.sqrt(P/((lam-1)*k0*L))/(T0) #macro XS
sig_t0 = Sig_t0/num_dens # micro XS
# number of regions in the problem
N = 4
infdim = 50.0 # length at which the reflective boundary conditions will be to simulate infiniteness in YZ dimension

# generate one group cross section data
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0.0, 20.0e6]) #  groups in eV

# create temperature data for cross sections
temps = [float(T) for T in range(T0,Tmax+1)]
xsdata = openmc.XSdata('slab_xs', energy_groups=groups, temperatures=temps,num_delayed_groups=0) # ISOTROPIC
xsdata.order = 0

# Above isotropic, below S2. Need s2 for paper, but curious about isotropic results
xsdata_s2 = openmc.XSdata('slab_xs_s2', energy_groups=groups, temperatures=temps, representation='angle', num_delayed_groups=0) # ANGUlAR DIST
xsdata_s2.scatter_format = 'tabular' #set a pmf where -1 and 1 are only options for mu
xsdata_s2.order = 2 # two points, -1 and 1

# populate XS data
for T in range(T0,Tmax+1):
    Sig_t = (Sig_t0 * T0) / T
    Sig_s = s*Sig_t
    nu_Sig_f = f*Sig_t
    # isotropic
    xsdata.set_total(np.array([Sig_t]),temperature=T)
    xsdata.set_scatter_matrix(np.array([Sig_s]),temperature=T) # TODO figure out scatter matrix
    xsdata.set_nu_fission(np.array([nu_Sig_f]),temperature=T) # TODO this vs set_fission
    # s2
    xsdata_s2.set_total(np.array([Sig_t]),temperature=T)
    xsdata_s2.set_scatter_matrix(np.array([Sig_s]),temperature=T) # TODO figure out scatter matrix
    xsdata_s2.set_nu_fission(np.array([nu_Sig_f]),temperature=T) # TODO this vs set_fission

# export xsdata
one_g_cross_sections_file = openmc.MGXSLibrary(groups) # initialize the library
# one_g_cross_sections_file.add_xsdata(xsdata) # add benchmark XS data
one_g_cross_sections_file.add_xsdata(xsdata_s2) # add benchmark XS data
one_g_cross_sections_file.export_to_hdf5('one_gxs.h5') # write to disk

# create macroscopic object and export materials
slab = openmc.Material(1, "slab")
slab.set_density('macro',1.)
slab.add_macroscopic('slab_xs_s2')

materials = openmc.Materials([slab])
materials.cross_sections = 'one_gxs.h5'
materials.export_to_xml()

# TODO I think this is obsolete now with the addition of the Macroscopic object
# create material for slab
# slab = openmc.Material(1, "slab",T0)
# slab.set_density('atom/b-cm',num_dens) #probably want this one we don't want to load the nuclear data for A=180
# mats = openmc.Materials([slab])
# mats.export_to_xml()

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

# settings
settings = openmc.Settings()
batches = 100
inactive = 30
particles = 2000

settings.batches = batches
settings.inactive = inactive
settings.particles = particles
settings.output = {'tallies': True}

# Create a uniform spatial source distribution over fissionable zones
bounds = [-L, -infdim, -infdim, L, infdim, infdim]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.Source(space=uniform_dist)
settings.temperature = {'default': T0,
                        'method': 'interpolation',
                        'range': (0.0, 900.0)} # good to load all temperatures you could encounter in multiphysics
settings.export_to_xml()

# create MGXS tallies
mgxs_tallies = openmc.Tallies()
