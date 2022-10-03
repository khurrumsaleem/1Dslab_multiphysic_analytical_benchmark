import openmc
import openmc.mgxs as mgxs
import numpy as np
# problem parameters
T0 = 293
L = 100.0 # TODO determine which L, initial or equilibrium
infdim = 50.0 # length at which the reflective boundary conditions will be to simulate infiniteness in YZ dimension
rho = 1.2 # g/cc
N = 4 # number of regions in the problem
N_A = 6.022e23 # Avagadro's number
A = 180 # mass number for slab material
num_dens = rho*N_A/A

# create material for slab
slab = openmc.Material(1, "slab",T0)
slab.set_density('atom/b-cm',num_dens) #probably want this one we don't want to load the nuclear data for A=180 
mats = openmc.Materials([slab])
mats.export_to_xml()

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
root_universe = openmc.Universe(name='root universe', cells=cells)
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

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-L, -infdim, -infdim, L, infdim, infdim]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.Source(space=uniform_dist)
settings.temperature = {'default': T0,
                        'method': 'interpolation',
                        'range': (0.0, 900.0)} # good to load all temperatures you could encounter in multiphysics
settings.export_to_xml()

# generate one group cross sections
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0.0, 20.0e6]) #  groups in eV

# reaction types
total = mgxs.TotalXS(domain=root_universe, energy_groups=groups)
absorption = mgxs.AbsorptionXS(domain=root_universe, energy_groups=groups)
scattering = mgxs.ScatterXS(domain=root_universe, energy_groups=groups)
fission = mgxs.FissionXS(domain=root_universe, energy_groups=groups, nu=True)

# create MGXS tallies 
mgxs_tallies = openmc.Tallies()
mgxs_tallies.append(total.tallies['flux'])
mgxs_tallies.append(absorption.tallies['flux'])
mgxs_tallies.append(scattering.tallies['flux'])
mgxs_tallies.append(fission.tallies['flux'])
mgxs_tallies.export_to_xml()