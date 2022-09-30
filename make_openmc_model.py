import openmc
import openmc.mgxs as mgxs
import numpy as np
# problem parameters
T0 = 293
L = 100.0 # TODO determine which L, initial or equilibrium
infdim = 50.0 # length at which the reflective boundary conditions will be to simulate infiniteness in YZ dimension
rho = 1.2 # g/cc
N = 4 # number of regions in the problem

# create material for slab
slab = openmc.Material(1, "slab",T0)
slab.set_density('g/cm3',rho)
mats = openmc.Materials([slab])
mats.export_to_xml()

# define mesh tally MAYBE DONT SINCE EACH MGXS object has a .tallies to use instead?
mesh = openmc.RegularMesh()
mesh.dimension = (N,1,1) # N slices in the x dimension, 1 in y and z dimensions
mesh.lower_left = (-L/2,-infdim,-infdim)
mesh.upper_right = (L/2,infdim,infdim)
# mesh_filter = openmc.MeshFilter(mesh) TBD if we want this
mesh.build_cells(['vacuum','vacuum','reflective','reflective','reflective','reflective']) # xmin , x max, ymin, ymax, zmin, zmax TODO are the internal ones transmission by default?
mesh.export_to_xml() # TODO how to export mesh version of geometry is this enough?

# # TODO determine if mesh is better and we create our own cells but save the below just incase
# # boundary planes of the problem
# x_0 = openmc.XPlane(x0 = -L/2, boundary_type='vacuum')
# x_N = openmc.XPlane(x0 = L/2, boundary_type='vacuum')
# y_minus = openmc.YPlane(y0 = -infdim, boundary_type='reflective')
# y_plus = openmc.YPlane(y0 = infdim, boundary_type='reflective')
# z_minus = openmc.ZPlane(z0 = -infdim, boundary_type='reflective')
# z_plus = openmc.ZPlane(z0 = infdim, boundary_type='reflective')


# x_1 = openmc.XPlane(x0 = -L/4 , boundary_type='transmission')
# x_2 = openmc.XPlane(x0 = 0 , boundary_type='transmission')
# x_3 = openmc.XPlane(x0 = L/4, boundary_type='transmission')


# # create regions
# region_1 = +x_0 & -x_1 & +y_minus & -y_plus & +z_minus & -z_plus
# region_2 = +x_1 & -x_2 & +y_minus & -y_plus & +z_minus & -z_plus
# region_3 = +x_2 & -x_3 & +y_minus & -y_plus & +z_minus & -z_plus
# region_4 = +x_3 & -x_N & +y_minus & -y_plus & +z_minus & -z_plus

# # create cells
# cell_1 = openmc.Cell(cell_id=1, fill=uo2, region = region_1 )
# cell_2 = openmc.Cell(cell_id=2, fill=uo2, region = region_2 )
# cell_3 = openmc.Cell(cell_id=3, fill=uo2, region = region_3 )
# cell_4 = openmc.Cell(cell_id=4, fill=uo2, region = region_4 )

# # define universe of cells
# root_universe = openmc.Universe(name='root universe', cells=[cell_1,cell_2,cell_3,cell_4])
# geom = openmc.Geometry([cell_1,cell_2,cell_3,cell_4])
# geom.export_to_xml()

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
bounds = [-L -infdim, -infdim, L, infdim, infdim]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.Source(space=uniform_dist)
settings.temperature = {'default': T0,
                        'method': 'interpolation',
                        'range': (294.0, 1600.0)} # good to load all temperatures you could encounter in multiphysics
settings.export_to_xml()

# generate one group cross sections
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0.0, 20.0e6]) #  groups in eV

# domain (openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh) – The domain for spatial homogenization
# domain_type ({'material', 'cell', 'distribcell', 'universe', 'mesh'}) – The domain type for spatial homogenization

# Question: does specifying domain = root_universe command openmc to compute the XS in each cell?
total = mgxs.TotalXS(domain=root_universe, groups=groups)
absorption = mgxs.AbsorptionXS(domain=root_universe, groups=groups)
scattering = mgxs.ScatterXS(domain=root_universe, groups=groups)
fission = mgxs.FissionXS(domain=root_universe, groups=groups, nu=True)

# Question how do we get our library in here, especially since I didn't provide element information in the materials




# tally = openmc.Tally(tally_id=1,name="mesh_tally")
# tally.filters = [mesh_filter]
# tally.scores = ["flux"]
# tallies = openmc.Tallies([tally])
# tallies.export_to_xml()

mgxs_tallies = openmc.Tallies()
mgxs_tallies += total.tallies.values()
mgxs_tallies += absorption.tallies.values()
mgxs_tallies += scattering.tallies.values()
mgxs_tallies += fission.tallies.values()
mgxs_tallies.export_to_xml()
