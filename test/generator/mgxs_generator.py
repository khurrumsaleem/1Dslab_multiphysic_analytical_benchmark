import openmc
import openmc.lib
import openmc.mgxs as mgxs
import numpy as np
import h5py


# create OpenMC model
model = openmc.Model()

Tmin = 293
T0 = 303 # temp in range that should be set-able
Tmax = 313

# create macroscopic object and export material with temperature XS library generated below
fuel = openmc.Material(material_id=1,name="fuel")
fuel.add_element('U', 1., enrichment=3.0)
fuel.set_density('macro',18.9)
materials = openmc.Materials([fuel])
model.materials = materials

# define geometry
H = 100.0 # cm
xplane0 = openmc.XPlane(x0  = 0.0, boundary_type='vacuum')
xplane1 = openmc.XPlane(x0 = H, boundary_type='vacuum')
yplane0 = openmc.YPlane(y0 = 0.0, boundary_type='vacuum')
yplane1 = openmc.YPlane(y0 = H, boundary_type='vacuum')
zplane0 = openmc.ZPlane(z0 = 0.0, boundary_type='vacuum')
zplane1 = openmc.ZPlane(z0 = H, boundary_type='vacuum')
cube_region = +xplane0 & -xplane1 & +yplane0 & -yplane1 & +zplane0 & -zplane1
root_cell = openmc.Cell(cell_id=0, region=cube_region, fill=fuel)
geometry = openmc.Geometry([root_cell])
model.geometry = geometry

# define settings
settings = openmc.Settings()
batches = 200
inactive = 50
particles = 5000
settings.batches = batches
settings.inactive = inactive
settings.particles = particles
settings.temperature = {'default': T0,
                        'method': 'interpolation',
                        'range': (Tmin, Tmax)} # good to load all temperatures you could encounter in multiphysics
bounds = [0, 0, 0, H, H, H]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.Source(space=uniform_dist)
settings.output = {'tallies': True,'summary':True}
model.settings = settings

# generate multigroup objects and tallies
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0.0, 20.0e6]) #  groups in eV

total = mgxs.TotalXS(domain=root_cell, energy_groups=groups)
absorption = mgxs.AbsorptionXS(domain=root_cell, energy_groups=groups)
scattering = mgxs.ScatterXS(domain=root_cell, energy_groups=groups)
fission = mgxs.FissionXS(domain=root_cell, energy_groups=groups)
chi = mgxs.Chi(domain=root_cell, energy_groups=groups)

tallies = openmc.Tallies()
tallies += total.tallies.values()
tallies += absorption.tallies.values()
tallies += scattering.tallies.values()
tallies += fission.tallies.values()
tallies += chi.tallies.values()
model.tallies = tallies

# run model
model.run()

