import openmc
import openmc.lib
import openmc.mgxs as mgxs
import numpy as np
import h5py

Tmin = 293
T0 = 303 # temp in range that should be set-able
Tmax = 313

# create macroscopic object and export material with temperature XS library generated below
cold_fuel = openmc.Material(material_id=1,name="cold_fuel")
hot_fuel = openmc.Material(material_id=2,name="hot_fuel")
# default names for the cross section set used in the CE file to generate the data are set1,set2,...
# there are two cells in the model and we computed XS for each
cold_fuel.add_macroscopic('set1')
hot_fuel.add_macroscopic('set2')
materials = openmc.Materials([cold_fuel,hot_fuel])
materials.cross_sections = 'onegxs.h5'
materials.export_to_xml()

# define geometry
H = 100.0 # cm
xplane0 = openmc.XPlane(x0  = 0.0, boundary_type='vacuum')
xmidplane = openmc.XPlane(x0  = H/2, boundary_type='transmission')
xplane1 = openmc.XPlane(x0 = H, boundary_type='vacuum')
yplane0 = openmc.YPlane(y0 = 0.0, boundary_type='vacuum')
yplane1 = openmc.YPlane(y0 = H, boundary_type='vacuum')
zplane0 = openmc.ZPlane(z0 = 0.0, boundary_type='vacuum')
zplane1 = openmc.ZPlane(z0 = H, boundary_type='vacuum')
cold_region = +xplane0 & -xmidplane & +yplane0 & -yplane1 & +zplane0 & -zplane1
hot_region = +xmidplane & -xplane1 & +yplane0 & -yplane1 & +zplane0 & -zplane1
cold_cell = openmc.Cell(cell_id=0, region=cold_region, fill=cold_fuel)
cold_cell.temperature = Tmin
hot_cell = openmc.Cell(cell_id=1, region=hot_region, fill=hot_fuel)
hot_cell.temperature = Tmax
geometry = openmc.Geometry([cold_cell,hot_cell])
geometry.export_to_xml()

# define settings
settings = openmc.Settings()
batches = 20
inactive = 10
particles = 5000
settings.energy_mode = 'multi-group'
settings.batches = batches
settings.inactive = inactive
settings.particles = particles
settings.temperature = {'default': T0,
                        'method': 'interpolation',
                        'range': (Tmin, Tmax)} # good to load all temperatures you could encounter in multiphysics
bounds = [0, 0, 0, H, H, H]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.Source(space=uniform_dist)
settings.export_to_xml()

# initialize an openmc simulation
openmc.lib.init()
# C-API call to set (valid) temperature of cell that results in error
openmc.lib.Cell.set_temperature(self=openmc.lib.cells[0],T=T0,set_contained=True)
# openmc.lib.run()
# openmc.lib.finalize()