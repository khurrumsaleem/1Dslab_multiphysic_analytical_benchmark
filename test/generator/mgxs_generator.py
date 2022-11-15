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
batches = 20
inactive = 10
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

# create mgxs library object
onegxs_lib = openmc.mgxs.Library(geometry) # initialize the library
onegxs_lib.energy_groups = groups
onegxs_lib.mgxs_types = ['total', 'absorption', 'nu-fission', 'fission',
                       'nu-scatter matrix', 'multiplicity matrix', 'chi']
# Specify a "cell" domain type for the cross section tally filters
onegxs_lib.domain_type = "material"
# Specify the cell domains over which to compute multi-group cross sections
onegxs_lib.domains = geometry.get_all_materials().values()
# Set the Legendre order to 3 for P3 scattering
onegxs_lib.legendre_order = 3
# check library and build
onegxs_lib.check_library_for_openmc_mgxs()
onegxs_lib.build_library()

# run model
sp = model.run()

# load statepoint
statepoint = openmc.StatePoint(sp, autolink=True)

# merge tallies and load from statepoint
onegxs_lib.load_from_statepoint(statepoint)
onegxs_file = onegxs_lib.create_mg_library(xs_type='macro')
# onegxs_file.export_to_hdf5()
onegxs_lib.build_hdf5_store(filename='onegxs.h5',directory="/home/lgross/1Dslab_multiphysic_analytical_benchmark/test")