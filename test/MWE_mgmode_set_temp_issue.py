import openmc
import openmc.lib
import openmc.mgxs as mgxs
import numpy as np
import h5py

Tmin = 293
T0 = 303 # temp in range that should be set-able
Tmax = 313

# create macroscopic object and export material with temperature XS library generated below
fuel = openmc.Material(material_id=1,name="fuel")
fuel.add_element('U', 1., enrichment=3.0)
fuel.set_density('macro',18.9)
materials = openmc.Materials([fuel])
materials.cross_sections = 'onegxs.h5'
materials.export_to_xml()

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
geometry.export_to_xml()

# define settings
settings = openmc.Settings()
batches = 200
inactive = 50
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

# generate one group library
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0.0, 20.0e6]) #  groups in eV

sp = openmc.StatePoint('./generator/statepoint.200.h5', autolink=True)
print(sp.get_tally())

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# # METHOD 1  https://nbviewer.org/github/openmc-dev/openmc-notebooks/blob/main/mgxs-part-i.ipynb
# # error:  ERROR: Attribute "energy_groups" does not exist in object / ... confusing becasue I set them right below here
# total = mgxs.TotalXS(domain=root_cell, energy_groups=groups)
# absorption = mgxs.AbsorptionXS(domain=root_cell, energy_groups=groups)
# scattering = mgxs.ScatterXS(domain=root_cell, energy_groups=groups)
# fission = mgxs.FissionXS(domain=root_cell, energy_groups=groups)
# chi = mgxs.Chi(domain=root_cell, energy_groups=groups)

# total.load_from_statepoint(sp)
# absorption.load_from_statepoint(sp)
# scattering.load_from_statepoint(sp)
# fission.load_from_statepoint(sp)
# chi.load_from_statepoint(sp)

# total.print_xs()

# total.build_hdf5_store(filename='onegxs.h5',directory="/home/lgross/1Dslab_multiphysic_analytical_benchmark/test",append=True)
# absorption.build_hdf5_store(filename='onegxs.h5',directory="/home/lgross/1Dslab_multiphysic_analytical_benchmark/test",append=True)
# scattering.build_hdf5_store(filename='onegxs.h5',directory="/home/lgross/1Dslab_multiphysic_analytical_benchmark/test",append=True)
# fission.build_hdf5_store(filename='onegxs.h5',directory="/home/lgross/1Dslab_multiphysic_analytical_benchmark/test",append=True)
# chi.build_hdf5_store(filename='onegxs.h5',directory="/home/lgross/1Dslab_multiphysic_analytical_benchmark/test",append=True)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# METHOD 2 https://nbviewer.org/github/openmc-dev/openmc-notebooks/blob/main/mg-mode-part-ii.ipynb

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
# merge tallies and load from statepoint
onegxs_lib.load_from_statepoint(sp)
onegxs_file = onegxs_lib.create_mg_library(xs_type='macro')
onegxs_file.export_to_hdf5()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# METHOD 3  https://nbviewer.org/github/openmc-dev/openmc-notebooks/blob/main/mgxs-part-iii.ipynb
onegxs_lib.build_hdf5_store(filename='onegxs.h5',directory="/home/lgross/1Dslab_multiphysic_analytical_benchmark/test")
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #

# initialize an openmc simulation
openmc.lib.init()
# C-API call to set (valid) temperature of cell that results in error
openmc.lib.Cell.set_temperature(self=openmc.lib.cells[0],T=T0,set_contained=True)
# openmc.lib.run()
# openmc.lib.finalize()