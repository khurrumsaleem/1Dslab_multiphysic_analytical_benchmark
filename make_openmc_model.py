import openmc
import openmc.mgxs as mgxs

# problem parameters
T0 = 300
L = 1.0 # TODO determine which L, initial or equilibrium
infdim = 5.0 # length at which the reflective boundary conditions will be to simulate infiniteness in YZ dimension
rho = 10.0 # g/cc
N = 4 # number of regions in the problem

# create materials
# did Greisheimer and Kooreman use an actual material that would make this easier?
# if not maybe UO2?
uo2 = openmc.Material(1, "uo2")
uo2.add_element('U', 1.0, enrichment=5.0) # this will be used in computing lambda/the eigenvalue
uo2.add_element('O', 2.0)
uo2.set_density('g/cm3',rho)
mats = openmc.Materials([uo2])
mats.export_to_xml()

# TODO maybe change x_N to x_4 or get loop to work lol
# boundary planes of the problem
x_0 = openmc.XPlane(x0 = -L/2, boundary_type='vacuum')
x_N = openmc.XPlane(x0 = L/2, boundary_type='vacuum')
y_minus = openmc.YPlane(y0 = -infdim, boundary_type='reflective')
y_plus = openmc.YPlane(y0 = infdim, boundary_type='reflective')
z_minus = openmc.ZPlane(z0 = -infdim, boundary_type='reflective')
z_plus = openmc.ZPlane(z0 = infdim, boundary_type='reflective')


# TODO see if there's a way to use loops to do this, getting stuck on variable number of variables
# maybe precreate the correct number and loop to update, but not sure on this
# x planes to create N x 1 x 1 mesh, make
# using N = 4, so there need to be 3 internal planes
# most likely need more than 4 for the real deal, which motivates a loop
x_1 = openmc.XPlane(x0 = -L/4 , boundary_type='transmission')
x_2 = openmc.XPlane(x0 = 0 , boundary_type='transmission')
x_3 = openmc.XPlane(x0 = L/4, boundary_type='transmission')


# create regions
region_1 = +x_0 & -x_1 & +y_minus & -y_plus & +z_minus & -z_plus
region_2 = +x_1 & -x_2 & +y_minus & -y_plus & +z_minus & -z_plus
region_3 = +x_2 & -x_3 & +y_minus & -y_plus & +z_minus & -z_plus
region_4 = +x_3 & -x_N & +y_minus & -y_plus & +z_minus & -z_plus


# create cells
cell_1 = openmc.Cell(cell_id=1, fill=uo2, region = region_1 )
cell_2 = openmc.Cell(cell_id=2, fill=uo2, region = region_2 )
cell_3 = openmc.Cell(cell_id=3, fill=uo2, region = region_3 )
cell_4 = openmc.Cell(cell_id=4, fill=uo2, region = region_4 )

geom = openmc.Geometry([cell_1,cell_2,cell_3,cell_4])
geom.export_to_xml()

# define mesh tally
mesh = openmc.RegularMesh()
mesh.dimension = (N,1,1) # N slices in the x dimension, 1 in y and z dimensions
mesh.lower_left = (-L/2,-infdim,-infdim)
mesh.upper_right = (L/2,infdim,infdim)
mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(tally_id=1,name="mesh_tally")
tally.filters = [mesh_filter]
tally.scores = ["flux"]
tallies = openmc.Tallies([tally])
tallies.export_to_xml()


# settings
settings = openmc.Settings()
# more settings

# set temperatures (perhaps initial temp to T_0 in every cell and have it update)
settings.temperature = {'default': 600.0,
                        'method': 'interpolation',
                        'range': (294.0, 1600.0)} # good to load all temperatures you could encounter in multiphysics


settings.export_to_xml()