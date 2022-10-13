import openmc
import openmc.mgxs as mgxs
import numpy as np
import h5py

# problem physical parameters
T0 = 293
Tmax = 900
L0 = 100
L = 106.47 # equilibrium length from paper
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
nu = 30/11 # n per fission, this value comes from assuming there is no non-fisssion absorption and using the provided ratios
lam = 0.5*(1+np.sqrt(1+(16*q*q*phi0*phi0)/(P*P))) # eigenvalue solution
Sig_t0 = np.sqrt(P/((lam-1)*k0*L))/(T0) #macro XS
sig_t0 = Sig_t0/num_dens # micro XS
# number of regions in the problem
N = 4
infdim = 50.0 # length at which the reflective boundary conditions will be to simulate infiniteness in YZ dimension

model = openmc.Model()

# generate one group cross section data
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0.0, 20.0e6]) #  groups in eV

# create temperature data for cross sections
# isotropic angular data
temps = [float(T) for T in range(T0,Tmax+1)]
xsdata = openmc.XSdata('slab_xs', energy_groups=groups, temperatures=temps, num_delayed_groups=0)
xsdata.order = 0

# populate XS data for each temperature
for T in range(T0,Tmax+1):
    Sig_t = (Sig_t0 * T0) / T
    Sig_s = s*Sig_t
    nu_Sig_f = f*Sig_t
    Sig_a = nu_Sig_f / nu # assuming no non-fisssion absorption, Sig_a = Sig_f or Sig_A = Sig_t * f / nu
    # isotropic
    xsdata.set_total(np.array([Sig_t]),temperature=T)
    xsdata.set_scatter_matrix(np.array([[[Sig_s]]]),temperature=T)
    xsdata.set_absorption(np.array([Sig_a]),temperature=T)
    xsdata.set_nu_fission(np.array([nu_Sig_f]),temperature=T)

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

model.materials = materials

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

model.geometry = geom

# create MGXS tallies
mesh_filter = openmc.MeshFilter(mesh)
tally = openmc.Tally(tally_id=1, name="mesh_tally")
tally.filters = [mesh_filter]
tally.scores = ['flux']
mgxs_tallies = openmc.Tallies([tally])
mgxs_tallies.export_to_xml()

model.tallies = mgxs_tallies

# settings
settings = openmc.Settings()
batches = 150
inactive = 50
particles = 3000
# settings.run_mode = 'fixed source' # TODO run in fixed source, see what tracks look like (ensures patch and source angular distribtuion is working
settings.energy_mode = 'multi-group'
settings.batches = batches
settings.inactive = inactive
settings.particles = particles
settings.output = {'tallies': True}

# Create a uniform spatial source distribution over fissionable zones  and a bidirectional angular distribution
bounds = [-L, -infdim, -infdim, L, infdim, infdim]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
mu = openmc.stats.Discrete([0.0],[1.0])
phi = openmc.stats.Discrete([0,np.pi],[0.5,0.5])
bidirectional_x = openmc.stats.PolarAzimuthal(mu=mu,phi=phi)
settings.source = openmc.Source(space=uniform_dist,angle=bidirectional_x)
settings.temperature = {'default': T0,
                        'method': 'interpolation',
                        'range': (0.0, 900.0)} # good to load all temperatures you could encounter in multiphysics
settings.max_tracks = 100
settings.export_to_xml()

model.settings = settings

statepoint_file = model.run(tracks=True)