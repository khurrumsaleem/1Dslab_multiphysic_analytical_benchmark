# GLOBAL VARS
# problem physical parameters
T0 = 293
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
sig_t0 = np.sqrt(P/((lam-1)*k0*L))/(num_dens*T0)
# number of regions in the problem
N = 4
infdim = 50.0 # length at which the reflective boundary conditions will be to simulate infiniteness in YZ dimension

  [Mesh]
    [centered_mesh]
      type=GeneratedMeshGenerator
      dim = 3
      xmax = ${fparse -1*L/2}
      xmin = ${fparse L/2}
      ymax = ${fparse -1*infdim}
      ymin = ${fparse infdim}
      zmax = ${fparse -1*infdim}
      zmin = ${fparse infdim}
      nx = ${fparse N}
      ny = 1
      nz = 1
      length_unit = 'cm'
    []
  []

[AuxVariables]
   # always set
   [temp]
       family = MONOMIAL
       order = CONSTANT
   []
   [heat_source]
       family = MONOMIAL
       order = CONSTANT
   []
   [fission_tally_std_dev]
       family = MONOMIAL
       order = CONSTANT
   []
[]


# May want this for special parameters in the paper
# Otherwise
# This AuxVariable and AuxKernel is only here to get the postprocessors
# to evaluate correctly. This can be deleted after MOOSE issue #17534 is fixed.
[AuxVariables]
    [dummy]
    []
  []

[AuxKernels]
  [dummy]
    type = ConstantAux
    variable = dummy
    value = 0.0
  []
[]

# TODO may be useful depending on if MOOSE or OpenMC runs first
# [ICs]
#   [temp]
#       type = ConstantIC
#       variable = temp
#       value = 293.0
#   []
# []

# we probably want the mesh to be the same bewtween MOOSE and OpenMC
# but allow both of them to change due to thermal expansion. TODO look into this
[Problem]
  type = OpenMCCellAverageProblem
  initial_properties = xml
  verbose = true
  tally_type = mesh
  solid_cell_level = 0
  solid_blocks = ANY_BLOCK_ID
  # scaling TODO just do everything in cm since scaling may be messed up rn
  power = P
[]

[Executioner]
  type=Steady
  nl_forced_its = 5
  verbose = true
[]


#TODO determine if we need these. I think probably not?
# [MultiApps]
#     [app]
#         type = TransientMultiApp
#         app_type = CardinalApp
#         input_files = 'sub.i'
#         execute_on = timestep_end
#         sub_cylcing = true
#     []
# []


# [Transfers]
#     [get_temperature_from_sub]
#         type = MultiAppNearestNodeTransfer
#         source_variable = sub_temp
#         variable = temp
#         from_multi_app = app
#     []
#     [send_heat_source_to_sub]
#         type =MultiAppNearestNodeTransfer
#         source_variable = heat_source
#         variable = openmc_heat_source
#         to_multi_app = app
#         from_postprrocessor_to_be_preserved = total_openmc_source
#         to_postprocessor_to_be_preserved = total_received_source
#     []
# []

[Postprocessors]
    [total_openmc_source]
        type = ElementIntegralVariablePostprocessor
        variable = heat_source
    []
[]

[Outputs]
    exodus = true
[]