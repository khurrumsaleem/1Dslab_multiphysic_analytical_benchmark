# GLOBAL VARS
# problem physical parameters
T0 = 293
L0 = 100
L = 106.47 # equilibrium length from paper (TODO perhaps use formula)
rho = 1.2 # g/cc
P = 1.0e22 # eV/s
q = 1e8 # eV
k0= 1.25e19 # eV/(s-cm-K^2) k(T) = k0 T(x)
N = 4 # number of regions in the problem
infdim = 50.0 # length at which the reflective boundary conditions will be to simulate infiniteness in YZ dimension
eV_to_J = 1.602e-19 # J per eV
[Mesh]
  [centered_mesh]
    type = FileMeshGenerator
    file = mesh_in.e
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

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = temp
  []
  [heat_source]
    type = CoupledForce
    variable = temp
    v = heat_source
  []
[]

[Functions]
  [conductivity]
    type = ParsedFunction
    value = ${fparse k0 * temp}
  []
[]

[Materials]
  type = HeatConductionMaterial
  prop_names = 'thermal_conductivity'
  temp = temp
  thermal_conductivity_temperature_function = conductivity
  block = ANY_BLOCK_ID
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

[Problem]
  type = OpenMCCellAverageProblem
  initial_properties = xml
  verbose = true
  tally_type = mesh
  solid_cell_level = 0
  solid_blocks = ANY_BLOCK_ID
  # scaling TODO just do everything in cm since scaling may be messed up rn
  power = ${fparse P*eV_to_J} # in W
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