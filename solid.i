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

[Variables]
  [temp]
    initial_condition = ${T0}
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

[AuxKernels]
  [dummy]
    type = ConstantAux
    variable = dummy
    value = 0.0
  []
[]


[Functions]
  [conductivity]
    type = ParsedFunction
    value = ${fparse k0 * temp} # TODO unit conversion check
  []
[]

[Materials]
  type = HeatConductionMaterial
  temp = temp
  thermal_conductivity_temperature_function = conductivity
  block = ANY_BLOCK_ID
[]

[BCs]
  [left_convective_BC]
      type = ConvectiveFluxFunction
      T_infinity = T0
      boundary = left
      coefficient = ${fparse h/k0}
      variable =temp
  []
  [righ_convective_BC]
      type = ConvectiveFluxFunction
      T_infinity = T0
      boundary = right
      coefficient = ${fparse h/k0}
      variable =temp
  []
[]

[MultiApps]
  [openmc]
      type = TransientMultiApp
      app_type = CardinalApp
      input_files = 'openmc.i'
      execute_on = timestep_end
  []
[]

[Transfers]
  [heat_source_from_openmc]
    type = MultiAppNearestNodeTransfer
    from_multi_app = openmc
    variable = heat_source
    source_variable = heat_source
    from_postprocessors_to_be_preserved = heat_source
    to_postprocessors_to_be_preserved = source_integral
  []
  [temp_to_openmc]
    type = MultiAppMeshFunctionTransfer
    to_multi_app = openmc
    variable = temp
    source_variable = temp
  []
[]

[Executioner]
  type=Steady
  nl_abs_tol = 1e-8
  nl_forced_its = 5
  verbose = true
[]

[Outputs]
  exodus = true
[]

[Postprocessors]
  [total_openmc_source]
      type = ElementIntegralVariablePostprocessor
      variable = heat_source
  []
[]