# GLOBAL VARS
# problem physical parameters
T0 = 293
L = 106.47 # equilibrium length from paper (TODO perhaps use formula)
P = 1.0e22 # eV/s
q = 1e8 # eV
k0= 1.25e19 # eV/(s-cm-K^2) k(T) = k0 T(x)
phi0 = 2.5e14 # 1/s-cm^2 flux at the origin
eV_to_J = 1.602e-19 # J per eV
lam = ${fparse 0.5*(1+sqrt(1+(16*q*q*phi0*phi0)/(P*P)))} # eigenvalue solution
h = ${fparse 1/(sqrt(L*(lam-1)/(k0*P)) - (2*T0)/(P)) }

[Mesh]
  [centered_mesh]
    type = FileMeshGenerator
    file = mesh_50_in.e
  []
[]

[Variables]
  [temp]
  []
[]

[ICs]
  [temp]
    type = ConstantIC
    variable = temp
    value = ${T0}
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

# This AuxVariable and AuxKernel is only here to get the postprocessors
# to evaluate correctly. This can be deleted after MOOSE issue #17534 is fixed.
[AuxVariables]
  [heat_source]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = ${fparse q*eV_to_J}
  []
[]

[Functions]
  [conductivity]
    type = ParsedFunction
    value = '${fparse k0*eV_to_J*100} * t' # here 't' means temp. the 100 and ev_to_J is to get units of W/m-K
    # value = ${fparse k0*eV_to_J} TODO see if neg temp is caused by thiisi
  []
[]

[Materials]
  [thermal_parameters]
    type = HeatConductionMaterial
    temp = temp
    thermal_conductivity_temperature_function = conductivity
  []
[]

[BCs]
  # [left_convective_BC]
  #     type = ConvectiveFluxFunction
  #     T_infinity = ${T0}
  #     boundary = left
  #     coefficient = ${h}
  #     variable = temp
  # []
  [test_BC]
    type = DirichletBC
    variable = temp
    value = 310
    boundary = left
  []
  [righ_convective_BC]
      type = ConvectiveFluxFunction
      T_infinity = ${T0}
      boundary = right
      coefficient = ${h}
      variable = temp
  []
[]

[Executioner]
  type = Transient
  nl_abs_tol = 1e-5
  nl_rel_tol = 1e-7
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  verbose = true
[]

[Outputs]
  exodus = true
[]


[Postprocessors]
  [source_integral]
      type = ElementIntegralVariablePostprocessor
      variable = heat_source
      execute_on = transfer
  []
[]