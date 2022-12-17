# GLOBAL VARS
# problem physical parameters
T0 = 293
P = 1.0e22 # eV/s
q = 1e8 # eV
L = 106.47 # equilibrium length
k0= 1.25e19 # eV/(s-cm-K^2) k(T) = k0 T(x)
phi0 = 2.5e14 # 1/s-cm^2 flux at the origin
eV_to_J = 1.602e-19 # J per eV
lam = ${fparse 0.5*(1+sqrt(1+(16*q*q*phi0*phi0)/(P*P)))} # eigenvalue solution
Sig_t0 = ${fparse sqrt(P/((lam-1)*k0*L))/(T0)}

[Mesh]
  [centered_mesh]
    type = FileMeshGenerator
    file = mesh_50_in.e
  []
[]

[AuxVariables]
  [cell_id]
    family = MONOMIAL
    order = CONSTANT
  []
  [cell_instance]
    family = MONOMIAL
    order = CONSTANT
  []
  [temp_error]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [cell_id]
    type = CellIDAux
    variable = cell_id
  []
  [cell_instance]
    type = CellInstanceAux
    variable = cell_instance
  []
  [analytical temp]
    type = FunctionAux
    variable = analytical_temp
    function = analytical_temp
  []
  [temp_error_computer]
    type = FunctionAux
    variable = temp_error
    function = temp_error_formula
  []
[]

[Functions]
  [analytical_temp]
    type = ParsedFunction
    value = '${fparse Sig_t0*L*sqrt( (q*phi0*L/P)*(*phi0*L/P) - (lam -1)*x*x)}'
  []
  [temp_error]
    type = ParsedFunction
    value = ${fparse analytical_temp-temp}$
  []
[]

[ICs]
  [temp]
    type = ConstantIC
    variable = temp
    value = ${T0}
  []
  [heat_source_IC]
    type = ConstantIC
    variable = heat_source
    value = ${fparse q*eV_to_J/(L*1*1)} # W/cm^3
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  initial_properties = xml
  verbose = true
  tally_type = mesh
  tally_name = heat_source
  tally_score = kappa_fission
  solid_cell_level = 1
  solid_blocks = ANY_BLOCK_ID
  mesh_template = mesh_50_in.e
  power = ${fparse P*eV_to_J} # convert from eV/s to W
[]

[Executioner]
  type = Transient
  check_aux = true
  verbose = true
  num_steps = 4
[]

[MultiApps]
  [solid]
      type = TransientMultiApp
      app_type = CardinalApp
      input_files = 'solid_50.i'
      execute_on = timestep_end
      sub_cycling = false
      interpolate_transfers = false
  []
[]

[Transfers]
  [heat_source_to_solid]
    type = MultiAppMeshFunctionTransfer
    to_multi_app = solid
    variable = heat_source
    source_variable = heat_source
    from_postprocessors_to_be_preserved = heat_source
    to_postprocessors_to_be_preserved = source_integral
  []
  [temp_from_solid]
    type = MultiAppMeshFunctionTransfer
    from_multi_app = solid
    variable = temp
    source_variable = temp
  []
[]

[Outputs]
  exodus = true
  csv = true
[]

[Postprocessors]
  [heat_source]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
  []
  [max_tally_rel_err]
    type = FissionTallyRelativeError
  []
  [max_heat_source]
    type = ElementExtremeValue
    variable = heat_source
  []
[]

[VectorPostprocessors]
  [temp]
    type = ElementValueSampler
    variable = 'temp'
    sort_by = x
    execute_on = timestep_end
  []
  [temp_error]
    variable = 'temp_error'
    sort_by = x
    execute_on = timestep_end
  []
[]