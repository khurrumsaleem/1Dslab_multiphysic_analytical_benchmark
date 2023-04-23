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
    file = mesh_1000_in.e
  []
[]

[AuxVariables]
  [temp_analytical]
    family = MONOMIAL
    order = CONSTANT
  []
  [temp_error]
    family = MONOMIAL
    order = CONSTANT
  []
  [dummy_zero]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [temp_error_computer]
    type = ParsedAux
    variable = temp_error
    function = 'temp_analytical - temp'
    args = 'temp_analytical temp'
  []
[]

[Functions]
  [analytical_temp_formula]
    type = ParsedFunction
    vars = 'Sig_t0    T0    L    q    phi0    lam    P'
    vals = '${Sig_t0} ${T0} ${L} ${q} ${phi0} ${lam} ${P}'
    value = 'Sig_t0*T0*sqrt((q*L*phi0/P)*(q*L*phi0/P) - (lam -1)*x*x)'
  []
[]

[ICs]
  [temp]
    type = ConstantIC
    variable = temp
    value = ${T0}
  []
  [kappa_fission_IC]
    type = ConstantIC
    variable = kappa_fission
    value = ${fparse q*eV_to_J/(L*1*1)} # W/cm^3
  []
  [temp_analytical]
    type = FunctionIC
    variable = temp_analytical
    function = analytical_temp_formula
  []
  [dummy_zero]
    type = ConstantIC
    variable = dummy_zero
    value = 0
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  initial_properties = xml
  verbose = true
  tally_type = cell
  tally_blocks = '0'
  tally_score = 'kappa_fission flux'
  source_rate_normalization = 'kappa_fission'
  solid_cell_level = 1
  solid_blocks = ANY_BLOCK_ID
  inactive_batches = 50
  batches = 150
  power = ${fparse P*eV_to_J} # convert from eV/s to W
  relaxation = robbins_monro
  export_properties = true
[]

[Executioner]
  type = Transient
  check_aux = true
  verbose = true
  num_steps = 200
[]

[MultiApps]
  [solid]
      type = TransientMultiApp
      app_type = CardinalApp
      input_files = 'solid_1000.i'
      execute_on = timestep_end
      sub_cycling = false
      interpolate_transfers = false
  []
[]

[Transfers]
  [kappa_fission_to_solid]
    type = MultiAppMeshFunctionTransfer
    to_multi_app = solid
    variable = kappa_fission
    source_variable = kappa_fission
    from_postprocessors_to_be_preserved = kappa_fission
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
  [final]
    type = CSV
    execute_on = 'FINAL'
  []
[]

[Postprocessors]
  [kappa_fission]
    type = ElementIntegralVariablePostprocessor
    variable = kappa_fission
  []
  [L2_distance_analytic_temp_to_field]
    type = ElementL2Error
    variable = temp
    function = analytical_temp_formula
    execute_on = timestep_end
  []
  [L2_norm_temp_analytical]
    type = ElementL2Error
    variable = dummy_zero
    function = analytical_temp_formula
    execute_on = timestep_end
  []
  [ratio_diff_to_analytical_norm]
    type = ParsedPostprocessor
    function = 'L2_distance_analytic_temp_to_field / L2_norm_temp_analytical'
    pp_names = 'L2_distance_analytic_temp_to_field L2_norm_temp_analytical'
    execute_on = timestep_end
  []
[]

[VectorPostprocessors]
  [temp]
    type = ElementValueSampler
    variable = 'temp'
    sort_by = x
    execute_on = 'FINAL'
  []
  [temp_error]
    type = ElementValueSampler
    variable = 'temp_error'
    sort_by = x
    execute_on = 'FINAL'
  []
  [flux]
    type = ElementValueSampler
    variable = 'flux'
    sort_by = x
    execute_on = 'FINAL'
  []
[]