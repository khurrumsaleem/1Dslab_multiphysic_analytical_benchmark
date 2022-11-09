# GLOBAL VARS
# problem physical parameters
T0 = 293
P = 1.0e22 # eV/s
q = 1e8 # eV
eV_to_J = 1.602e-19 # J per eV

[Mesh]
  [centered_mesh]
    type = FileMeshGenerator
    file = mesh_50_in.e
  []
[]

[AuxVariables]
  # always set
  [temp]
      family = MONOMIAL
      order = CONSTANT
      initial_condition = ${T0}
  []
  [heat_source]
      family = MONOMIAL
      order = CONSTANT
      initial_condition = ${fparse q*eV_to_J}
  []
  [fission_tally_std_dev]
      family = MONOMIAL
      order = CONSTANT
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  initial_properties = xml
  verbose = true
  tally_type = mesh
  tally_score = kappa_fission
  solid_cell_level = 0
  solid_blocks = ANY_BLOCK_ID
  mesh_template = mesh_50_in.e
  power = ${fparse P*eV_to_J} # convert from eV/s to W
[]

[Executioner]
  type = Transient
  # nl_abs_tol = 1e-6
  # nl_rel_tol = 1e-11
  # num_steps = 5
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
  steady_state_detection = true
  steady_state_tolerance = 1e-3
  # check_aux = true
  verbose = true
[]

[MultiApps]
  [solid]
      type = TransientMultiApp
      app_type = CardinalApp
      input_files = 'solid_50.i'
      execute_on = timestep_end
      sub_cycling = true
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
[]