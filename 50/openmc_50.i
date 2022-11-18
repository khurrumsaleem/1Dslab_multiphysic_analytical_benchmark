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
  [cell_id]
    family = MONOMIAL
    order = CONSTANT
  []
  [cell_instance]
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
[]

[ICs]
  [temp]
    type = ConstantIC
    variable = temp
    value = ${T0}
  []
  [heat_source]
    type = ConstantIC
    variable = heat_source
    value = ${fparse q*eV_to_J}
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  initial_properties = xml
  verbose = true
  tally_type = mesh
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
  num_steps = 2
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
[]