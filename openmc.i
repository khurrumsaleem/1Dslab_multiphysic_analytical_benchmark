# GLOBAL VARS
# problem physical parameters
P = 1.0e22 # eV/s
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

[Problem]
  type = OpenMCCellAverageProblem
  initial_properties = xml
  verbose = true
  tally_type = mesh
  tally_score = kappa_fission
  solid_cell_level = 0
  solid_blocks = ANY_BLOCK_ID
  power = ${fparse P*eV_to_J} # convert from eV/s to W
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