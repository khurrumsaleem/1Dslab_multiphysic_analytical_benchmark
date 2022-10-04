# GLOBAL VARS
P = 100 # W
T0 = 293 # K sink temp
L0 = 100.0
# L compute from formula or allow expansion from MOOSE to figure it out
# could be cool to do both and see how much they agree
infdim = 50 # cm
rho = 1.2 # g/cc TODO make sure this unit is used
N = 4 # number of mesh regions

  [Mesh]
    [centered_mesh]
      type=GeneratedMeshGenerator
      dim = 3
      xmax = ${fparse -1*L0/2}
      xmin = ${fparse L0/2}
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

# [AuxVariables]
#    # always set
#    [temp]
#        family = MONONIAL
#        order = CONSTANT
#    []
#
#    [heat_source]
#        family = MONOMIAL
#        order = constant
#    []
#    # if you set fluid
#    [density]
#        family = MONOMIAL
#        order = CONSTANT
#    []
#    [fission_tally_std_dev]
#        family = MONOMIAL
#        order = CONSTANT
# []


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