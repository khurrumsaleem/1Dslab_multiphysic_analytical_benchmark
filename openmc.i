
# GLOBAL VARS
P = 100
[Mesh]
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

# we probably want the mesh to be the same bewtween MOOSE and OpenMC
# but allow both of them to change due to thermal expansion. TODO look into this
[Problem]
type = OpenMCCellAverageProblem
initial_properties = xml
verbose = true
tally_type = mesh
mesh_template = 'openmc_in.e'
solid_cell_level = 0
solid_blocks = ''
# scaling TODO just do everything in cm since scaling may be messed up rn
power = P
[]

# maybe more steps for iteration for flux/XS
[Executioner]
    type=Transient
    num_steps = 1
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