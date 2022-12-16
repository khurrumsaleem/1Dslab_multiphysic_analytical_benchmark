# Geometry variables
# NOTE THE GEOMETRY IS IN CM TO MATCH OPENMC AND THUS
# ALL PHYSICAL PARAMETERS NEED TO USE CM FOR LENGTH
L = 106.47 # cm equilibrium length from paper
N = 50 # number of regions in the problem
infdim = 0.5 # cm need y=z=1 length for y,z contribution factor of 1 to the power integral

[Mesh]
  [centered_mesh]
    type=GeneratedMeshGenerator
    dim = 3
    xmax = ${fparse L/2}
    xmin = ${fparse -L/2}
    ymax = ${fparse infdim}
    ymin = ${fparse -infdim}
    zmax = ${fparse infdim}
    zmin = ${fparse -infdim}
    nx = ${fparse N}
    ny = 1
    nz = 1
    length_unit = 'cm'
  []
[]