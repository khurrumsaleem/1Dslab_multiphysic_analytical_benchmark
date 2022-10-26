# Geometry variables
L = 106.47 # equilibrium length from paper (TODO perhaps use formula)
N = 250 # number of regions in the problem
infdim = 50.0 # length at which the reflective boundary conditions will be to simulate infiniteness in YZ dimension

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