""" demo file for flooding of variable bedding """

from __future__ import division

from firedrake import *
from flooddrake import *


n = 20
mesh = SquareMesh(n, n, 50)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
v_mv = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu * v_mv

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).assign(0.0)

# setup bed
bed = Function(V)
x = SpatialCoordinate(V.mesh())
bed.sub(0).interpolate(conditional(x[0] < 15, (15 - x[0]), conditional(x[0] > 25, 1.5 * (x[0] - 25), 0.0)))


# setup actual depth
w = g.assign(g - bed)

bedFile = File('bed.pvd')
bedFile.write(bed.sub(0))

# setup source (is only a depth function)
source = Function(v_h).assign(0.075)

# parameters
parameters["flooddrake"].update({"eps2": 1e-12,
                                 "ubnd2": 1e0})

# timestep
solution = Timestepper(V, bed, source, 0.025)

solution.stepper(0, 30, w, 0.025)


# plot actual depth
depthFile = File('depth.pvd')
h = solution.w.sub(0)
depth = Function(v_h).assign(h)
depthFile.write(depth)
