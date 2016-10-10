""" demo file for flooding of variable bedding """

from __future__ import division

from firedrake import *
from flooddrake import *


n = 40
mesh = SquareMesh(n, n, 1)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
v_mv = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu * v_mv

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).assign(0.05)

# setup bed
bed = Function(V)
x = SpatialCoordinate(V.mesh())
bed.sub(0).interpolate(conditional(x[0] > (50 / 50), (1.5 / 50) * (15 / 50), conditional(x[0] < 15 / 50, (10 / 50) * ((15 / 50) - x[0]), conditional(x[0] > (25 / 50), (15 / 50) * (x[0] - (25 / 50)), 0.0))))


# setup actual depth
w = g.assign(g - bed)

bedFile = File('bed.pvd')
bedFile.write(bed.sub(0))

# setup source (is only a depth function)
source = Function(v_h).assign(0.00025)

# parameters
parameters["flooddrake"].update({"eps2": 8e-3})

# timestep
solution = Timestepper(V, bed, source, 0.5)

solution.stepper(0, 100, w, 0.5)


# plot actual depth
depthFile = File('depth.pvd')
h = solution.w.sub(0)
depth = Function(v_h).assign(h)
depthFile.write(depth)
