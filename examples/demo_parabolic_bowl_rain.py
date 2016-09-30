""" demo file for simple 2d shallow water equations for rainfall into a parabolic bowl"""

from __future__ import division

from firedrake import *
from flooddrake import *

# Meshsize
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
g.sub(0).assign(0.4)

# setup bed
bed = Function(V)
bed.sub(0).interpolate((4.0 / (50.0 ** 2)) * (pow(x[0] - 25, 2) + pow(x[1] - 25, 2)))

# setup actual depth
w = g.assign(g - bed)

# setup source (is only a depth function)
source = Function(v_h).assign(0.005)

# parameters
parameters["flooddrake"].update({"eps2": 1e-1})

# timestep
solution = Timestepper(V, bed, source, 0.5)

solution.stepper(0, 75, w, 1)
