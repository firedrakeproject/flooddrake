""" demo file for simple 2d shallow water equations flooding on a slope """

from __future__ import division

from firedrake import *
from flooddrake import *

# Meshsize
n = 40
mesh = SquareMesh(n, n, 50)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
v_mv = FunctionSpace(mesh, "DG", 1)
V = v_h*v_mu*v_mv

# parameters
parameters["flooddrake"].update({"eps2": 2e-1})

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).interpolate(conditional(x[0] < 20, 10.0/9.8, parameters["flooddrake"]["eps2"]))

# setup bed
bed = Function(V)

# setup state
state = State(V, g, bed)

# setup source (is only a depth function)
source = Function(v_h)

# timestep
solution = Timestepper(V, state.bed, source, 0.05)

solution.stepper(0, 5, state.w, 0.5)
