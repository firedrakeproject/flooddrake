""" demo file for simple 2d shallow water equations for dam break"""

from __future__ import division

from firedrake import *
from flooddrake import *

# Meshsize
n = 10
mesh = UnitSquareMesh(n, n)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
v_mv = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu * v_mv

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).interpolate(conditional(
    pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2) < 0.05, 1.0, 0.8))

# setup bed
bed = Function(V)

# setup state
state = State(V, g, bed)

# setup source (is only a depth function)
source = Function(v_h)

# timestep
solution = Timestepper(V, state.bed, source, 0.1)

solution.stepper(0, 0.75, state.w, 0.1)
