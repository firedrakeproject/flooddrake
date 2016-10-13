""" demo file for simple 1d shallow water equations for dam break """

from __future__ import division

from firedrake import *
from flooddrake import *

# Meshsize
n = 10
mesh = UnitIntervalMesh(n)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).interpolate(conditional(pow(x[0] - 0.8, 2) < 0.01, 2, 1))

# setup bed
bed = Function(V)
bed.sub(0).interpolate(2 * x[0])

# setup state
state = State(V, g, bed)

# setup source (is only a depth function)
source = Function(v_h)

# timestep
solution = Timestepper(V, state.bed, source, 0.025)

solution.stepper(0, 0.75, state.w, 0.025)
