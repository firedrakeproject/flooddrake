""" demo file for simple 1d shallow water equations with steady rain """

from __future__ import division

from firedrake import *
from flooddrake import *

# Meshsize
n = 20
mesh = UnitIntervalMesh(n)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu

# parameters
parameters["flooddrake"].update({"eps1": 1e-4})

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).assign(0.1)

# setup bed
bed = Function(V)
x = SpatialCoordinate(V.mesh())
bed.sub(0).interpolate(conditional(x[0] > 0.75, 2 * (1-x[0]), 2 * abs(x[0]-0.5)))

# setup state
state = State(V, g, bed)

# setup source (is only a depth function)
# constant rain
source = Function(v_h).assign(0.001)

# timestep
solution = Timestepper(V, state.bed, source, 0.25)

solution.stepper(0, 50, state.w, 0.25)
