""" demo file for simple 1d shallow water equations with an emerged bump - well balanced.
Courtesy of S. Tomasoni """

from __future__ import division

from firedrake import *
from flooddrake import *


# Meshsize
mesh = IntervalMesh(300, 0, 25)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).assign(0.1)

# setup bed
bed = Function(V)
bed.sub(0).interpolate(conditional((x[0] - 8) * (x[0] - 12) < 0, 0.2 - 0.05 * (x[0] - 10) * (x[0] - 10), 0.0))

# parameters
parameters["flooddrake"].update({"eps1": 1e-20})

# setup state
state = State(V, g, bed)

# setup source (is only a depth function)
source = Function(v_h)

# timestep
solution = Timestepper(V, state.bed, source, 0.5)

solution.stepper(0, 4, state.w, 0.5)
