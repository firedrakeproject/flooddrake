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

# parameters
parameters["flooddrake"].update({"eps2": 8e-3})

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).assign(0.05)

# setup bed
bed = Function(V)
x = SpatialCoordinate(V.mesh())
bed.sub(0).interpolate(conditional(x[0] > (50 / 50), (1.5 / 50) * (15 / 50), conditional(x[0] < 15 / 50, (10 / 50) * ((15 / 50) - x[0]), conditional(x[0] > (25 / 50), (15 / 50) * (x[0] - (25 / 50)), 0.0))))

# setup state
state = State(V, g, bed)

# setup source (is only a depth function)
source = Function(v_h).assign(0.00025)

# timestep
solution = Timestepper(V, state.bed, source, 0.5)

solution.stepper(0, 100, state.w, 0.5)
