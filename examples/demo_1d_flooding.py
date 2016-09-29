""" demo file for simple 2d shallow water equations flooding on a slope """

from __future__ import division

from firedrake import *
from flooddrake import *

# Meshsize
n = 40
mesh = IntervalMesh(n, 50)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu

# parameters
parameters["flooddrake"].update({"eps1": 2e-1})

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).interpolate(conditional(x[0] < 20, 10.0/9.8, parameters["flooddrake"]["eps1"]))

# setup bed
bed = Function(V)

# setup actual depth
w = g.assign(g - bed)

# setup source (is only a depth function)
source = Function(v_h)

# timestep
solution = Timestepper(V, bed, source, 0.25)

solution.stepper(0, 10, w, 0.5)
