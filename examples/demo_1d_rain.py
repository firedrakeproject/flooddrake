""" demo file for simple 1d shallow water equations with steady rain """

from __future__ import division

from firedrake import *
from flooddrake import *

# Meshsize
n = 30
mesh = UnitIntervalMesh(n)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).assign(0.0)

# setup bed
bed = Function(V)
x = SpatialCoordinate(V.mesh())
bed.sub(0).interpolate(conditional(x[0] > 0.75, 2 * (1-x[0]), 2 * abs(x[0]-0.5)))

# setup actual depth
w = g.assign(g - bed)

# parameters
parameters["flooddrake"].update({"eps1": 1e-6,
                                 "ubnd1": 1e2})

# setup source (is only a depth function)
# constant rain
source = Function(v_h).assign(0.05)  # realisatic rainfall - 60 mm h^-1

# timestep
solution = Timestepper(V, bed, source, 0.025)

solution.stepper(0, 10, w, 0.1)
