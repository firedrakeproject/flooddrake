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

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).interpolate(conditional((x[0] + x[1]) < 20, 5.5, 5))

# setup bed
bed = Function(V)
bed.sub(0).interpolate(conditional((x[0] + x[1]) > (100.0 * 1.6), 3.75 * 3.2, (7.5 / 100.0) * (x[0] + x[1])))

# setup actual depth
w = g.assign(g - bed)

# setup source (is only a depth function)
source = Function(v_h)

# parameters
parameters["flooddrake"].update({"eps2": 3e-1})

# timestep
solution = Timestepper(V, bed, source, 0.5)

solution.stepper(0, 50, w, 0.5)
