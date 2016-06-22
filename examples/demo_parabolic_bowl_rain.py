""" demo file for simple 2d shallow water equations for parabolic bowl with steady rain """

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

# for slope limiter
v_hcg = FunctionSpace(mesh, "CG", 1)
v_mucg = FunctionSpace(mesh, "CG", 1)
v_mvcg = FunctionSpace(mesh, "CG", 1)
VCG = v_hcg * v_mucg * v_mvcg

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).assign(0.5)

# setup bed
bed = Function(V)
bed.sub(0).interpolate(2 * (pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)))

# setup actual depth
w = g.assign(g - bed)

# setup source (is only a depth function)
source = Function(v_h).assign(0.2)

# timestep
solution = Timestepper(V, VCG, bed, source, float(0.00625 / 2.0))

solution.stepper(0, 2, w)
