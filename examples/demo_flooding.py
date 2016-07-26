""" demo file for simple 2d shallow water equations flooding on a slope """

from __future__ import division

from firedrake import *
from flooddrake import *

# Meshsize
n = 10
mesh = SquareMesh(n, n, 10)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
v_mv = FunctionSpace(mesh, "DG", 1)
V = v_h*v_mu*v_mv

# for slope limiter
v_hcg = FunctionSpace(mesh, "CG", 1)
v_mucg = FunctionSpace(mesh, "CG", 1)
v_mvcg = FunctionSpace(mesh, "CG", 1)
VCG = v_hcg*v_mucg*v_mvcg

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).interpolate(conditional(x[0] < 4, 1.0/9.8, 0.0))

# setup bed
bed = Function(V)

# setup actual depth
w = g.assign(g - bed)

# setup source (is only a depth function)
source = Function(v_h)

# timestep
solution = Timestepper(V, VCG, bed, source, 0.00625)

solution.stepper(0, 50, w, 0.025)
