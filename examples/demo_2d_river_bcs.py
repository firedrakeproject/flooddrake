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
g.sub(0).assign(0.5)

# setup bed
bed = Function(V)

# setup actual depth
w = g.assign(g - bed)

# setup source (is only a depth function)
source = Function(v_h)

# boundary conditions - river flow through 1, 2 and solid wall 3, 4.
boundary_w1 = Function(V)
boundary_w1.sub(1).assign(0.25)
boundary_w2 = Function(V)
boundary_w2.sub(1).assign(0.25)

boundary_conditions = [BoundaryConditions(1, option='river', value=boundary_w1),
                       BoundaryConditions(2, option='river', value=boundary_w2)]

# timestep
solution = Timestepper(V, bed, source, 0.1, boundary_conditions=boundary_conditions)

solution.stepper(0, 50, w, 0.1)
