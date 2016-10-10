""" demo file for simple 1d shallow water equations with steady rain """

from __future__ import division

from firedrake import *
from flooddrake import *

# Meshsize
n = 50
mesh = UnitIntervalMesh(n)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).assign(0.05)

# setup bed
bed = Function(V)
x = SpatialCoordinate(V.mesh())
bed.sub(0).assign(0)

# setup actual depth
w = g.assign(g - bed)

# setup boundary river - inflow through river, and outflow
boundary_w1 = Function(V)
boundary_w1.sub(0).assign(0.05 - bed.sub(0))
boundary_w1.sub(1).assign(0.05)  # river inflow
boundary_conditions = [BoundaryConditions(1, option='inflow', value=boundary_w1),
                       BoundaryConditions(2, option='outflow')]

# setup source (is only a depth function)
# constant rain
source = Function(v_h)

# timestep
solution = Timestepper(V, bed, source, 0.25, boundary_conditions=boundary_conditions)

solution.stepper(0, 50, w, 0.0625)
