""" demo file for simple 1d shallow water equations with steady rain """

from __future__ import division

import math
import random
import numpy as np

from firedrake import *
from flooddrake import *

# Meshsize
n = 15
mesh = UnitIntervalMesh(n)

# mixed function space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu

# for slope limiter
v_hcg = FunctionSpace(mesh, "CG", 1)
v_mucg = FunctionSpace(mesh, "CG", 1)
VCG = v_hcg * v_mucg

# setup free surface depth
g = Function(V)
g.sub(0).assign(0.5)

# setup bed
bed = Function(V)
x = SpatialCoordinate(V.mesh())
bed.sub(0).interpolate(2 * x[0])

# setup actual depth
w = g.assign(g - bed)

# setup source (is only a depth function)
# constant rain
source = Function(v_h).assign(0.05)

# timestep
solution = Timestepper(V, VCG, bed, source, 0.025)

solution.stepper(0, 2, w)
