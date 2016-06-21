
from __future__ import division  # Get proper divison

import math
import random

import numpy as np


from firedrake import *

from flooddrake import *

""" demo file for simple 1d shallow water equations with steady rain """


# Meshsize

n = 15
mesh = UnitIntervalMesh(n)

# mixed function space
X = FunctionSpace(mesh, "DG", 1)
Y = FunctionSpace(mesh, "DG", 1)
V = X * Y


# for slope limiter
XCG = FunctionSpace(mesh, "CG", 1)
YCG = FunctionSpace(mesh, "CG", 1)
VCG = XCG * YCG 


# setup free surface depth
g = interpolate(Expression(
    ["0.5", 0]), V)

# setup bed

# pointless trivial second dimension
bed = interpolate(Expression(["x[0]*2", 0]), V)


# setup actual depth
w = g.assign(g - bed)

# setup source (is only a depth function)
# constant rain
source = interpolate(Expression("0.05"),X)

# timestep

solution = Timestepper(V, VCG, bed, source, 0.025)

solution.stepper(0, 2, w)
