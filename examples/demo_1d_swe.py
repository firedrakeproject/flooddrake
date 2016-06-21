
from __future__ import division  # Get proper divison

import math
import random

import numpy as np


from firedrake import *

from flooddrake import *

""" demo file for simple 2d shallow water equations """


# Meshsize

n = 10
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
    ['pow(x[0]-0.5,2)< 0.02 ? 1 : (pow(x[0]-0.5,2) < 0.02 ? -1.0 : 0.8)', 0]), V)

# setup bed

# pointless trivial second dimension
bed = interpolate(Expression(["x[0]*2", 0]), V)


# setup actual depth
w = g.assign(g - bed)


# timestep

solution = Timestepper(V, VCG, bed, 0.025)

solution.stepper(0, 0.75, w)
