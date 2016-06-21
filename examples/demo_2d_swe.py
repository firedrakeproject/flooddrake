
from __future__ import division  # Get proper divison

import math
import random

import numpy as np


from firedrake import *

from flooddrake import *

""" demo file for simple 2d shallow water equations for dam break"""


# Meshsize

n = 10
mesh = UnitSquareMesh(n, n)

# mixed function space
X = FunctionSpace(mesh, "DG", 1)
Y = FunctionSpace(mesh, "DG", 1)
Z = FunctionSpace(mesh, "DG", 1)
V = X * Y * Z


# for slope limiter
XCG = FunctionSpace(mesh, "CG", 1)
YCG = FunctionSpace(mesh, "CG", 1)
ZCG = FunctionSpace(mesh, "CG", 1)
VCG = XCG * YCG * ZCG


# setup free surface depth
g = interpolate(Expression(
    ['pow(x[0]-0.5,2) + pow(x[1]-0.5,2)< 0.05 ? 1 : (pow(x[0]-0.5,2) + pow(x[1]-0.5,2)< 0.05 ? -1.0 : 0.8)', 0, 0]), V)

# setup bed

# pointless trivial second dimension
bed = interpolate(Expression(["0", 0, 0]), V)


# setup actual depth
w = g.assign(g - bed)

# setup source (is only a depth function)
source = interpolate(Expression("0"),X)



# timestep

solution = Timestepper(V, VCG, bed, source, 0.025)

solution.stepper(0, 0.75, w)
