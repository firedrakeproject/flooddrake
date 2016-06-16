
from __future__ import division  # Get proper divison

import math
import random

import numpy as np


from firedrake import *

from flooddrake import *

""" demo file for simple 2d shallow water equations - from Ern et al 2011"""


# Meshsize

n = 8
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
    ['x[0] < 0.2 ? (3.0/9.8) : (x[0] < 0.2 ? -1.0 : 0.0)', 0, 0]), V)

# setup bed

bed = interpolate(Expression(["0", 0, 0]), V)


# setup actual depth
w = g.assign(g - bed)


# timestep

solution = Timestepper(V, VCG, bed, 0.0125)

solution.stepper(0, 2, w)
