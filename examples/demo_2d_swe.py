
from __future__ import division # Get proper divison

import math
import random

import numpy as np 


from firedrake import *
parameters["reorder_meshes"] = False

from flooddrake import *

""" demo file for simple 2d shallow water equations """



# Meshsize

n=10
mesh=UnitSquareMesh(n,n)

# mixed function space
X=FunctionSpace(mesh,"DG",1)
Y=FunctionSpace(mesh,"DG",1)
Z=FunctionSpace(mesh,"DG",1)
V=X*Y*Z


# for slope limiter
XCG=FunctionSpace(mesh,"CG",1)
YCG=FunctionSpace(mesh,"CG",1)
ZCG=FunctionSpace(mesh,"CG",1)
VCG=XCG*YCG*ZCG




# setup free surface depth
g=Function(V)
g.interpolate(Expression(['pow(x[0]-0.5,2) + pow(x[1]-0.5,2)< 0.05 ? 1 : (pow(x[0]-0.5,2) + pow(x[1]-0.5,2)< 0.05 ? -1.0 : 0.8)',0,0]))

# setup bed
bed=Function(V)
bed.interpolate(Expression(["pow(x[0]-0.5,2)*0",0,0])) # pointless trivial second dimension



# setup actual depth
w=g.assign(g-bed)


# timestep

solution = Timestepper(V,VCG,bed,0.025)

solution.stepper(0,0.75,w)





