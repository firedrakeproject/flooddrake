""" demo file for simple 2d shallow water equations for flooding of river with random bumpy land"""

from __future__ import division

from firedrake import *
from flooddrake import *

import numpy as np


n = 75
mesh = SquareMesh(n, n, 50)


# setup distribution of rainfall
def rainfall_distribution():
    rv = np.exp(np.random.normal(0, 0.25)) * 0.001
    return rv


# mixed functio n space
v_h = FunctionSpace(mesh, "DG", 1)
v_mu = FunctionSpace(mesh, "DG", 1)
v_mv = FunctionSpace(mesh, "DG", 1)
V = v_h * v_mu * v_mv

# mixed function space
v_hcg = FunctionSpace(mesh, "CG", 1)
v_mucg = FunctionSpace(mesh, "CG", 1)
v_mvcg = FunctionSpace(mesh, "CG", 1)
Vcg = v_hcg * v_mucg * v_mvcg

# setup free surface depth
g = Function(V)
x = SpatialCoordinate(V.mesh())
g.sub(0).assign(0.75)

# setup bed
bed = Function(V)
x = SpatialCoordinate(V.mesh())
H = 0.5
scale = 0.3
e = (1.0 / 50.0) * 2.0
theta = np.pi / 20.0
L = 25.0
bed.sub(0).interpolate(H + ((e * abs(x[0] - 25)) / cos(theta)) +
                       (scale * (sin(4 * np.pi * x[0] / L) * sin(4 * np.pi * x[1] / L))))

# parameters
parameters["flooddrake"].update({"eps2": 8e-2})

w = g.assign(g - bed)

# setup actual depth
w.sub(0).dat.data[np.where(w.sub(0).dat.data[:] <= 0)[0]] = parameters["flooddrake"]["eps2"]

# setup source (is only a depth function)
source = Function(v_h).assign(rainfall_distribution())

# timestep
solution = Timestepper(V, bed, source, 2.0)

solution.stepper(0, 200, w, 2.0)
