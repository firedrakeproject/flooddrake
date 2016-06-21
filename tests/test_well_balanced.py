
from __future__ import division  # Get proper divison

import math
import random

import numpy as np


from firedrake import *

from flooddrake import *

# test well balanced condition for none flat bedding


def test_well_balanced():

    n = 5
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
    g = interpolate(Expression(['1', 0, 0]), V)  # uniform depth

    # setup bed
    bed = interpolate(Expression(["pow(x[0]-0.5,2)*4", 0, 0]), V)

    # setup actual depth
    w = g.assign(g - bed)
    
    # source term
    source = interpolate(Expression("0"),X)

    w_start = Function(V).assign(w)

    # timestep
    t_end = 0.01
    solution = Timestepper(V, VCG, bed, source, Courant=0.025)
    w_end = solution.stepper(0, t_end, w)

    h_start, mu_start, mv_start = split(w_start)
    h_end, mu_end, mv_end = split(w_end)

    depth_start = Function(X).project(h_start)
    depth_end = Function(X).project(h_end)

    depth_diff = np.max(np.abs(depth_start.dat.data - depth_end.dat.data))

    assert depth_diff <= 1e-4


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
