""" test well balancedness """

from __future__ import division

import math
import random
import numpy as np

from firedrake import *
from flooddrake import *


def test_well_balanced():

    n = 5
    mesh = UnitSquareMesh(n, n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    v_mv = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu * v_mv

    # for slope limiter
    v_hcg = FunctionSpace(mesh, "CG", 1)
    v_mucg = FunctionSpace(mesh, "CG", 1)
    v_mvcg = FunctionSpace(mesh, "CG", 1)
    VCG = v_hcg * v_mucg * v_mvcg

    # setup free surface depth
    g = Function(V)
    g.sub(0).assign(1)

    # setup bed
    bed = Function(V)
    x = SpatialCoordinate(V.mesh())
    bed.sub(0).interpolate(4 * pow(x[0] - 0.5, 2))

    # setup actual depth
    w = g.assign(g - bed)

    # source term
    source = Function(v_h)

    w_start = Function(V).assign(w)

    # timestep
    t_end = 0.01
    solution = Timestepper(V, VCG, bed, source, Courant=0.025)
    w_end = solution.stepper(0, t_end, w)

    h_start, mu_start, mv_start = split(w_start)
    h_end, mu_end, mv_end = split(w_end)

    depth_start = Function(v_h).project(h_start)
    depth_end = Function(v_h).project(h_end)

    depth_diff = np.max(np.abs(depth_start.dat.data - depth_end.dat.data))

    assert depth_diff <= 1e-4


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
