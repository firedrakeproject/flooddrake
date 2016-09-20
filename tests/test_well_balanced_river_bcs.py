""" test well balancedness with rive bcs """

from __future__ import division

import numpy as np

from firedrake import *
from flooddrake import *


def test_well_balanced():

    n = 40
    mesh = UnitIntervalMesh(n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu

    # setup free surface depth
    g = Function(V)
    x = SpatialCoordinate(V.mesh())
    g.sub(0).assign(0.4)

    # setup bed
    bed = Function(V)
    bed.sub(0).interpolate(conditional(x[0] > 0.75, 2 * (1-x[0]), 2 * abs(x[0]-0.5)))

    # setup actual depth
    w = g.assign(g - bed)

    # source term
    source = Function(v_h)

    w_start = Function(V).assign(w)
    SM = SlopeModification(V)
    ds = SM.Modification(w_start)
    depth_start = Function(v_h).project(ds.sub(0))

    # boundary w
    boundary_w1 = Function(V)
    boundary_w1.sub(0).assign(0.4 - bed.sub(0))
    boundary_conditions = [BoundaryConditions(1, option='inflow', value=boundary_w1),
                           BoundaryConditions(2, option='outflow')]

    # timestep
    t_end = 0.5
    solution = Timestepper(V, bed, source, 0.025, boundary_conditions=boundary_conditions)
    w_end = solution.stepper(0, t_end, w, 0.025)

    h_start, mu_start = split(w_start)
    h_end, mu_end = split(w_end)

    depth_end = Function(v_h).project(h_end)

    depth_diff = np.max(np.abs(depth_start.dat.data - depth_end.dat.data))

    assert depth_diff <= 1e-4


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
