""" test well balancedness """

from __future__ import division

import numpy as np

from firedrake import *
from flooddrake import *


def test_balanced_under_source():

    n = 40
    mesh = UnitIntervalMesh(n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu

    # setup free surface depth
    g = Function(V)
    g.sub(0).assign(0.5)

    # setup bed
    bed = Function(V)

    # setup state
    state = State(V, g, bed)

    # source term
    source = Function(v_h).assign(0.05)

    # timestep
    t_end = 0.01
    solution = Timestepper(V, state.bed, source, 0.025)
    w_end = solution.stepper(0, t_end, state.w, 0.025)

    h_end, mu_end = split(w_end)

    depth_end = Function(v_h).project(h_end)

    # Check max and min depths difference is less than tolerance - balanced
    depth_diff = np.max(depth_end.dat.data) - np.min(depth_end.dat.data)

    assert depth_diff <= 1e-4


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
