""" test slope modification """

from __future__ import division

import numpy as np

from firedrake import *
from flooddrake import *


def test_slope_modification():

    n = 15
    mesh = PeriodicUnitSquareMesh(n, n)

    # mixed functionspace
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    v_mv = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu * v_mv

    # setup initial condition -> here all -1's -> these should change to zeros
    w = Function(V)
    w.sub(0).assign(-1)
    w.sub(1).assign(-1)
    w.sub(2).assign(-1)

    SM = SlopeModification(V)

    # slope modification
    W = SM.Modification(w)
    assert np.max(np.abs(W.dat.data[0])) < 1e-10
    assert np.max(np.abs(W.dat.data[1])) < 1e-10
    assert np.max(np.abs(W.dat.data[2])) < 1e-10

    # now setup a different initial condition -> here everything should stay
    # same within numerical error
    w = Function(V)
    w.sub(0).assign(1)
    w.sub(1).assign(-0.2)
    w.sub(2).assign(-0.2)

    # slope modification
    W = SM.Modification(w)
    assert np.max(np.abs(W.dat.data[0] - 1)) < 1e-10
    assert np.max(np.abs(W.dat.data[1] + 0.2)) < 1e-10
    assert np.max(np.abs(W.dat.data[2] + 0.2)) < 1e-10


def test_slope_modification_mean_preserving():

    n = 15
    mesh = PeriodicUnitSquareMesh(n, n)

    # mixed functionspace
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    v_mv = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu * v_mv

    # mixed functionspace
    v_hav = FunctionSpace(mesh, "DG", 0)
    v_muav = FunctionSpace(mesh, "DG", 0)
    v_mvav = FunctionSpace(mesh, "DG", 0)
    VAV = v_hav * v_muav * v_mvav

    # setup initial condition
    w = Function(V)
    x = SpatialCoordinate(V.mesh())
    w.sub(0).interpolate(x[0] * x[1])

    cell_av = Function(VAV).project(w)

    SM = SlopeModification(V)

    # slope modification
    W = SM.Modification(w)

    new_cell_av = Function(VAV).project(W)

    assert np.max(
        np.abs(
            new_cell_av.dat.data[0] -
            cell_av.dat.data[0])) < 1e-10


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
