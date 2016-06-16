
from __future__ import division  # Get proper divison

import math
import random

import numpy as np


from firedrake import *

from flooddrake import *

# test slope modification


def test_slope_modification():

    n = 15
    mesh = PeriodicUnitSquareMesh(n, n)

    # mixed functionspace
    X = FunctionSpace(mesh, "DG", 1)
    Y = FunctionSpace(mesh, "DG", 1)
    Z = FunctionSpace(mesh, "DG", 1)
    V = X * Y * Z

    # setup initial condition -> here all -1's -> these should change to zeros
    w = interpolate(Expression([-1, -1, -1]), V)

    # slope modification
    W = SlopeModification(w)
    assert np.max(np.abs(W.dat.data[0])) < 1e-10
    assert np.max(np.abs(W.dat.data[1])) < 1e-10
    assert np.max(np.abs(W.dat.data[2])) < 1e-10

    # now setup a different initial condition -> here everything should stay
    # same within numerical error
    w = interpolate(Expression([1, -1, -1]), V)

    # slope modification
    W = SlopeModification(w)
    assert np.max(np.abs(W.dat.data[0] - 1)) < 1e-10
    assert np.max(np.abs(W.dat.data[1] + 1)) < 1e-10
    assert np.max(np.abs(W.dat.data[2] + 1)) < 1e-10


# test slope modification mean preserving property for nonnegative depth

def test_slope_modification_mean_preserving():

    n = 15
    mesh = PeriodicUnitSquareMesh(n, n)

    # mixed functionspace
    X = FunctionSpace(mesh, "DG", 1)
    Y = FunctionSpace(mesh, "DG", 1)
    Z = FunctionSpace(mesh, "DG", 1)
    V = X * Y * Z

    # mixed functionspace for average
    XAV = FunctionSpace(mesh, "DG", 0)
    YAV = FunctionSpace(mesh, "DG", 0)
    ZAV = FunctionSpace(mesh, "DG", 0)
    VAV = XAV * YAV * ZAV

    # setup initial condition -> here all -1's -> these should change to zeros
    w = interpolate(Expression(['x[0]*x[1]', 0, 0]), V)

    cell_av = Function(VAV).project(w)

    # slope modification
    W = SlopeModification(w)

    new_cell_av = Function(VAV).project(W)

    assert np.max(
        np.abs(
            new_cell_av.dat.data[0] -
            cell_av.dat.data[0])) < 1e-10


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
