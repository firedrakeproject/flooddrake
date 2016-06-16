
from __future__ import division  # Get proper divison

import math
import random

import numpy as np


from firedrake import *

from flooddrake import *

# test slope limiter


def test_slope_limiter():

    n = 15
    mesh = PeriodicUnitSquareMesh(n, n)

    # mixed functionspace
    X = FunctionSpace(mesh, "DG", 1)
    Y = FunctionSpace(mesh, "DG", 1)
    Z = FunctionSpace(mesh, "DG", 1)
    V = X * Y * Z

    # mixed functionspace for slope limiter
    XCG = FunctionSpace(mesh, "CG", 1)
    YCG = FunctionSpace(mesh, "CG", 1)
    ZCG = FunctionSpace(mesh, "CG", 1)
    VCG = XCG * YCG * ZCG

    # setup initial condition
    w = interpolate(Expression([1, 1, 1]), V)

    b = interpolate(Expression(['2*x[0]', 0, 0]), V)

    w.assign(w - b)

    b_, b1, b2 = split(b)

    # slope limiting
    W = SlopeLimiter(w, b_, VCG)

    # check that it's invariant as all cell averages are same.
    assert np.max(np.abs(W.dat.data[0] - w.dat.data[0])) < 1e-10


# test slope limiter mean preserving property

def test_slope_limiter_mean_preserving():

    n = 15
    mesh = PeriodicUnitSquareMesh(n, n)

    # mixed functionspace
    X = FunctionSpace(mesh, "DG", 1)
    Y = FunctionSpace(mesh, "DG", 1)
    Z = FunctionSpace(mesh, "DG", 1)
    V = X * Y * Z

    # mixed functionspace for slope limiter
    XCG = FunctionSpace(mesh, "CG", 1)
    YCG = FunctionSpace(mesh, "CG", 1)
    ZCG = FunctionSpace(mesh, "CG", 1)
    VCG = XCG * YCG * ZCG

    # mixed functionspace for average
    XAV = FunctionSpace(mesh, "DG", 0)
    YAV = FunctionSpace(mesh, "DG", 0)
    ZAV = FunctionSpace(mesh, "DG", 0)
    VAV = XAV * YAV * ZAV

    # setup initial condition
    w = interpolate(Expression(['x[0]*x[1]', 0, 0]), V)

    b = interpolate(Expression(['x[0]', 0, 0]), V)

    w.assign(w - b)

    cell_av = Function(VAV).project(w)

    b_, b1, b2 = split(b)

    # slope limiting
    W = SlopeLimiter(w, b_, VCG)

    new_cell_av = Function(VAV).project(w)

    # check that it's mean preserving.
    assert np.max(
        np.abs(
            cell_av.dat.data[0] -
            new_cell_av.dat.data[0])) < 1e-10


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
