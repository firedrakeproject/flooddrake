
from __future__ import division  # Get proper divison

import math
import random

import numpy as np


from firedrake import *

from flooddrake import *

# test conservation of mass on flat 0 bed with both flat source and
# un-flat source


def test_conservation_mass_2d_flat_source():

    n = 10
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
        ['pow(x[0]-0.5,2) + pow(x[1]-0.5,2)< 0.05 ? 0.8 : (pow(x[0]-0.5,2) + pow(x[1]-0.5,2)< 0.05 ? -1.0 : 0.8)', 0, 0]), V)

    # setup bed
    bed = interpolate(Expression(["0", 0, 0]), V)

    # setup actual depth
    w = g.assign(g - bed)

    w_start = Function(V).assign(w)

    # timestep
    t_end = 0.01
    solution = Timestepper(V, VCG, bed, Courant=0.025)
    w_end = solution.stepper(0, t_end, w)

    h_start, mu_start, mv_start = split(w_start)
    h_end, mu_end, mv_end = split(w_end)

    depth_start = Function(X).project(h_start)
    depth_end = Function(X).project(h_end)

    mass_diff = np.abs(assemble(depth_start * dx) - assemble(depth_end * dx))

    assert mass_diff <= 1e-5


def test_conservation_mass_2d_unflat_source():

    n = 10
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
        ['pow(x[0]-0.5,2) + pow(x[1]-0.5,2)< 0.05 ? 1 : (pow(x[0]-0.5,2) + pow(x[1]-0.5,2)< 0.05 ? -1.0 : 0.8)', 0, 0]), V)

    # setup bed
    bed = interpolate(Expression(["0", 0, 0]), V)

    # setup actual depth
    w = g.assign(g - bed)

    w_start = Function(V).assign(w)

    # timestep
    t_end = 0.01
    solution = Timestepper(V, VCG, bed, Courant=0.025)
    w_end = solution.stepper(0, t_end, w)

    h_start, mu_start, mv_start = split(w_start)
    h_end, mu_end, mv_end = split(w_end)

    depth_start = Function(X).project(h_start)
    depth_end = Function(X).project(h_end)

    mass_diff = np.abs(assemble(depth_start * dx) - assemble(depth_end * dx))

    assert mass_diff <= 1e-5


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
