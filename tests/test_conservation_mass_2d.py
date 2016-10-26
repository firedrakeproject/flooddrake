""" test conservation of mass in 2d dam break problem """

from __future__ import division

import numpy as np

from firedrake import *
from flooddrake import *


def test_conservation_mass_2d_flat_source():

    n = 5
    mesh = UnitSquareMesh(n, n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    v_mv = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu * v_mv

    # setup free surface depth
    g = Function(V)
    g.sub(0).assign(0.8)

    # setup bed
    bed = Function(V)

    # setup state
    state = State(V, g, bed)

    # source term
    source = Function(v_h)

    w_start = Function(V).assign(state.w)

    # timestep
    t_end = 0.01
    solution = Timestepper(V, state.bed, source, 0.025)
    w_end = solution.stepper(0, t_end, state.w, 0.025)

    h_start, mu_start, mv_start = split(w_start)
    h_end, mu_end, mv_end = split(w_end)

    depth_start = Function(v_h).project(h_start)
    depth_end = Function(v_h).project(h_end)

    mass_diff = np.abs(assemble(depth_start * dx) - assemble(depth_end * dx))

    assert mass_diff <= 1e-4


def test_conservation_mass_2d_unflat_source():

    n = 5
    mesh = UnitSquareMesh(n, n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    v_mv = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu * v_mv

    # setup free surface depth
    g = Function(V)
    x = SpatialCoordinate(V.mesh())
    g.sub(0).interpolate(conditional(
        pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2) < 0.05, 1.0, 0.8))

    # setup bed
    bed = Function(V)

    # setup state
    state = State(V, g, bed)

    # source term
    source = Function(v_h)

    w_start = Function(V).assign(state.w)

    # timestep
    t_end = 0.01
    solution = Timestepper(V, state.bed, source, 0.025)
    w_end = solution.stepper(0, t_end, state.w, 0.025)

    h_start, mu_start, mv_start = split(w_start)
    h_end, mu_end, mv_end = split(w_end)

    depth_start = Function(v_h).project(h_start)
    depth_end = Function(v_h).project(h_end)

    mass_diff = np.abs(assemble(depth_start * dx) - assemble(depth_end * dx))

    assert mass_diff <= 1e-4


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
