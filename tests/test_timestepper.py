""" test timestepper """

from __future__ import division

import numpy as np

from firedrake import *
from flooddrake import *


def test_timestepper_1():

    n = 5
    mesh = UnitIntervalMesh(n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu

    # setup free surface depth
    g = Function(V)
    g.sub(0).assign(0.4)

    # setup bed
    bed = Function(V)
    bed.sub(0).assign(0.1)

    # setup state
    state = State(V, g, bed)

    # source term
    source = Function(v_h)

    # timestep - make sure it the timestep won't be a factor of final time to
    # see the function deal with it
    t_end = 0.1
    solution = Timestepper(V, state.bed, source, 0.025)
    solution.stepper(0, t_end, state.w, 0.025)

    assert solution.t == t_end


def test_timestepper_2():

    n = 5
    mesh = UnitIntervalMesh(n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu

    # setup free surface depth
    g = Function(V)
    g.sub(0).assign(0.4)

    # setup bed
    bed = Function(V)
    bed.sub(0).assign(0.1)

    # setup state
    state = State(V, g, bed)

    # source term
    source = Function(v_h)

    # timestep - make sure it the timestep won't be a factor of final time to
    # see the function deal with it
    t_end = 0.1
    solution = Timestepper(V, state.bed, source, 0.025)
    solution.stepper(0, t_end, state.w, 0.025)

    assert solution.t == t_end

    assert solution.c == np.floor(t_end / 0.025) + 1


def test_timestepper_3():

    n = 5
    mesh = UnitIntervalMesh(n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu

    # setup free surface depth
    g = Function(V)
    g.sub(0).assign(0.4)

    # setup bed
    bed = Function(V)
    bed.sub(0).assign(0.1)

    # setup state
    state = State(V, g, bed)

    # source term
    source = Function(v_h)

    # timestep - make sure it the timestep won't be a factor of final time to
    # see the function deal with it
    t_end = 0.1
    solution = Timestepper(V, state.bed, source, 0.025)
    solution.stepper(0, t_end, state.w, 0.015)

    assert solution.t == t_end

    assert solution.c == np.floor(t_end / 0.015) + 1


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
