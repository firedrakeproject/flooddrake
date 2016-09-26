""" test well balancedness with rive bcs """

from __future__ import division

from firedrake import *
from flooddrake import *


def test_default_boundaries_1d_inflow():

    n = 10
    mesh = UnitIntervalMesh(n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu

    # setup bed
    bed = Function(V)

    # source term
    source = Function(v_h)

    # boundary - only give one boundary
    boundary_w = Function(V)
    marker = 1
    boundary_conditions = [BoundaryConditions(marker, option='inflow', value=boundary_w)]

    # timestep
    solution = Timestepper(V, bed, source, 0.025, boundary_conditions=boundary_conditions)

    # check bcs
    for i in range(2):
        if solution.boundary_conditions[i].marker == marker:
            assert solution.boundary_conditions[i].option == 'inflow'
        else:
            assert solution.boundary_conditions[i].option == 'solid wall'
            assert solution.boundary_conditions[i].value is None


def test_default_boundaries_1d_outflow():

    n = 10
    mesh = UnitIntervalMesh(n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu

    # setup bed
    bed = Function(V)

    # source term
    source = Function(v_h)

    # boundary - only give one boundary
    marker = 1
    boundary_conditions = [BoundaryConditions(marker, option='outflow')]

    # timestep
    solution = Timestepper(V, bed, source, 0.025, boundary_conditions=boundary_conditions)

    # check bcs
    for i in range(2):
        if solution.boundary_conditions[i].marker == marker:
            assert solution.boundary_conditions[i].option == 'outflow'
        else:
            assert solution.boundary_conditions[i].option == 'solid wall'
            assert solution.boundary_conditions[i].value is None


def test_default_boundaries_2d():

    n = 10
    mesh = UnitSquareMesh(n, n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    v_mv = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu * v_mv

    # setup bed
    bed = Function(V)

    # source term
    source = Function(v_h)

    # boundary - only give one boundary
    boundary_w = Function(V)
    marker = 2
    boundary_conditions = [BoundaryConditions(marker, option='inflow', value=boundary_w)]

    # timestep
    solution = Timestepper(V, bed, source, 0.025, boundary_conditions=boundary_conditions)

    # check bcs
    for i in range(4):
        if solution.boundary_conditions[i].marker == marker:
            assert solution.boundary_conditions[i].option == 'inflow'
        else:
            assert solution.boundary_conditions[i].option == 'solid wall'
            assert solution.boundary_conditions[i].value is None


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
