""" test timestepper """

from __future__ import division

from firedrake import *
from flooddrake import *


def test_timestepper():

    n = 5
    mesh = UnitIntervalMesh(n)

    # mixed function space
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    V = v_h * v_mu

    # for slope limiter
    v_hcg = FunctionSpace(mesh, "CG", 1)
    v_mucg = FunctionSpace(mesh, "CG", 1)
    VCG = v_hcg * v_mucg

    # setup free surface depth
    g = Function(V)
    g.sub(0).assign(0.4)

    # setup bed
    bed = Function(V)
    bed.sub(0).assign(0.1)

    # setup actual depth
    w = g.assign(g - bed)

    # source term
    source = Function(v_h)

    # timestep - make sure it the timestep won't be a factor of final time to
    # see the function deal with it
    t_end = 0.1
    solution = Timestepper(V, VCG, bed, source, Courant=0.015)
    solution.stepper(0, t_end, w, 0.025)

    assert solution.t == t_end


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
