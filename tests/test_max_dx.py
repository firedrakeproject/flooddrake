""" test adaptive timestepping """

from __future__ import division

import numpy as np

from firedrake import *
from flooddrake import *


def test_max_cell_length_1d():

    n = 10.0
    Dx = 1.0 / n
    mesh = UnitIntervalMesh(n)

    max_cell_length = MaxDx(mesh)

    actual_hyp = np.sqrt((Dx ** 2))

    assert np.all(max_cell_length.dat.data == actual_hyp)


def test_max_cell_length_2d():

    n = 10.0
    Dx = 1.0 / n
    mesh = UnitSquareMesh(n, n)

    max_cell_length = MaxDx(mesh)

    actual_hyp = np.sqrt((Dx ** 2) + (Dx ** 2))

    assert np.all((max_cell_length.dat.data - actual_hyp) <= 1e-5)


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
