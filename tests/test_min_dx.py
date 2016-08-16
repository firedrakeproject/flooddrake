""" test minimum dx function """

from __future__ import division

import numpy as np

from firedrake import *
from flooddrake import *


def test_min_cell_length_1d():

    n = 10.0
    Dx = 1.0 / n
    mesh = UnitIntervalMesh(n)

    min_cell_length = MinDx(mesh)

    actual_hyp = np.sqrt((Dx ** 2))

    assert np.max(np.abs(min_cell_length.dat.data - actual_hyp)) < 1e-8


def test_min_cell_length_2d():

    n = 10.0
    Dx = 1.0 / n
    mesh = UnitSquareMesh(n, n)

    min_cell_length = MinDx(mesh)

    actual_hyp = np.sqrt((Dx ** 2))

    assert np.max(np.abs(min_cell_length.dat.data - actual_hyp)) < 1e-8


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
