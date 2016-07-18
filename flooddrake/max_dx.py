""" finds the maximum edge length of a cell in a mesh and puts it into DG0 function """

from flooddrake import *

import numpy as np


def MaxDx(mesh):
    """ Finds the maximum cell edge length for each cell in a DG0 function

        :param mesh: :class:`Mesh` to find max cell edge length of
        :type mesh: :class:`Mesh`

    """

    if mesh.geometric_dimension() == 2:

        max_cell_length = Function(FunctionSpace(mesh, 'DG', 0))
        max_cell_length.interpolate(MaxCellEdgeLength(mesh))

    if mesh.geometric_dimension() == 1:

        max_cell_length = Function(FunctionSpace(mesh, 'DG', 0))
        max_cell_length.assign((np.max(mesh.coordinates.dat.data) - np.min(mesh.coordinates.dat.data)) / CellVolume(mesh).ufl_domain().topology.num_cells())

    return max_cell_length
