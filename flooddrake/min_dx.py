""" finds the minimum edge length of a cell in a mesh """

from flooddrake import *

from mpi4py import MPI


def MinDx(mesh):
    """ Finds the minimum cell edge length for each cell in a DG0 function

        :param mesh: :class:`Mesh` to find min cell edge length of
        :type mesh: :class:`Mesh`

    """

    if mesh.geometric_dimension() == 2:

        min_cell_length = Function(FunctionSpace(mesh, 'DG', 0))
        min_cell_length.interpolate(MinCellEdgeLength(mesh))

    if mesh.geometric_dimension() == 1:

        min_cell_length = Function(FunctionSpace(mesh, 'DG', 0))
        min_cell_length.interpolate(CellVolume(mesh))

    return min_cell_length
