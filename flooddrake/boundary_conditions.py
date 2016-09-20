""" boundary conditions """

from __future__ import division
from __future__ import absolute_import


# boundary condition options
options = ['solid wall', 'river']


class BoundaryConditions(object):

    """ Implementation of a weakly imposed boundary conditions for the boundary flux in flooddrake

        :param marker: the marker of the boundary (e.g. 1, 2, 3, 4) for 2d domain
        :type marker: int

        :param option: boundary condition option
        :type option: str (either 'solid wall' or 'river')

        :param value: state vector at marked boundary
        :type value: :class:`Function` (None if option='solid wall')

        See help(mesh) for details on markers. E.g. The boundary markers for a UnitSquareMesh
        are numbered as follows:

            * 1: plane x == 0
            * 2: plane x == 1
            * 3: plane y == 0
            * 4: plane y == 1

    """

    def __init__(self, marker, option='solid wall', value=None):

        self.option = option
        self.value = value
        self.marker = marker

        if self.option not in options:
            raise ValueError('bc option must be either solid wall or river')

        if self.option == 'river':
            if self.value is None:
                raise ValueError('river bc option needs w specified at boundary')

        super(BoundaryConditions, self).__init__()
