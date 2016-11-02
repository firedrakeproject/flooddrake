""" boundary conditions """

from __future__ import division
from __future__ import absolute_import


# boundary condition options
options = ['solid wall', 'inflow', 'outflow']

# horizontal directions
directions = ['both', 'x', 'y']


class BoundaryConditions(object):

    """ Implementation of a weakly imposed boundary conditions for the boundary flux in flooddrake

        :param marker: the marker of the boundary (e.g. 1, 2, 3, 4) for 2d domain
        :type marker: int

        :param option: boundary condition option
        :type option: str (either 'solid wall', 'inflow' or 'outflow')

        :param value: state vector at marked boundary
        :type value: :class:`Function` (None if option='solid wall')

        See help(mesh) for details on markers. E.g. The boundary markers for a UnitSquareMesh
        are numbered as follows:

            * 1: plane x == 0
            * 2: plane x == 1
            * 3: plane y == 0
            * 4: plane y == 1

    """

    def __init__(self, marker, option='solid wall', value=None, direction='both'):

        self.option = option
        self.value = value
        self.marker = marker
        self.direction = direction

        if self.option not in options:
            raise ValueError('bc option must be either solid wall, inflow or outflow')

        # set any value to None if option is not inflow
        if self.option == 'outflow' or self.option == 'solid wall':
            self.value = None

        if self.option == 'inflow':
            if self.value is None:
                raise ValueError('inflow bc option needs w specified at boundary')

        # check that one of directions is given
        if self.direction not in directions:
            raise ValueError('horizontal direction of condition must either be both, x or y')

        super(BoundaryConditions, self).__init__()
