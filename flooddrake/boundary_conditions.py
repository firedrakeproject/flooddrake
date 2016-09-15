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

    """

    def __init__(self, marker, option='solid wall', value=None):

        self.option = option
        self.value = value
        self.marker = marker

        if type(marker) is not int:
            raise TypeError('marker must be integer denoting boundary numbers')

        if self.option not in options:
            raise ValueError('bc option must be either solid wall or river')

        if self.option == 'river':
            if self.value is None:
                raise ValueError('river bc option needs w specified at boundary')

        if option == 'solid wall':
            self._list = None

        super(BoundaryConditions, self).__init__()
