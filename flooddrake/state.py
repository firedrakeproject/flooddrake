""" defines a state object for the state vector and topography """

from __future__ import division

from __future__ import absolute_import

from firedrake import *

import numpy as np


class State(object):

    def __init__(self, V, init_state_vector, bed):

        """ Creates (and initializes) a state object given the state vector and bed topography

            :arg V: The :class:`FunctionSpace` that the topography and state vector live on
            :type V: :class:`FunctionSpace`.

            :arg init_state_vec: The initial state vector (2 components in 1D, 3 in 2D)
            :type init_state_vec: :class:`Function:

            :arg bed: The topography of the problem
            :type bed: :class:`Function`

        """

        # check arguments
        assert isinstance(init_state_vector, function.Function) is True
        assert isinstance(bed, function.Function) is True
        assert bed.function_space() == init_state_vector.function_space()

        # get dimension
        self.dim = V.mesh().geometric_dimension()

        # check if dimensions with Functions agree
        assert len(init_state_vector) == self.dim + 1
        assert len(init_state_vector) == self.dim + 1

        # define free surface height
        self.w = Function(V).assign(init_state_vector - bed)

        # cancel out negative depths (below bed) -> and set to negigible depth (slope mod)
        if self.dim == 1:
            ind = np.where(self.w.sub(0).dat.data[:] <= 0)[0]
            self.w.sub(0).dat.data[ind] = parameters["flooddrake"]["eps1"]
        if self.dim == 2:
            ind = np.where(self.w.sub(0).dat.data[:] <= 0)[0]
            self.w.sub(0).dat.data[ind] = parameters["flooddrake"]["eps2"]

        self.bed = bed

        super(State, self).__init__()
