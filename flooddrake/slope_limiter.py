from __future__ import division
from __future__ import absolute_import

from firedrake import *


class SlopeLimiter(object):

    def __init__(self, b, V):

        self.V = V
        self.b = b

        if self.V.mesh().geometric_dimension() == 2:

            # split function spaces
            self.v, self.vu, self.vv = self.V.split()

        if self.V.mesh().geometric_dimension() == 1:

            # split function spaces
            self.v, self.vu = self.V.split()

        # define functions
        self.H = Function(self.v)
        self.v_func = Function(self.v)
        self.b_ = Function(self.v).project(self.b)

        # check that function spaces match
        if self.v != self.vu:
            raise ValueError('Two function spaces from momementum and depth dont match')

        # define the vertex based limiter objects
        self.SL = VertexBasedLimiter(self.v)

        super(SlopeLimiter, self).__init__()

    def Limiter(self, w):
        """ Slope limiter the prevention of shocks. This is from Kuzmin (2011)

            :param w: state vector

        """

        # Carry out limiting on the free surface depth
        h = w.sub(0)

        self.H.assign(h)
        self.v_func.assign(self.H + self.b_)

        self.SL.apply(self.v_func)
        w.sub(0).assign(self.v_func - self.b_)

        # Carry out limiting on the second component
        mu = w.sub(1)

        self.SL.apply(mu)

        if self.V.mesh().geometric_dimension() == 2:

            # Carry out limiting on the third component
            mv = w.sub(2)

            self.SL.apply(mv)

        return w
