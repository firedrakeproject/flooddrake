from __future__ import division
from __future__ import absolute_import


def Source(source, t, func):
    """ Varies the source term in the DG shallow water implementation over time according to a user defined function

        :param source: The source field
        :type source: :class:`Function'

        :param time: Time
        :type time: double

        :param func: Time varying function to scale source term with. Input is time, Output scalar.
        :type func: str


    """

    return source.assign(func(t)*source)
