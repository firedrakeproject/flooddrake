from __future__ import division
from __future__ import absolute_import


class ModelParameters(object):

    """ Model parameters

    """

    def __init__(self):

        self.g = 9.81  # gravity

        super(ModelParameters, self).__init__()

        # could even have a 'epsilon' depth in here, because this is used in
        # different classes!
