from __future__ import division
from __future__ import absolute_import

from firedrake import *


parameters.add(Parameters("flooddrake",
                          gravity=9.8,
                          eps1=1e-6,
                          eps2=1e-6))
