
from __future__ import division # Get proper divison

import math
import random

import numpy as np 


from firedrake import *
parameters["reorder_meshes"] = False

from flooddrake import *

# test slope modification


def test_slope_modification():
    
    n=15
    mesh=PeriodicUnitSquareMesh(n,n)
    
    # mixed functionspace
    X=FunctionSpace(mesh,"DG",1)
    Y=FunctionSpace(mesh,"DG",1)
    Z=FunctionSpace(mesh,"DG",1)
    V=X*Y*Z
    
    # setup initial condition -> here all -1's -> these should change to zeros
    w=Function(V)
    w.interpolate(Expression([-1,-1,-1]))
    
    # slope modification
    W=SlopeModification(w,V)
    if np.all(W.dat.data[0]==0)!=1:
        raise AssertionError('failed')
    if np.all(W.dat.data[1]==0)!=1:
        raise AssertionError('failed')
    if np.all(W.dat.data[2]==0)!=1:
        raise AssertionError('failed')
    
    # now setup a different initial condition -> here everything should stay same within numerical error
    w.interpolate(Expression([1,-1,-1]))
    
    # slope modification
    W=SlopeModification(w,V)
    if np.max(np.abs(W.dat.data[0]-1))>1e-10:
        raise AssertionError('failed')
    if np.max(np.abs(W.dat.data[1]+1))>1e-10:
        raise AssertionError('failed')
    if np.max(np.abs(W.dat.data[2]+1))>1e-10:
        raise AssertionError('failed')



if __name__ == "__main__":
    test_slope_modification()






