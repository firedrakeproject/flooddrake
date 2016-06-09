
from __future__ import division # Get proper divison

import math
import random
import numpy as np 



from firedrake import *
parameters["reorder_meshes"] = False


def SlopeModification(w,V):
    
    """ Slope modification for the prevention of non negative flows in some verticies of cells. This is from Ern et al (2011)
    
    	:param w: State vector function
    	
    	:param V: :class:`MixedFunctionSpace'
    
    """
    
    
    # make input functions
    h, mu, mv = split(w)
    v, vu, vv = split(V)
    v_func=Function(v).project(h)
    v_u_func=Function(vu).project(mu)
    v_v_func=Function(vv).project(mv)
    
    # now define kernel
    slope_modification_kernel = """ float scale=0, E=0, EPSILON=0;
    for(int i=0;i<vert_cell.dofs;i++){
    new_cell[0][0]+=vert_cell[i][0];
    scale=scale+1.0;
    }
    new_cell[0][0]=new_cell[0][0]/scale;
    if (new_cell[0][0]<=E){
    for(int i=0;i<new_vert_cell.dofs;i++){
    new_vert_cell[i][0]= EPSILON;
    new_vert_u_cell[i][0]=0;
    new_vert_v_cell[i][0]=0;
    }
    }
    if (new_cell[0][0]>E){
    for(int i=0;i<new_vert_cell.dofs;i++){
    if (vert_cell[i][0]>E){
    new_vert_cell[i][0]=vert_cell[i][0];
    new_vert_u_cell[i][0]=vert_u_cell[i][0];
    new_vert_v_cell[i][0]=vert_v_cell[i][0];
    }
    if (vert_cell[i][0]<=E){
    new_vert_cell[i][0]=EPSILON;
    new_vert_u_cell[i][0]=0;
    new_vert_v_cell[i][0]=0;
    }
    }
    }
    """
    
    new_v_u_func=Function(vu).assign(0)
    new_v_v_func=Function(vv).assign(0)
    new_v_func=Function(v).assign(0)
    new_func=Function(FunctionSpace(v.mesh(),'DG',0))
    
    # par loop
    par_loop(slope_modification_kernel,dx,{"new_vert_v_cell":(new_v_v_func,RW),"new_vert_u_cell":(new_v_u_func,RW),"new_vert_cell":(new_v_func,RW),"new_cell":(new_func,RW),"vert_cell":(v_func,READ),"vert_u_cell":(v_u_func,READ),"vert_v_cell":(v_v_func,READ)})
    
    # project back to w
    W=Function(w.function_space()).assign(0)
    W.dat.data[0][:]+=new_v_func.dat.data[:].astype(float)
    W.dat.data[1][:]+=new_v_u_func.dat.data[:].astype(float)
    W.dat.data[2][:]+=new_v_v_func.dat.data[:].astype(float)
    
    return W


