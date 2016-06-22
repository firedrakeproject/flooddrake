from __future__ import division
from __future__ import absolute_import

from firedrake import *


def SlopeLimiter(w, b, VCG):
    """ Slope limiter the prevention of shocks. This is from Kuzmin (2011)

        :param w: state vector

        :param b: bed :class:`Function`

        :param VCG: Continuous :class:`MixedFunctionSpace`


    """

    V = w.function_space()

    if V.mesh().geometric_dimension() == 2:

        # split functions
        h, mu, mv = split(w)

        # split function spaces
        v, vu, vv = split(V)
        v_cg, vu_cg, vv_cg = split(VCG)

    if V.mesh().geometric_dimension() == 1:

        # split functions
        h, mu = split(w)

        # split function spaces
        v, vu = split(V)
        v_cg, vu_cg = split(VCG)

    # make a function for bedding
    b_ = Function(v).project(b)

    # get free surface depth function
    v_func = Function(v).project(h + b_)

    W = Function(V).assign(w)

    # define the cg function to find maximum and minimum values for each vertex
    v_max_cg = Function(v_cg).assign(-1000000)
    v_min_cg = Function(v_cg).assign(1000000)

    # cell average functions of depth and free surface
    c = Function(FunctionSpace(v.mesh(), 'DG', 0)).project(v_func)

    # Kernel for searching for min and maxes of cell averages
    vert_max_min_kernel = """
        #define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
        #define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
        for(int i=0;i<vert_cell_max.dofs;i++){
        vert_cell_max[i][0]=MAX(vert_cell_max[i][0],cell_av[0][0]);
        vert_cell_min[i][0]=MIN(vert_cell_min[i][0],cell_av[0][0]);
        }
        """

    # slope limiting kernel
    slope_limiter_kernel = """ float u_min=1000000.0, u_max=-10000000.0, EPSILON=0;
        #define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
        #define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
        for(int i=0;i<vert_cell_min.dofs;i++){
        u_min=MIN(vert_cell_min[i][0],u_min);
        u_max=MAX(vert_cell_max[i][0],u_max);
        }
        float alpha=1.0;
        for(int j=0;j<vert_cell_dg.dofs;j++){
        if(vert_cell_dg[j][0]>cell_av[0][0]){
        alpha=MIN(MIN(1,(u_max-cell_av[0][0])/(vert_cell_dg[j][0]-cell_av[0][0])),alpha);
        }
        if(vert_cell_dg[j][0]==cell_av[0][0]){
        alpha=MIN(1,alpha);
        }
        if(vert_cell_dg[j][0]<cell_av[0][0]){
        alpha=MIN(MIN(1,(u_min-cell_av[0][0])/(vert_cell_dg[j][0]-cell_av[0][0])),alpha);
        }
        }
        float a=0;
        for(int i=0;i<vert_cell_dg.dofs;i++){
        if (cell_depth[i][0]>EPSILON){
        a=a+1;
        }
        }
        for(int i=0;i<vert_cell_dg.dofs;i++){
        if (a==vert_cell_dg.dofs){
        vert_cell_dg[i][0]=cell_av[0][0] +(alpha*(vert_cell_dg[i][0]-cell_av[0][0]));
        }
        }
        """

    par_loop(vert_max_min_kernel, dx, {
        "cell_av": (c, READ),
        "vert_cell_min": (v_min_cg, RW),
        "vert_cell_max": (v_max_cg, RW)
    })

    par_loop(slope_limiter_kernel, dx, {
        "vert_cell_dg": (v_func, RW),
        "vert_cell_min": (v_min_cg, READ),
        "vert_cell_max": (v_max_cg, READ),
        "cell_av": (c, READ),
        "cell_depth": (h, READ)
    })

    # limited depth to state vector function
    W.sub(0).assign(v_func - b_)

    return W
