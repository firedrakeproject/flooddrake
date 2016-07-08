from __future__ import division
from __future__ import absolute_import

from firedrake import *


class SlopeLimiter(object):

    def __init__(self, b, V, VCG):

        self.V = V
        self.VCG = VCG
        self.b = b

        if self.V.mesh().geometric_dimension() == 2:

            # split function spaces
            self.v, self.vu, self.vv = split(self.V)
            self.v_cg, self.vu_cg, self.vv_cg = split(self.VCG)

        if self.V.mesh().geometric_dimension() == 1:

            # split function spaces
            self.v, self.vu = split(self.V)
            self.v_cg, self.vu_cg = split(self.VCG)

        # define the cg function to find maximum and minimum values for each vertex
        self.v_max_cg = Function(self.v_cg)
        self.v_min_cg = Function(self.v_cg)

        # define bed and free surface
        self.b_ = Function(self.v).project(self.b)
        self.H = Function(self.v)

        # cell average functions of depth and free surface
        self.c = Function(FunctionSpace(self.v.mesh(), 'DG', 0))
        self.v_func = Function(self.v)
        self.Project = Projector(self.v_func, self.c)

        # Kernel for searching for min and maxes of cell averages
        self.vert_max_min_kernel = """
        #define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
        #define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
        for(int i=0;i<vert_cell_max.dofs;i++){
            vert_cell_max[i][0]=MAX(vert_cell_max[i][0],cell_av[0][0]);
            vert_cell_min[i][0]=MIN(vert_cell_min[i][0],cell_av[0][0]);
        }
        """

        # slope limiting kernel
        self.slope_limiter_kernel = """ float u_min=1000000.0, u_max=-10000000.0;
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
        for(int i=0;i<vert_cell_dg.dofs;i++){
            vert_cell_dg[i][0]=cell_av[0][0] +(alpha*(vert_cell_dg[i][0]-cell_av[0][0]));
        }
        """

        super(SlopeLimiter, self).__init__()

    def Limiter(self, w):
        """ Slope limiter the prevention of shocks. This is from Kuzmin (2011)

            :param w: state vector

            :param b: bed :class:`Function`

            :param VCG: Continuous :class:`MixedFunctionSpace`


        """

        h = w.sub(0)

        self.H.assign(h)
        self.v_func.assign(self.H + self.b_)

        self.Project.project()

        self.v_min_cg.assign(1000000)
        self.v_max_cg.assign(-1000000)

        par_loop(self.vert_max_min_kernel, dx, {
            "cell_av": (self.c, READ),
            "vert_cell_min": (self.v_min_cg, RW),
            "vert_cell_max": (self.v_max_cg, RW)
        })

        par_loop(self.slope_limiter_kernel, dx, {
            "vert_cell_dg": (self.v_func, RW),
            "vert_cell_min": (self.v_min_cg, READ),
            "vert_cell_max": (self.v_max_cg, READ),
            "cell_av": (self.c, READ)
        })

        # limited depth to state vector function
        w.sub(0).assign(self.v_func - self.b_)

        return w
