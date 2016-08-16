""" adaptive time stepping for dg shallow water equation model """

from __future__ import division
from __future__ import absolute_import

from mpi4py import MPI

from flooddrake import *
from flooddrake.parameters import ModelParameters
from flooddrake.min_dx import MinDx

import numpy as np


class AdaptiveTimestepping(object):

    def __init__(self, V, max_timestep):

        self.V = V
        self.mesh = self.V.mesh()
        self.N = FacetNormal(self.mesh)

        self.max_timestep = max_timestep
        self.p = self.V.ufl_element().degree()

        self.gravity = ModelParameters().g

        self.c_w_s = Function(FunctionSpace(self.mesh, 'DG', 0))

        # find min cell edge lengths
        self.min_lengths = MinDx(self.V.mesh())

        self.max_wave_speed_kernel_2d = """ const double g=9.8; float max=-10000000, wave_speed=0; int a=0;
        for(int i=0;i<vert_u_cell.dofs;i++){
            if (vert_cell[i][0]<=0){
                wave_speed=-10000000;
                a=a+1;
            }
            if (vert_cell[i][0]>0){
                wave_speed=sqrt(pow((vert_u_cell[i][0]/vert_cell[i][0]),2)+pow((vert_v_cell[i][0]/vert_cell[i][0]),2))+sqrt(g*vert_cell[i][0]);
            }
            max=fmax(wave_speed,max);
        }
        if (a==vert_u_cell.dofs){
            cell_wave_speed[0][0]=10000000;
        }
        if (a<vert_u_cell.dofs){
            cell_wave_speed[0][0]=cell_lengths[0][0]/max;
        }
        """

        self.max_wave_speed_kernel_1d = """ const double g=9.8; float max=-10000000, wave_speed=0; int a=0;
        for(int i=0;i<vert_u_cell.dofs;i++){
            if (vert_cell[i][0]<=0){
                wave_speed=-10000000;
                a=a+1;
            }
            if (vert_cell[i][0]>0){
                wave_speed=fabs(vert_u_cell[i][0]/vert_cell[i][0])+sqrt(g*vert_cell[i][0]);
            }
            max=fmax(wave_speed,max);
        }
        if (a==vert_u_cell.dofs){
            cell_wave_speed[0][0]=10000000;
        }
        if (a<vert_u_cell.dofs){
            cell_wave_speed[0][0]=cell_lengths[0][0]/max;
        }
        """

        super(AdaptiveTimestepping, self).__init__()

    def FindTimestep(self, w):
        """ Finds the CFL criterion satisfying timestep of the DG flooddrake model

            :param w: state vector

        """

        if self.V.mesh().geometric_dimension() == 2:

            # split functions
            h, mu, mv = split(w)

        if self.V.mesh().geometric_dimension() == 1:

            # split functions
            h, mu = split(w)

        # find minimum cfl timestep
        self.c_w_s.assign(0)

        if self.V.mesh().geometric_dimension() == 2:

            par_loop(self.max_wave_speed_kernel_2d, dx, {"cell_wave_speed": (self.c_w_s, RW),
                                                         "vert_cell": (h, READ),
                                                         "vert_u_cell": (mu, READ),
                                                         "vert_v_cell": (mv, READ),
                                                         "cell_lengths": (self.min_lengths, READ)})

        if self.V.mesh().geometric_dimension() == 1:

            par_loop(self.max_wave_speed_kernel_1d, dx, {"cell_wave_speed": (self.c_w_s, RW),
                                                         "vert_cell": (h, READ),
                                                         "vert_u_cell": (mu, READ),
                                                         "cell_lengths": (self.min_lengths, READ)})

        cfl_timestep = self.c_w_s.comm.allreduce(self.c_w_s.dat.data_ro.min(),
                                                 MPI.MIN)

        delta_t = (1.0 / ((2.0 * self.p) + 1)) * cfl_timestep

        return np.min([delta_t * 0.5, self.max_timestep])
