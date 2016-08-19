from __future__ import division
from __future__ import absolute_import

from firedrake import *


class SlopeModification(object):

    def __init__(self, V):

        self.V = V

        if self.V.mesh().geometric_dimension() == 2:

            # split function spaces
            self.v, self.vu, self.vv = split(self.V)

        if self.V.mesh().geometric_dimension() == 1:

            # split function spaces
            self.v, self.vu = split(self.V)

        self.nf = Function(self.V)

        if self.V.mesh().geometric_dimension() == 2:
            self.new_v_func, self.new_v_u_func, self.new_v_v_func = split(self.nf)

        if self.V.mesh().geometric_dimension() == 1:
            self.new_v_func, self.new_v_u_func = split(self.nf)

        # to remove when user specifies params
        eps1 = parameters["flooddrake"]["eps1"]
        eps2 = parameters["flooddrake"]["eps2"]
        bnd1 = parameters["flooddrake"]["ubnd1"]
        bnd2 = parameters["flooddrake"]["ubnd2"]

        self.slope_modification_2d_kernel = """ double new_cell = 0; const double E=%(epsilon)s;  const double UB=%(ubnd)s; int j;
        for(int i=0;i<vert_cell.dofs;i++){
            new_cell+=vert_cell[i][0];
        }
        new_cell=new_cell/vert_cell.dofs;
        if (new_cell<=0){
            for(int i=0;i<new_vert_cell.dofs;i++){
                new_vert_cell[i][0]=0;
                new_vert_u_cell[i][0]=0;
                new_vert_v_cell[i][0]=0;
            }
        }
        if (new_cell>0){
            int c=0;
            for(int i=0;i<new_vert_cell.dofs;i++){
                if (vert_cell[i][0]>E){
                    new_vert_cell[i][0]=vert_cell[i][0];
                    new_vert_u_cell[i][0]=vert_u_cell[i][0];
                    new_vert_v_cell[i][0]=vert_v_cell[i][0];
                }
                if (vert_cell[i][0]<=E){
                    c=c+1;
                    j=i;
                }
            }
            if (c==0){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    if (sqrt(pow((new_vert_u_cell[i][0]/new_vert_cell[i][0]),2)+pow((new_vert_v_cell[i][0]/new_vert_cell[i][0]),2))>UB){
                        new_vert_u_cell[i][0]=0;
                        new_vert_v_cell[i][0]=0;
                    }
                }
            }
            if (c==1){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    new_vert_cell[i][0]=(new_cell/(new_cell-vert_cell[j][0]))*(vert_cell[i][0]-vert_cell[j][0]);
                    if (new_vert_cell[i][0]>0){
                        if (sqrt(pow((new_vert_u_cell[i][0]/new_vert_cell[i][0]),2)+pow((new_vert_v_cell[i][0]/new_vert_cell[i][0]),2))>UB){
                            new_vert_u_cell[i][0]=0;
                            new_vert_v_cell[i][0]=0;
                        }
                    }
                }
                new_vert_u_cell[j][0]=0;
                new_vert_v_cell[j][0]=0;
            }
            if (c==2){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    if (vert_cell[i][0]<=E){
                        new_vert_cell[i][0]=0;
                        new_vert_u_cell[i][0]=0;
                        new_vert_v_cell[i][0]=0;
                    }
                    if (vert_cell[i][0]>E){
                        new_vert_cell[i][0]=new_cell*vert_cell.dofs;
                        if (sqrt(pow((new_vert_u_cell[i][0]/new_vert_cell[i][0]),2)+pow((new_vert_v_cell[i][0]/new_vert_cell[i][0]),2))>UB){
                            new_vert_u_cell[i][0]=0;
                            new_vert_v_cell[i][0]=0;
                        }
                    }
                }
            }
        }
        """

        self.slope_modification_1d_kernel = """ double new_cell = 0; const double E=%(epsilon)s; const double UB=%(ubnd)s; int j;
        for(int i=0;i<vert_cell.dofs;i++){
            new_cell+=vert_cell[i][0];
        }
        new_cell=new_cell/vert_cell.dofs;
        if (new_cell<=0){
            for(int i=0;i<new_vert_cell.dofs;i++){
                new_vert_cell[i][0]=0;
                new_vert_u_cell[i][0]=0;
            }
        }
        if (new_cell>0){
            int c=0;
            for(int i=0;i<new_vert_cell.dofs;i++){
                if (vert_cell[i][0]>E){
                    new_vert_cell[i][0]=vert_cell[i][0];
                    new_vert_u_cell[i][0]=vert_u_cell[i][0];
                }
                if (vert_cell[i][0]<=E){
                    c=c+1;
                    j=i;
                }
            }
            if (c==0){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    if (new_vert_u_cell[i][0]/new_vert_cell[i][0]>UB){
                        new_vert_u_cell[i][0]=0;
                    }
                    if ((new_vert_u_cell[i][0]/new_vert_cell[i][0])<-UB){
                        new_vert_u_cell[i][0]=0;
                    }
                }
            }
            if (c==1){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    new_vert_cell[i][0]=(new_cell/(new_cell-vert_cell[j][0]))*(vert_cell[i][0]-vert_cell[j][0]);
                    if (new_vert_cell[i][0]>0){
                        if (new_vert_u_cell[i][0]/new_vert_cell[i][0]>UB){
                            new_vert_u_cell[i][0]=0;
                        }
                        if ((new_vert_u_cell[i][0]/new_vert_cell[i][0])<-UB){
                            new_vert_u_cell[i][0]=0;
                        }
                    }
                }
                new_vert_u_cell[j][0]=0;
            }
            if (c==2){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    if (vert_cell[i][0]<=E){
                        new_vert_cell[i][0]=0;
                        new_vert_u_cell[i][0]=0;
                    }
                    if (vert_cell[i][0]>E){
                        new_vert_cell[i][0]=new_cell*vert_cell.dofs;
                        if (new_vert_u_cell[i][0]/new_vert_cell[i][0]>UB){
                            new_vert_u_cell[i][0]=0;
                        }
                        if ((new_vert_u_cell[i][0]/new_vert_cell[i][0])<-UB){
                            new_vert_u_cell[i][0]=0;
                        }
                    }
                }
            }
        }
        """

        # replace parameter strings
        self.slope_modification_1d_kernel = self.slope_modification_1d_kernel % {"epsilon": eps1,
                                                                                 "ubnd": bnd1}
        self.slope_modification_2d_kernel = self.slope_modification_2d_kernel % {"epsilon": eps2,
                                                                                 "ubnd": bnd2}

        super(SlopeModification, self).__init__()

    def Modification(self, w):
        """ Slope modification for the prevention of non negative flows in some verticies of cells. This is from Ern et al (2011)

            :param w: state vector

        """

        if self.V.mesh().geometric_dimension() == 2:

            # split functions
            h, mu, mv = split(w)

        if self.V.mesh().geometric_dimension() == 1:

            # split functions
            h, mu = split(w)

        self.nf.assign(0)

        # par loop

        if self.V.mesh().geometric_dimension() == 2:
            par_loop(self.slope_modification_2d_kernel, dx, {
                "new_vert_v_cell": (self.new_v_v_func, RW),
                "new_vert_u_cell": (self.new_v_u_func, RW),
                "new_vert_cell": (self.new_v_func, RW),
                "vert_cell": (h, READ),
                "vert_u_cell": (mu, READ),
                "vert_v_cell": (mv, READ)
            })

        if self.V.mesh().geometric_dimension() == 1:
            par_loop(self.slope_modification_1d_kernel, dx, {
                "new_vert_u_cell": (self.new_v_u_func, RW),
                "new_vert_cell": (self.new_v_func, RW),
                "vert_cell": (h, READ),
                "vert_u_cell": (mu, READ)
            })

        w.assign(self.nf)

        return w
