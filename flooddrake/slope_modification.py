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

        self.slope_modification_2d_kernel = """ double new_cell = 0; const double E=1e-6;  const double UB=1e0; int j;
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
                        if (new_vert_u_cell[i][0]<0){
                            new_vert_u_cell[i][0]=0;
                        }
                        if (new_vert_u_cell[i][0]>=0){
                            new_vert_u_cell[i][0]=0;
                        }
                        if (new_vert_v_cell[i][0]<0){
                            new_vert_v_cell[i][0]=0;
                        }
                        if (new_vert_v_cell[i][0]>=0){
                            new_vert_v_cell[i][0]=0;
                        }
                    }
                }
            }
            if (c==1){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    new_vert_cell[i][0]=(new_cell/(new_cell-vert_cell[j][0]))*(vert_cell[i][0]-vert_cell[j][0]);
                    if (new_vert_cell[i][0]>0){
                        if (sqrt(pow((new_vert_u_cell[i][0]/new_vert_cell[i][0]),2)+pow((new_vert_v_cell[i][0]/new_vert_cell[i][0]),2))>UB){
                            if (new_vert_u_cell[i][0]<0){
                                new_vert_u_cell[i][0]=0;
                            }
                            if (new_vert_u_cell[i][0]>=0){
                                new_vert_u_cell[i][0]=0;
                            }
                            if (new_vert_v_cell[i][0]<0){
                                new_vert_v_cell[i][0]=0;
                            }
                            if (new_vert_v_cell[i][0]>=0){
                                new_vert_v_cell[i][0]=0;
                            }
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
                        new_vert_cell[i][0]=new_cell*3;
                        if (sqrt(pow((new_vert_u_cell[i][0]/new_vert_cell[i][0]),2)+pow((new_vert_v_cell[i][0]/new_vert_cell[i][0]),2))>UB){
                            if (new_vert_u_cell[i][0]<0){
                                new_vert_u_cell[i][0]=0;
                            }
                            if (new_vert_u_cell[i][0]>=0){
                                new_vert_u_cell[i][0]=0;
                            }
                            if (new_vert_v_cell[i][0]<0){
                                new_vert_v_cell[i][0]=0;
                            }
                            if (new_vert_v_cell[i][0]>=0){
                                new_vert_v_cell[i][0]=0;
                            }
                        }
                    }
                }
            }
        }
        """

        self.slope_modification_1d_kernel = """ double new_cell = 0; const double E=1e-6; const double UB=1e0; int j;
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
                        new_vert_u_cell[i][0]=UB*new_vert_cell[i][0];
                    }
                    if ((new_vert_u_cell[i][0]/new_vert_cell[i][0])<-UB){
                        new_vert_u_cell[i][0]=(-UB)*new_vert_cell[i][0];
                    }
                }
            }
            if (c==1){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    new_vert_cell[i][0]=(new_cell/(new_cell-vert_cell[j][0]))*(vert_cell[i][0]-vert_cell[j][0]);
                    if (new_vert_cell[i][0]>0){
                        if (new_vert_u_cell[i][0]/new_vert_cell[i][0]>UB){
                            new_vert_u_cell[i][0]=UB*new_vert_cell[i][0];
                        }
                        if ((new_vert_u_cell[i][0]/new_vert_cell[i][0])<-UB){
                            new_vert_u_cell[i][0]=(-UB)*new_vert_cell[i][0];
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
                        new_vert_cell[i][0]=new_cell*3;
                        if (new_vert_u_cell[i][0]/new_vert_cell[i][0]>UB){
                            new_vert_u_cell[i][0]=UB*new_vert_cell[i][0];
                        }
                        if ((new_vert_u_cell[i][0]/new_vert_cell[i][0])<-UB){
                            new_vert_u_cell[i][0]=(-UB)*new_vert_cell[i][0];
                        }
                    }
                }
            }
        }
        """

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
