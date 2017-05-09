from __future__ import division
from __future__ import absolute_import

from firedrake import *


class SlopeModification(object):

    def __init__(self, V):

        self.V = V

        if self.V.mesh().geometric_dimension() == 2:

            # split function spaces
            self.v, self.vu, self.vv = self.V.split()

        if self.V.mesh().geometric_dimension() == 1:

            # split function spaces
            self.v, self.vu = self.V.split()

        self.nf = Function(self.V)

        if self.V.mesh().geometric_dimension() == 2:
            self.new_v_func, self.new_v_u_func, self.new_v_v_func = split(self.nf)

        if self.V.mesh().geometric_dimension() == 1:
            self.new_v_func, self.new_v_u_func = split(self.nf)

        # to remove when user specifies params
        eps1 = parameters["flooddrake"]["eps1"]
        eps2 = parameters["flooddrake"]["eps2"]

        # add error if p>1 - in future make an option where we can project to p=1 from p>1
        if self.V.ufl_element().degree() > 1:
            raise ValueError('Slope Modification technique cannot be used with DGp, p > 1')

        if self.V.ufl_element().degree() == 1:
            self.slope_modification_2d_kernel = """ double new_cell = 0; const double E=%(epsilon)s; int a=0, n1, n2, n3, a1=0, a3=0, flag1, flag2, flag3, deltau, deltav, npos;
            #define STEP(X) (((X) <= 0) ? 0 : 1)
            for(int i=0;i<vert_cell.dofs;i++){
                new_cell+=vert_cell[i][0];
            }
            new_cell=new_cell/vert_cell.dofs;
            for(int i=0;i<new_vert_cell.dofs;i++){
                if (vert_cell[i][0]>E){
                    a=a+1;
                }
            }
            if (a==new_vert_cell.dofs){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    new_vert_cell[i][0]=vert_cell[i][0];
                }
            }
            if (new_cell<=E){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    new_vert_cell[i][0]=new_cell;
                }
            }
            if (new_cell>E){
                if (a<new_vert_cell.dofs){
                    for(int i=1;i<new_vert_cell.dofs;i++){
                        if (vert_cell[0][0]>=vert_cell[i][0]){
                            a1=a1+1;
                        }
                    }
                    for(int i=0;i<new_vert_cell.dofs-1;i++){
                        if (vert_cell[2][0]>=vert_cell[i][0]){
                            a3=a3+1;
                        }
                    }
                    if (a1==2){
                        n3=0;
                        if (a3==2){
                            n2=2;
                            n1=1;
                        }
                        if (a3==1){
                            n2=2;
                            n1=1;
                        }
                        if (a3==0){
                            n1=2;
                            n2=1;
                        }
                    }
                    if (a1==0){
                        n1=0;
                        if (a3==2){
                            n3=2;
                            n2=1;
                        }
                        if (a3==1){
                            n2=2;
                            n3=1;
                        }
                        if (a3==0){
                            n2=2;
                            n3=1;
                        }
                    }
                    if (a1==1){
                        n2=0;
                        if (a3==0){
                            n3=1;
                            n1=2;
                        }
                        if (a3==2){
                            n3=2;
                            n1=1;
                        }
                        if (a3==1){
                            if (vert_cell[0][0]>=vert_cell[2][0]){
                                n3=1;
                                n1=2;
                            }
                            if (vert_cell[0][0]<vert_cell[2][0]){
                                n3=2;
                                n1=1;
                            }
                        }
                    }
                    new_vert_cell[n1][0]=E;
                    new_vert_cell[n2][0]=fmax(E,vert_cell[n2][0]-(new_vert_cell[n1][0]-vert_cell[n1][0])/2.0);
                    new_vert_cell[n3][0]=vert_cell[n3][0]-(new_vert_cell[n1][0]-vert_cell[n1][0])-(new_vert_cell[n2][0]-vert_cell[n2][0]);
                }
            }
            flag1=STEP(new_vert_cell[0][0]-E);
            flag2=STEP(new_vert_cell[1][0]-E);
            flag3=STEP(new_vert_cell[2][0]-E);
            npos=flag1+flag2+flag3;
            deltau=(vert_u_cell[0][0]*(1-flag1))+(vert_u_cell[1][0]*(1-flag2))+(vert_u_cell[2][0]*(1-flag3));
            deltav=(vert_v_cell[0][0]*(1-flag1))+(vert_v_cell[1][0]*(1-flag2))+(vert_v_cell[2][0]*(1-flag3));
            if (npos>0){
                new_vert_u_cell[0][0]=flag1*(vert_u_cell[0][0]+(deltau/npos));
                new_vert_u_cell[1][0]=flag2*(vert_u_cell[1][0]+(deltau/npos));
                new_vert_u_cell[2][0]=flag3*(vert_u_cell[2][0]+(deltau/npos));
                new_vert_v_cell[0][0]=flag1*(vert_v_cell[0][0]+(deltav/npos));
                new_vert_v_cell[1][0]=flag2*(vert_v_cell[1][0]+(deltav/npos));
                new_vert_v_cell[2][0]=flag3*(vert_v_cell[2][0]+(deltav/npos));
            }
            if (npos==0){
                new_vert_u_cell[0][0]=0;
                new_vert_u_cell[1][0]=0;
                new_vert_u_cell[2][0]=0;
                new_vert_v_cell[0][0]=0;
                new_vert_v_cell[1][0]=0;
                new_vert_v_cell[2][0]=0;
            }
            for(int i=0;i<new_vert_cell.dofs;i++){
                if  (new_vert_cell[i][0]<=0){
                    new_vert_cell[i][0]=0;
                }
            }
            """

            self.slope_modification_1d_kernel = """ double new_cell = 0; const double E=%(epsilon)s; int a=0, n1, n2, flag1, flag2, deltau, npos;
            #define STEP(X) (((X) <= 0) ? 0 : 1)
            for(int i=0;i<vert_cell.dofs;i++){
                new_cell+=vert_cell[i][0];
            }
            new_cell=new_cell/vert_cell.dofs;
            for(int i=0;i<new_vert_cell.dofs;i++){
                if (vert_cell[i][0]>E){
                    a=a+1;
                }
            }
            if (a==new_vert_cell.dofs){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    new_vert_cell[i][0]=vert_cell[i][0];
                }
            }
            if (new_cell<=E){
                for(int i=0;i<new_vert_cell.dofs;i++){
                    new_vert_cell[i][0]=new_cell;
                }
            }
            if (new_cell>E){
                if (a<new_vert_cell.dofs){
                    if (vert_cell[0][0]>=vert_cell[1][0]){
                        n1=1;
                        n2=0;
                    }
                    if (vert_cell[0][0]<vert_cell[1][0]){
                        n1=0;
                        n2=1;
                    }
                    new_vert_cell[n1][0]=E;
                    new_vert_cell[n2][0]=vert_cell[n2][0]-(new_vert_cell[n1][0]-vert_cell[n1][0]);
                }
            }
            flag1=STEP(new_vert_cell[0][0]-E);
            flag2=STEP(new_vert_cell[1][0]-E);
            npos=flag1+flag2;
            deltau=(vert_u_cell[0][0]*(1-flag1))+(vert_u_cell[1][0]*(1-flag2));
            if (npos>0){
                new_vert_u_cell[0][0]=flag1*(vert_u_cell[0][0]+(deltau/npos));
                new_vert_u_cell[1][0]=flag2*(vert_u_cell[1][0]+(deltau/npos));
            }
            if (npos==0){
                new_vert_u_cell[0][0]=0;
                new_vert_u_cell[1][0]=0;
            }
            for(int i=0;i<new_vert_cell.dofs;i++){
                if  (new_vert_cell[i][0]<=0){
                    new_vert_cell[i][0]=0;
                }
            }
            """

        if self.V.ufl_element().degree() == 0:
            self.slope_modification_2d_kernel = """ const double E=%(epsilon)s;
            if (vert_cell[0][0]>E){
                new_vert_cell[0][0]=vert_cell[0][0];
                new_vert_u_cell[0][0]=vert_u_cell[0][0];
                new_vert_v_cell[0][0]=vert_v_cell[0][0];
            }
            if (vert_cell[0][0]<=E){
                new_vert_cell[0][0]=E;
                new_vert_u_cell[0][0]=0;
                new_vert_v_cell[0][0]=0;
            }
            """

            self.slope_modification_1d_kernel = """ const double E=%(epsilon)s;
            if (vert_cell[0][0]>E){
                new_vert_cell[0][0]=vert_cell[0][0];
                new_vert_u_cell[0][0]=vert_u_cell[0][0];
            }
            if (vert_cell[0][0]<=E){
                new_vert_cell[0][0]=E;
                new_vert_u_cell[0][0]=0;
            }
            """

        # replace parameter strings
        self.slope_modification_1d_kernel = self.slope_modification_1d_kernel % {"epsilon": eps1}
        self.slope_modification_2d_kernel = self.slope_modification_2d_kernel % {"epsilon": eps2}

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
