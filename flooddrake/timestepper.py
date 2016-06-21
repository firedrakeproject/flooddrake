
from __future__ import division  # Get proper divison
from __future__ import absolute_import

import math
import random
from flooddrake.slope_modification import SlopeModification
from flooddrake.slope_limiter import SlopeLimiter
from flooddrake.flux import Fluxes
from flooddrake.parameters import ModelParameters
import numpy as np


from firedrake import *

# should really change this child object to the flux, or solver functions,
# but for now just change to new type of class. When changing, change the
# object in class header and the super command.


class Timestepper(object):

    def __init__(self, V, VCG, bed, Courant=0.025):

        self.b = bed
        self.mesh = V.mesh()
        self.VCG = VCG

        # define solver
        self.v = TestFunction(V)

        n = np.power(
            (CellVolume(
                self.mesh).ufl_domain().topology.num_cells() /
                2.0),
            1.0 /
            CellVolume(
                self.mesh).ufl_domain().geometric_dimension())
        
        self.dt = Constant(Courant / n)
        
        self.V = V

        # define flux and state vectors
        if self.mesh.geometric_dimension() == 2:
            self.N = FacetNormal(self.mesh)
            self.b_, _1, _2 = split(self.b)
            self.v_h, self.v_mu, self.mv = split(self.V)
        if self.mesh.geometric_dimension() == 1:
            self.N = FacetNormal(self.mesh)[0]
            self.b_, _1 = split(self.b)
            self.v_h, self.v_mu = split(self.V)

        super(Timestepper, self).__init__()

    def __solver_variable(self, w):
        """ Solves a single timestep of the Firedrake DG model for either the one or two dimensional SWEs

                :param w: State vector function
                :type w: :class:`Function'

        """

        w_ = Function(self.V).assign(0)

        gravity = ModelParameters().g

        if self.mesh.geometric_dimension() == 1:

            # index functions and spaces
            h, mu = split(w)

            # define velocities
            velocity = conditional(h <= 0, zero(mu.ufl_shape), (mu * mu) / h)

            F = as_vector((mu, (velocity + ((gravity / 2) * h * h))))

        if self.mesh.geometric_dimension() == 2:

            # index functions and spaces
            h, mu, mv = split(w)

            # define velocities
            velocity = conditional(h <= 0, zero(mu.ufl_shape), (mu * mv) / h)
            velocity_u = conditional(h <= 0, zero(mu.ufl_shape), (mu * mu) / h)
            velocity_v = conditional(h <= 0, zero(mv.ufl_shape), (mv * mv) / h)

            F1 = as_vector(
                (mu, (velocity_u + ((gravity / 2) * (h * h))), velocity))
            F2 = as_vector(
                (mv, velocity, (velocity_v + ((gravity / 2) * (h * h)))))

        ######################## FLUX MODIFICATION #####################

        # height modification
        h_mod_plus = Max(0, h('+') - Max(0, self.b_('-') - self.b_('+')))
        h_mod_minus = Max(0, h('-') - Max(0, self.b_('+') - self.b_('-')))

        if self.mesh.geometric_dimension() == 2:

            velocity_u = conditional(h <= 0, zero(mu.ufl_shape), (mu) / h)
            velocity_v = conditional(h <= 0, zero(mv.ufl_shape), (mv) / h)

            # modified momentum
            mu_mod_plus = (h_mod_plus) * velocity_u('+')
            mu_mod_minus = (h_mod_minus) * velocity_u('-')
            mv_mod_plus = (h_mod_plus) * velocity_v('+')
            mv_mod_minus = (h_mod_minus) * velocity_v('-')

            # source modification
            delta_plus = as_vector((0, 
            	(gravity/2) * ((h_mod_plus * h_mod_plus) -
            		(h('+') * h('+'))) * self.N[0]('+'), 
            	(gravity / 2) * ((h_mod_plus * h_mod_plus) -
            		(h('+') * h('+'))) * self.N[1]('+')))
            
            delta_minus = as_vector((0, 
            	(gravity / 2) * ((h_mod_minus * h_mod_minus) -
            		(h('-') * h('-'))) * self.N[0]('+'), 
            	(gravity / 2) * ((h_mod_minus * h_mod_minus) -
            		(h('-') * h('-'))) * self.N[1]('+')))

            # set up modified state vectors
            w_plus = as_vector((h_mod_plus, mu_mod_plus, mv_mod_plus))
            w_minus = as_vector((h_mod_minus, mu_mod_minus, mv_mod_minus))

        if self.mesh.geometric_dimension() == 1:

            velocity = conditional(h <= 0, zero(mu.ufl_shape), mu / h)

            # modified momentum
            mu_mod_plus = (h_mod_plus) * velocity('+')
            mu_mod_minus = (h_mod_minus) * velocity('-')

            # source modification
            delta_plus = as_vector(
                (0, (gravity / 2) * ((h_mod_plus * h_mod_plus) - (h('+') * h('+'))) * self.N('+')))
            delta_minus = as_vector(
                (0, (gravity / 2) * ((h_mod_minus * h_mod_minus) - (h('-') * h('-'))) * self.N('+')))

            # set up modified state vectors
            w_plus = as_vector((h_mod_plus, mu_mod_plus))
            w_minus = as_vector((h_mod_minus, mu_mod_minus))

        # Define fluxes
        PosFlux = Fluxes(self.N('+'), self.V)
        pos_flux = PosFlux.Interior_Flux(w_plus, w_minus)
        NegFlux = Fluxes(self.N('-'), self.V)
        neg_flux = NegFlux.Interior_Flux(w_minus, w_plus)
        boundary_flux = PosFlux.Boundary_Flux(w)

        # Define source term
        if self.mesh.geometric_dimension() == 1:
            source = as_vector((0, gravity * h * self.b_.dx(0)))
        if self.mesh.geometric_dimension() == 2:
            source = as_vector(
                (0, gravity * h * self.b_.dx(0), gravity * h * self.b_.dx(1)))

        ####################### SOLVER ######################

        if self.mesh.geometric_dimension() == 2:
            L = (dot(self.v, ((w_ - w) / self.dt)) * dx - 
            	dot(self.v.dx(0), F1) * dx - 
            	dot(self.v.dx(1), F2) * dx + 
            	(dot(self.v('-'), neg_flux) + dot(self.v('+'), pos_flux)) * dS + 
            	dot(self.v, boundary_flux) * ds + 
            	(dot(source, self.v)) * dx - 
            	(dot(self.v('-'), delta_minus) + 
            	dot(self.v('+'), delta_plus)) * dS )
        
        if self.mesh.geometric_dimension() == 1:
            L = (dot(self.v, ((w_ - w) / self.dt)) * dx - 
                dot(self.v.dx(0), F) * dx +
                (dot(self.v('-'),  neg_flux) + dot(self.v('+'), pos_flux)) * dS + 
                dot(self.v, boundary_flux) * ds + 
                (dot(source, self.v)) * dx - 
                (dot(self.v('-'), delta_minus) + 
                dot(self.v('+'), delta_plus)) * dS )

        solve(
            L == 0,
            w_,
            nest=False,
            solver_parameters={
                'ksp_type': 'preonly',
                'pc_type': 'lu'})

        return w_

    def stepper(self, t_start, t_end, w):
        """ Timesteps the shallow water equations from t_start to t_end using a 3rd order SSP Runge-Kutta scheme

                :param t_start: start time
                :type t_start: float

                :param t_end: end time
                :type t_end: float

                :param w: Current state vector function

        """

        # used as the original state vector in each RK3 step
        w_old = Function(self.V).assign(w)

        # raise timestep error
        if (self.dt.dat.data > t_end):
            raise ValueError(
                'end time is less than one timestep or timestep is not factor of end time')

        Nt = int(t_end / self.dt.dat.data)

        ####################### RUNGE-KUTTA SCHEME ####################

        # initial slope modification
        w = SlopeModification(w)

        if self.mesh.geometric_dimension() == 1:

            # index functions and spaces
            h, mu = split(w)

        if self.mesh.geometric_dimension() == 2:

            # index functions and spaces
            h, mu, mv = split(w)


        hout = Function(self.v_h)
        # free surface depth
        hout_file = File("h.pvd")
        # bed depth
        bout = Function(self.v_h).project(self.b_)
        bout_file = File("b.pvd")

        hout_file.write(hout.project(h + self.b_))
        bout_file.write(bout)

        for i in range(Nt):

            A = h
            
            w_ = self.__solver_variable(w)

            w.assign(w_)

            # slope modification
            w = SlopeModification(w)

            w_ = self.__solver_variable(w)

            w.assign((3.0 / 4.0) * w_old + (1.0 / 4.0) * w_)

            # slope modification
            w = SlopeModification(w)

            w_ = self.__solver_variable(w)

            w.assign((1.0 / 3.0) * w_old + (2.0 / 3.0) * w_)

            # slope limiter for last RK step
            w = SlopeLimiter(w, self.b_, self.VCG)
            # slope modification
            w = SlopeModification(w)

            w_old.assign(w)

            # timstep complete

            if self.mesh.geometric_dimension() == 1:

                # index functions and spaces
                h, mu = split(w)

            if self.mesh.geometric_dimension() == 2:

                # index functions and spaces
                h, mu, mv = split(w)

            # This is to check conservation of mass - can delete when tests are
            # in place
            print '---- conservation of mass '
            print np.abs(assemble(h * dx) - assemble(A * dx))
            print '----'
            #

            hout_file.write(hout.project(h + self.b_))
            bout_file.write(bout)

        return w
