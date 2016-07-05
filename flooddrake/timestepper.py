from __future__ import division
from __future__ import absolute_import

from flooddrake.slope_modification import SlopeModification
from flooddrake.slope_limiter import SlopeLimiter
from flooddrake.flux import Interior_Flux, Boundary_Flux
from flooddrake.parameters import ModelParameters
import numpy as np

from firedrake import *

# should really change this child object to the flux, or solver functions,
# but for now just change to new type of class. When changing, change the
# object in class header and the super command.


class Timestepper(object):

    def __init__(self, V, VCG, bed, source, Courant=0.025, func=lambda x: 1):

        self.b = bed

        # currently source term only works for constant
        self.source_term = source
        self.func = func

        self.t = 0

        self.mesh = V.mesh()
        self.VCG = VCG

        # define solver
        self.v = TestFunction(V)

        n = np.power((CellVolume(self.mesh).ufl_domain().topology.num_cells()/2.0),
                     1.0/CellVolume(self.mesh).ufl_domain().geometric_dimension())

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

        self.gravity = ModelParameters().g

        super(Timestepper, self).__init__()

    def __update_slope_modification(self):
        """ Updates slope modifier

        """

        S = SlopeModification(self.w)

        self.w.assign(S)

    def __update_slope_limiter(self):
        """ Updates slope limiter

        """

        S = SlopeLimiter(self.w, self.b_, self.VCG)

        self.w.assign(S)

    def __solver_setup(self):
        """ Sets up the solver

        """

        self.w_ = Function(self.V)

        if self.mesh.geometric_dimension() == 1:

            # index functions and spaces
            self.h, self.mu = self.w.sub(0), self.w.sub(1)

            # define velocities
            velocity = conditional(self.h <= 0, zero(self.mu.ufl_shape), (self.mu * self.mu) / self.h)

            self.F = as_vector((self.mu, (velocity + ((self.gravity / 2) * self.h * self.h))))

        if self.mesh.geometric_dimension() == 2:

            # index functions and spaces
            self.h, self.mu, self.mv = self.w.sub(0), self.w.sub(1), self.w.sub(2)

            # define velocities
            velocity = conditional(self.h <= 0, zero(self.mu.ufl_shape), (self.mu * self.mv) / self.h)
            velocity_u = conditional(self.h <= 0, zero(self.mu.ufl_shape), (self.mu * self.mu) / self.h)
            velocity_v = conditional(self.h <= 0, zero(self.mv.ufl_shape), (self.mv * self.mv) / self.h)

            self.F1 = as_vector((self.mu, (velocity_u + ((self.gravity / 2) * (self.h * self.h))), velocity))
            self.F2 = as_vector((self.mv, velocity, (velocity_v + ((self.gravity / 2) * (self.h * self.h)))))

        # height modification
        h_mod_plus = Max(0, self.h('+') - Max(0, self.b_('-') - self.b_('+')))
        h_mod_minus = Max(0, self.h('-') - Max(0, self.b_('+') - self.b_('-')))

        if self.mesh.geometric_dimension() == 2:

            velocity_u = conditional(self.h <= 0, zero(self.mu.ufl_shape), (self.mu) / self.h)
            velocity_v = conditional(self.h <= 0, zero(self.mv.ufl_shape), (self.mv) / self.h)

            # modified momentum
            mu_mod_plus = (h_mod_plus) * velocity_u('+')
            mu_mod_minus = (h_mod_minus) * velocity_u('-')
            mv_mod_plus = (h_mod_plus) * velocity_v('+')
            mv_mod_minus = (h_mod_minus) * velocity_v('-')

            # source modification
            self.delta_plus = as_vector((0,
                                        (self.gravity / 2) * ((h_mod_plus * h_mod_plus)-(self.h('+') * self.h('+'))) * self.N[0]('+'),
                                        (self.gravity / 2) * ((h_mod_plus * h_mod_plus)-(self.h('+') * self.h('+'))) * self.N[1]('+')))

            self.delta_minus = as_vector((0,
                                         (self.gravity / 2) * ((h_mod_minus * h_mod_minus)-(self.h('-') * self.h('-'))) * self.N[0]('+'),
                                         (self.gravity / 2) * ((h_mod_minus * h_mod_minus)-(self.h('-') * self.h('-'))) * self.N[1]('+')))

            # set up modified state vectors
            self.w_plus = as_vector((h_mod_plus, mu_mod_plus, mv_mod_plus))
            self.w_minus = as_vector((h_mod_minus, mu_mod_minus, mv_mod_minus))

        if self.mesh.geometric_dimension() == 1:

            velocity = conditional(self.h <= 0, zero(self.mu.ufl_shape), self.mu / self.h)

            # modified momentum
            mu_mod_plus = (h_mod_plus) * velocity('+')
            mu_mod_minus = (h_mod_minus) * velocity('-')

            # source modification
            self.delta_plus = as_vector((0,
                                        (self.gravity / 2) * ((h_mod_plus * h_mod_plus) - (self.h('+') * self.h('+'))) * self.N('+')))
            self.delta_minus = as_vector((0,
                                         (self.gravity / 2) * ((h_mod_minus * h_mod_minus) - (self.h('-') * self.h('-'))) * self.N('+')))

            # set up modified state vectors
            self.w_plus = as_vector((h_mod_plus, mu_mod_plus))
            self.w_minus = as_vector((h_mod_minus, mu_mod_minus))

        # Define Fluxes - these get overidden every step
        self.PosFlux = Interior_Flux(self.N('+'), self.V, self.w_plus, self.w_minus)
        self.NegFlux = Interior_Flux(self.N('-'), self.V, self.w_minus, self.w_plus)
        self.BoundaryFlux = Boundary_Flux(self.V, self.w)

        # Define source term - these get overidden every step
        if self.mesh.geometric_dimension() == 1:
            self.source = as_vector((-self.source_term*self.func(self.t), self.gravity * self.h * self.b_.dx(0)))
        if self.mesh.geometric_dimension() == 2:
            self.source = as_vector((-self.source_term*self.func(self.t), self.gravity * self.h * self.b_.dx(0), self.gravity * self.h * self.b_.dx(1)))

        # solver - try to make this only once in the new version - by just assigning different values to all variable functions
        if self.mesh.geometric_dimension() == 2:
            self.L = (dot(self.v, ((self.w_ - self.w) / self.dt)) * dx -
                      dot(self.v.dx(0), self.F1) * dx -
                      dot(self.v.dx(1), self.F2) * dx +
                      (dot(self.v('-'), self.NegFlux) + dot(self.v('+'), self.PosFlux)) * dS +
                      dot(self.v, self.BoundaryFlux) * ds +
                      (dot(self.source, self.v)) * dx -
                      (dot(self.v('-'), self.delta_minus) +
                      dot(self.v('+'), self.delta_plus)) * dS)

        if self.mesh.geometric_dimension() == 1:
            self.L = (dot(self.v, ((self.w_ - self.w) / self.dt)) * dx -
                      dot(self.v.dx(0), self.F) * dx +
                      (dot(self.v('-'), self.NegFlux) + dot(self.v('+'), self.PosFlux)) * dS +
                      dot(self.v, self.BoundaryFlux) * ds +
                      (dot(self.source, self.v)) * dx -
                      (dot(self.v('-'), self.delta_minus) +
                      dot(self.v('+'), self.delta_plus)) * dS)

    def stepper(self, t_start, t_end, w):
        """ Timesteps the shallow water equations from t_start to t_end using a 3rd order SSP Runge-Kutta scheme

                :param t_start: start time
                :type t_start: float

                :param t_end: end time
                :type t_end: float

                :param w: Current state vector function

        """

        self.w = w

        # setup the solver
        self.__solver_setup()

        # used as the original state vector in each RK3 step
        self.w_old = Function(self.V).assign(self.w)

        # raise timestep error
        if (self.dt.dat.data > t_end):
            raise ValueError('end time is less than one timestep or timestep is not factor of end time')

        Nt = int(t_end / self.dt.dat.data)

        # initial slope modification
        self.__update_slope_modification()

        hout = Function(self.v_h)
        # free surface depth
        hout_file = File("h.pvd")
        # bed depth
        bout = Function(self.v_h).project(self.b_)
        bout_file = File("b.pvd")

        hout_file.write(hout.project(self.h + self.b_))
        bout_file.write(bout)

        self.t = 0

        for i in range(Nt):

            solve(self.L == 0, self.w_, nest=False,
                  solver_parameters={'ksp_type': 'preonly',
                                     'pc_type': 'lu'})

            self.w.assign(self.w_)

            # slope modification
            self.__update_slope_modification()

            solve(self.L == 0, self.w_, nest=False,
                  solver_parameters={'ksp_type': 'preonly',
                                     'pc_type': 'lu'})

            self.w.assign((3.0 / 4.0) * self.w_old + (1.0 / 4.0) * self.w_)

            # slope modification
            self.__update_slope_modification()

            solve(self.L == 0, self.w_, nest=False,
                  solver_parameters={'ksp_type': 'preonly',
                                     'pc_type': 'lu'})

            self.w.assign((1.0 / 3.0) * self.w_old + (2.0 / 3.0) * self.w_)

            # slope limiter for last RK step
            self.__update_slope_limiter()
            # slope modification
            self.__update_slope_modification()

            self.w_old.assign(self.w)

            # timstep complete

            hout_file.write(hout.project(self.h + self.b_))
            bout_file.write(bout)

            self.t += self.dt.dat.data

        return self.w
