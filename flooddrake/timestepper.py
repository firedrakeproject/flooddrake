from __future__ import division
from __future__ import absolute_import

from flooddrake.slope_modification import SlopeModification
from flooddrake.slope_limiter import SlopeLimiter
from flooddrake.flux import Interior_Flux, Boundary_Flux
from flooddrake.min_dx import MinDx
from flooddrake.adaptive_timestepping import AdaptiveTimestepping
from flooddrake.boundary_conditions import BoundaryConditions

from pyop2.profiling import timed_stage

from firedrake import *
from firedrake.logging import warning, RED

import numpy as np


# integer markers can be any of these types -> ds() only takes np.int as subdomain_id
marker_int_types = [np.int64, np.int32, np.int16, np.int8, np.int]


class Timestepper(object):

    def __init__(self, V, bed, source=None, MaxTimestep=0.025, func=lambda x: 1,
                 boundary_conditions=None, MinTimestep=1e-8):

        self.b = bed

        # currently source term only works for constant
        self.source_term = source
        self.func = func

        # initialize time attribute
        self.t = 0

        self.mesh = V.mesh()

        # define solver
        self.v = TestFunction(V)

        # Compute global minimum of cell edge length
        self.delta_x = MinDx(self.mesh)

        self.AT = AdaptiveTimestepping(V, MaxTimestep)
        self.dt = MaxTimestep
        self.MinTimestep = MinTimestep
        self.initial_dt = MaxTimestep
        self.Dt = Constant(self.dt)

        self.V = V

        # boundary conditions - default at solid wall for makers not given
        self.boundary_conditions = boundary_conditions
        if self.boundary_conditions is None:
            self.boundary_conditions = []
            for i in self.mesh.topology.exterior_facets.unique_markers:
                if type(i) in marker_int_types:
                    i = int(i)  # convert to np.int
                self.boundary_conditions.append(BoundaryConditions(i))
        else:
            if isinstance(self.boundary_conditions, list) is False:
                raise TypeError('boundary conditions need to be given as an list of conditions')
            markers = np.copy(self.mesh.topology.exterior_facets.unique_markers)
            # make a direction list for each element in markers
            if self.mesh.geometric_dimension() == 2:
                directions = []
                for m in markers:
                    directions.append(['x', 'y'])
            for bc in self.boundary_conditions:
                # check marker is within list of markers (and if all directions are set)
                if bc.marker not in markers:
                    raise ValueError('marker ' +
                                     str(bc.marker) +
                                     ' supplied to BoundaryConditions is ' +
                                     'not valid for mesh')
                if self.mesh.geometric_dimension() == 1:
                    markers = np.delete(markers, np.where(markers == bc.marker)[0])
                if self.mesh.geometric_dimension() == 2:
                    if bc.direction is 'both':
                        directions.pop(np.where(markers == bc.marker)[0])
                        markers = np.delete(markers, np.where(markers == bc.marker)[0])
                    else:
                        if bc.direction is 'x':
                            directions[np.where(markers == bc.marker)[0]].remove('x')
                        if bc.direction is 'y':
                            directions[np.where(markers == bc.marker)[0]].remove('y')
                        if len(directions[np.where(markers == bc.marker)[0]]) == 0:
                            directions.pop(np.where(markers == bc.marker)[0])
                            markers = np.delete(markers, np.where(markers == bc.marker)[0])
            # default markers / directions of markers that are still in list
            for i in markers:
                if type(i) in marker_int_types:
                    ii = int(i)  # convert to np.int
                if self.mesh.geometric_dimension() == 1:
                    self.boundary_conditions.append(BoundaryConditions(ii))
                if self.mesh.geometric_dimension() == 2:
                    if len(directions[np.where(markers == i)[0]]) == 2:
                        self.boundary_conditions.append(BoundaryConditions(ii))
                    if len(directions[np.where(markers == i)[0]]) == 1:
                        # now only do x or y direction default marker
                        self.boundary_conditions.append(BoundaryConditions(ii, direction=directions[np.where(markers == i)[0]][0]))

        # now get all markers
        self.mesh_markers = self.mesh.topology.exterior_facets.unique_markers

        # define flux and state vectors
        if self.mesh.geometric_dimension() == 2:
            self.N = FacetNormal(self.mesh)
            self.b_, _1, _2 = split(self.b)
            self.v_h, self.v_mu, self.mv = self.V.split()
        if self.mesh.geometric_dimension() == 1:
            self.N = FacetNormal(self.mesh)[0]
            self.b_, _1 = split(self.b)
            self.v_h, self.v_mu = self.V.split()

        # define default source
        if self.source_term is None:
            self.source_term = Function(self.v_h)

        self.gravity = parameters["flooddrake"]["gravity"]

        # negigible depth
        if self.mesh.geometric_dimension() == 2:
            self.E = parameters["flooddrake"]["eps2"]
        if self.mesh.geometric_dimension() == 1:
            self.E = parameters["flooddrake"]["eps1"]

        # set plotting scaling negigible depth constant
        self.plot_tol = 1.01

        self.SM = SlopeModification(self.V)
        self.SL = SlopeLimiter(self.b_, self.V)

        super(Timestepper, self).__init__()

    def __update_slope_modification(self):
        """ Updates slope modifier

        """

        self.w.assign(self.SM.Modification(self.w))

    def __update_slope_limiter(self):
        """ Updates slope limiter

        """

        self.w.assign(self.SL.Limiter(self.w))

    def __solver_setup(self):
        """ Sets up the solver

        """

        self.W = TrialFunction(self.V)

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

        # Define Boundary Fluxes - one for each boundary marker
        self.BoundaryFlux = []
        markers = np.copy(self.mesh_markers)
        self.__BCS = []
        for i in range(len(self.mesh_markers)):
            self.__BCS.append([])
        # iterate over all bcs and put them into one array for each marker
        for i in range(len(self.boundary_conditions)):
            self.__BCS[np.where(self.boundary_conditions[i].marker == markers)[0]].append(self.boundary_conditions[i])
        # compute the boundary flux for each marker
        for i in range(len(self.mesh_markers)):
            self.BoundaryFlux.append(Boundary_Flux(self.V, self.w, self.__BCS[i]))
            if self.mesh.geometric_dimension() == 2:
                if len(self.__BCS[i]) == 2 or self.__BCS[i][0].direction == 'both':
                    # check all markers are covered
                    markers = np.delete(markers, np.where(markers == self.__BCS[i][0].marker)[0])
            if self.mesh.geometric_dimension() == 1:
                markers = np.delete(markers, np.where(markers == self.__BCS[i][0].marker)[0])
        # check all boundaries have a flux
        assert len(markers) == 0

        # Define source term - these get overidden every step
        if self.mesh.geometric_dimension() == 1:
            self.source = as_vector((-self.source_term*self.func(self.t), self.gravity * self.h * self.b_.dx(0)))
        if self.mesh.geometric_dimension() == 2:
            self.source = as_vector((-self.source_term*self.func(self.t), self.gravity * self.h * self.b_.dx(0), self.gravity * self.h * self.b_.dx(1)))

        # solver - try to make this only once in the new version - by just assigning different values to all variable functions
        if self.mesh.geometric_dimension() == 2:
            self.L = - dot(self.v, self.W) * dx
            self.a = (- dot(self.v, self.w) * dx +
                      self.Dt * (- dot(self.v.dx(0), self.F1) * dx -
                                 dot(self.v.dx(1), self.F2) * dx +
                                 (dot(self.v('-'), self.NegFlux) +
                                  dot(self.v('+'), self.PosFlux)) * dS +
                                 (dot(self.source, self.v)) * dx -
                                 (dot(self.v('-'), self.delta_minus) +
                                  dot(self.v('+'), self.delta_plus)) * dS))
            # add in boundary fluxes
            for i in range(len(self.mesh_markers)):
                self.a += (self.Dt * (dot(self.v, self.BoundaryFlux[i]) *
                                      ds(self.__BCS[i][0].marker)))

        if self.mesh.geometric_dimension() == 1:
            self.L = - dot(self.v, self.W) * dx
            self.a = (- dot(self.v, self.w) * dx +
                      self.Dt * (- dot(self.v.dx(0), self.F) * dx +
                                 (dot(self.v('-'), self.NegFlux) +
                                  dot(self.v('+'), self.PosFlux)) * dS +
                                 (dot(self.source, self.v)) * dx -
                                 (dot(self.v('-'), self.delta_minus) +
                                  dot(self.v('+'), self.delta_plus)) * dS))
            # add in boundary fluxes
            for i in range(len(self.mesh_markers)):
                self.a += (self.Dt * (dot(self.v, self.BoundaryFlux[i]) *
                                      ds(self.__BCS[i][0].marker)))

        self.w_ = Function(self.V)

        self.problem = LinearVariationalProblem(self.L, self.a, self.w_)
        self.solver = LinearVariationalSolver(self.problem,
                                              solver_parameters={'ksp_type': 'preonly',
                                                                 'sub_pc_type': 'ilu',
                                                                 'pc_type': 'bjacobi',
                                                                 'mat_type': 'aij'})

    def stepper(self, t_start, t_end, w, t_visualization):
        """ Timesteps the shallow water equations from t_start to t_end using a 3rd order SSP Runge-Kutta scheme

                :param t_start: start time
                :type t_start: float

                :param t_end: end time
                :type t_end: float

                :param w: Current state vector function

                :param t_visualization: time interval to write the free surface water depth to a .pvd file
                :type t_visualization: float

        """

        self.w = w

        # setup the solver - this is a timed stage for profiling reasons!
        with timed_stage("Setup of forms and solver"):
            self.__solver_setup()

        # used as the original state vector in each RK3 step
        self.w_old = Function(self.V).assign(self.w)

        # initial slope modification
        self.__update_slope_modification()

        hout = Function(self.v_h)
        hout.rename("free surface depth")
        # free surface depth
        hout_file = File("h.pvd")
        # bed depth
        bout = Function(self.v_h).project(self.b_)
        bout.rename("topography")
        bout_file = File("b.pvd")

        self.Project = Projector(conditional(self.h <= (self.plot_tol * self.E), self.b_, self.h + self.b_), hout)
        self.Project.project()
        hout_file.write(hout)
        bout_file.write(bout)

        self.t = t_start

        # start counter of how many time dumps
        self.c = 1

        # warning boolean marker
        self.wmark = 0

        while self.t < t_end:

            # find new timestep
            with timed_stage("Finding adaptive time-step"):
                self.dt = self.AT.FindTimestep(self.w)
                self.Dt.assign(self.dt)

            # check that prescribed timestep doesn't fall below minimum timestep
            if self.dt < self.MinTimestep:
                if self.wmark is False:
                    warning(RED % "Minimum timestep has been reached. " +
                            "Simulation might become unstable.")
                    self.wmark = 1
                self.dt = self.MinTimestep
                self.Dt.assign(self.dt)

            # check if remaining time to next time dump is less than timestep
            # correct if neeeded
            if self.dt + self.t > self.c * t_visualization:
                self.dt = (self.c * t_visualization) - self.t
                self.Dt.assign(self.dt)

            # check if remaining time to end time is less than timestep
            # correct if needed
            if (self.t <= (self.c + 1) * t_visualization) and (self.dt + self.t > t_end):
                self.dt = t_end - self.t
                self.Dt.assign(self.dt)

            with timed_stage("Runge-Kutta time-stepping scheme"):
                self.solver.solve()

                self.w.assign(self.w_)

                # slope limiter
                with timed_stage("Slope limiting"):
                    self.__update_slope_limiter()
                # slope modification
                with timed_stage("Slope modification"):
                    self.__update_slope_modification()

                self.solver.solve()

                self.w.assign((3.0 / 4.0) * self.w_old + (1.0 / 4.0) * self.w_)

                # slope limiter
                with timed_stage("Slope limiting"):
                    self.__update_slope_limiter()
                # slope modification
                with timed_stage("Slope modification"):
                    self.__update_slope_modification()

                self.solver.solve()

                self.w.assign((1.0 / 3.0) * self.w_old + (2.0 / 3.0) * self.w_)

                # slope limiter
                with timed_stage("Slope limiting"):
                    self.__update_slope_limiter()
                # slope modification
                with timed_stage("Slope modification"):
                    self.__update_slope_modification()

                self.w_old.assign(self.w)

                # timstep complete - dump realisation if needed
                self.t += self.dt

            with timed_stage("Visualization"):
                if self.t == self.c * t_visualization:
                    self.Project.project()
                    hout_file.write(hout)
                    bout_file.write(bout)

                    self.c += 1

                    self.dt = self.initial_dt
                    self.Dt.assign(self.dt)

        # return timestep
        self.dt = self.initial_dt
        self.Dt.assign(self.dt)

        return self.w
