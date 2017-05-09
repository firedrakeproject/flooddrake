""" Convergence test case for simple 2d shallow water equations flooding on a slope.
Test case given in Ern et al. 2007 Ritter Solution.
P0 and P1 discontinuous functions considered.
"""

from __future__ import division

from firedrake import *
from flooddrake import *
from firedrake.mg.utils import get_level

import matplotlib.pyplot as plot


def error_of_flood_estimate(mesh, p, T):

    # mixed function space
    h, lvl = get_level(mesh)
    v_h = FunctionSpace(mesh, "DG", p)
    v_mu = FunctionSpace(mesh, "DG", p)
    v_mv = FunctionSpace(mesh, "DG", p)
    V = v_h*v_mu*v_mv

    # parameters
    if p == 0:
        parameters["flooddrake"].update({"eps2": 8e-3})
    if p == 1:
        parameters["flooddrake"].update({"eps2": 3.5e-3})

    # setup free surface depth
    g = Function(V)
    x = SpatialCoordinate(V.mesh())
    g.sub(0).interpolate(conditional(x[0] < 20, 1.0/9.8, parameters["flooddrake"]["eps2"]))

    # setup bed
    bed = Function(V)

    # setup state
    state = State(V, g, bed)

    # setup source (is only a depth function)
    source = Function(v_h)

    # timestep
    solution = Timestepper(V, state.bed, source, 0.05)

    solution.stepper(0, T, state.w, 0.5)

    return solution


# define mesh hierarchy
mesh = RectangleMesh(3, 3, 50, 40)
L = 4
mesh_hierarchy = MeshHierarchy(mesh, L)

# final time
T = 5

# preallocate error
error = np.zeros((2, L))

# find analytic solution
g = parameters["flooddrake"]["gravity"]
xi_0 = 1.0 / 9.8
x = SpatialCoordinate(mesh_hierarchy[-1])
ufl_expression = (conditional(x[0] < 20 - (T * sqrt(g * xi_0)), xi_0,
                  conditional(x[0] > 20 + (2 * T * sqrt(g * xi_0)), 0,
                  (1.0 / (9 * g * (T ** 2))) * pow(x[0] - 20 - (2 * T * sqrt(g * xi_0)), 2))))

# run convergence test over P0 and P1 functions
dx = np.zeros(L)
for p in range(2):

    # get function space of finest solution
    finest_fs = FunctionSpace(mesh_hierarchy[-1], 'DG', p)
    comparison_f = Function(finest_fs)

    # interpolate analytic solution and write to file
    finest_h = Function(finest_fs).interpolate(ufl_expression)
    analytic_hFile = File("analytic_h.pvd")
    analytic_hFile.write(finest_h)

    for i in range(L + 1 - p - 1):

        h = error_of_flood_estimate(mesh_hierarchy[i], p, T)

        if i < L:
            prolong(h.h, comparison_f)
        else:
            comparison_f.assign(h.h)

        error[p, i] = norm(comparison_f - finest_h)
        if p == 0:
            dx[i] = np.max(MinDx(mesh_hierarchy[i]).dat.data)

        print 'completed simulation on level ', i

# rate of decay
A0 = np.vstack([np.linspace(0, L - 1, L), np.ones(L)]).T
s0, res = np.linalg.lstsq(A0, -np.log(error[0, :]) / np.log(2))[0]
A1 = np.vstack([np.linspace(0, L - 2, L - 1), np.ones(L - 1)]).T
s1, res = np.linalg.lstsq(A1, -np.log(error[1, :-1]) / np.log(2))[0]

# plot convergence
plot.loglog(dx, error[0, :], 'r*-', linewidth=3)
plot.loglog(dx[:-1], 8e-1 * error[1, :-1], 'bo-', linewidth=3)
plot.loglog(dx, 2 ** - (s0 * np.linspace(0, L - 1, L)), 'k--', linewidth=3)
plot.loglog(dx[:-1], 2 ** - (s1 * np.linspace(0, L - 2, L - 1)), 'k-', linewidth=3)
plot.legend(['p0 convergence', 'p1 convergence',
             's=' + str(np.round(s0, 2)), 's=' + str(np.round(s1, 2))])
plot.xlabel('dx (m)')
plot.ylabel('norm of error')
plot.axis([np.ceil(dx[-1]), np.ceil(dx[0]), 1e-1, 1e1])
plot.show()
