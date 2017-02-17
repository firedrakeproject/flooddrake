""" Convergence test case for simple 2d shallow water equations flooding on a slope.
Test case given in Ern et al. 2007 Ritter Solution
"""

from __future__ import division

from firedrake import *
from flooddrake import *
from firedrake.mg.utils import get_level

import matplotlib.pyplot as plot


def error_of_flood_estimate(mesh, T):

    # mixed function space
    h, lvl = get_level(mesh)
    v_h = FunctionSpace(mesh, "DG", 1)
    v_mu = FunctionSpace(mesh, "DG", 1)
    v_mv = FunctionSpace(mesh, "DG", 1)
    V = v_h*v_mu*v_mv

    # parameters
    parameters["flooddrake"].update({"eps2": 4e-4 * (2 ** (-lvl))})

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
mesh = RectangleMesh(10, 10, 50, 40)
L = 3
mesh_hierarchy = MeshHierarchy(mesh, L)

T = 5

# get function space of finest solution
finest_fs = FunctionSpace(mesh_hierarchy[-1], 'DG', 1)
comparison_f = Function(finest_fs)

# find analytic solution
g = parameters["flooddrake"]["gravity"]
xi_0 = 1.0 / 9.8
x = SpatialCoordinate(mesh_hierarchy[-1])
ufl_expression = (conditional(x[0] < 20 - (T * sqrt(g * xi_0)), xi_0,
                  conditional(x[0] > 20 + (2 * T * sqrt(g * xi_0)), 0,
                  (1.0 / (9 * g * (T ** 2))) * pow(x[0] - 20 - (2 * T * sqrt(g * xi_0)), 2))))

finest_h = Function(finest_fs).interpolate(ufl_expression)

analytic_hFile = File("analytic_h.pvd")
analytic_hFile.write(finest_h)

# preallocate error
error = np.zeros(L + 1)

# run convergence test
for i in range(L + 1):
    h = error_of_flood_estimate(mesh_hierarchy[i], T)
    if i < L:
        prolong(h.h, comparison_f)
    else:
        comparison_f.assign(h.h)
    error[i] = norm(comparison_f - finest_h)
    print 'completed simulation on level ', i

# rate of decay
A = np.vstack([np.linspace(0, L, L + 1), np.ones(L + 1)]).T
s, res = np.linalg.lstsq(A, -np.log(error) / np.log(2))[0]

# plot convergence
plot.semilogy(np.linspace(0, L, L + 1), error, 'r*-', linewidth=3)
plot.semilogy(np.linspace(0, L, L + 1), 2 ** - (s * np.linspace(0, L, L + 1)), 'k--', linewidth=3)
plot.legend(['convergence', 's=' + str(np.round(s, 2))])
plot.xlabel('level')
plot.ylabel('norm of error')
plot.show()
