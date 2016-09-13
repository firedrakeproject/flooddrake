from __future__ import division
from __future__ import absolute_import

from firedrake import *


def Interior_Flux(N, V, wr, wl):
    """ Calculates the Interior fluxes between the positively restricted state vector function wr and negatively restricted state vector function wl

            :param wr: Positive / negative restriction of the state vector function
            :type wr: Restricted(:class:`Function')

            :param wl: Negative / positive restriction of the state vector function
            :type wl: Restricted(:class:`Function')

    """

    d = V.mesh().geometric_dimension()

    gravity = parameters["flooddrake"]["gravity"]

    # two different fluxes depending on dimension
    if d == 1:

        hr, mur = wr[0], wr[1]
        hl, mul = wl[0], wl[1]

        # Do HLLE flux
        vr = conditional(hr <= 0, zero(mur.ufl_shape), (mur / hr) * N)
        vl = conditional(hl <= 0, zero(mul.ufl_shape), (mul / hl) * N)

        h_star = 0.5 * (hl + hr)
        y = ((sqrt(hr) * vr) + (sqrt(hl) * vl)) / (sqrt(hl) + sqrt(hr))
        v_star = conditional(eq(0, (sqrt(hl) + sqrt(hr))), zero(y.ufl_shape), y)

        c_plus = Max(0, Max(vl + sqrt(gravity * hl), v_star + sqrt(gravity * h_star)))
        c_minus = Min(0, Min(vr - sqrt(gravity * hr), v_star - sqrt(gravity * h_star)))

        # Define F
        velocityr = conditional(hr <= 0, zero(mur.ufl_shape), (mur * mur) / hr)
        velocityl = conditional(hl <= 0, zero(mul.ufl_shape), (mul * mul) / hl)

        Fr = as_vector((mur, velocityr + (gravity / 2 * ((hr * hr)))))
        Fl = as_vector((mul, velocityl + (gravity / 2 * ((hl * hl)))))
        Wr = as_vector((hr, mur))
        Wl = as_vector((hl, mul))

        # Define Q
        I = as_matrix(((1, 0), (0, 1)))
        C = as_matrix(((0, 1), ((-v_star * v_star) + (gravity * h_star), (2 * v_star))))

        yl = ((c_plus + c_minus) / (c_plus - c_minus))
        yr = ((c_plus * c_minus) / (c_plus - c_minus))

        Q = (conditional(eq(c_plus - c_minus, 0), zero(yl.ufl_shape), yl) * C) - \
            (2 * conditional(eq(c_plus - c_minus, 0), zero(yr.ufl_shape), yr) * I)

        Flux = (0.5 * N * (Fr + Fl)) + (0.5 * dot(Q, Wr - Wl))

        return Flux

    if d == 2:

        hr, mur, mvr = wr[0], wr[1], wr[2]
        hl, mul, mvl = wl[0], wl[1], wl[2]

        # Do HLLC flux
        hl_zero = conditional(hl <= 0, 0, 1)
        ur = conditional(hr <= 0, zero(as_vector((mur / hr, mvr / hr)).ufl_shape),
                         hl_zero * as_vector((mur / hr, mvr / hr)))

        hr_zero = conditional(hr <= 0, 0, 1)
        ul = conditional(hl <= 0, zero(as_vector((mul / hl, mvl / hl)).ufl_shape),
                         hr_zero * as_vector((mul / hl, mvl / hl)))

        vr = dot(ur, N)
        vl = dot(ul, N)

        # set flux constants depending on wavelength
        c_minus = Min(vr - sqrt(gravity * hr), vl - sqrt(gravity * hl))
        c_plus = Min(vr + sqrt(gravity * hr), vl + sqrt(gravity * hl))

        # make sure we don't divide by 0 height
        y = (((0.5 * gravity * hr * hr) - (0.5 * gravity * hl * hl) +
             (hl * vl * (c_plus - vl)) - (hr * vr * (c_minus - vr))) /
             ((hl * (c_plus - vl)) - (hr * (c_minus - vr))))
        c_s = conditional(eq((hr * (c_minus - vr)), (hl * (c_plus - vl))), zero(y.ufl_shape), y)

        velocityl = conditional(hl <= 0, zero(mul.ufl_shape), (hr_zero * mul * mvl) / hl)
        velocity_ul = conditional(hl <= 0, zero(mul.ufl_shape), (hr_zero * mul * mul) / hl)
        velocity_vl = conditional(hl <= 0, zero(mvl.ufl_shape), (hr_zero * mvl * mvl) / hl)
        velocityr = conditional(hr <= 0, zero(mur.ufl_shape), (hl_zero * mur * mvr) / hr)
        velocity_ur = conditional(hr <= 0, zero(mur.ufl_shape), (hl_zero * mur * mur) / hr)
        velocity_vr = conditional(hr <= 0, zero(mvr.ufl_shape), (hl_zero * mvr * mvr) / hr)

        F1r = as_vector((mur,
                         velocity_ur + ((gravity / 2) * (hr * hr)),
                         velocityr))
        F2r = as_vector((mvr,
                         velocityr,
                         velocity_vr + ((gravity / 2) * (hr * hr))))

        F1l = as_vector((mul,
                         velocity_ul + ((gravity / 2) * (hl * hl)),
                         velocityl))
        F2l = as_vector((mvl,
                         velocityl,
                         velocity_vl + ((gravity / 2) * (hl * hl))))

        F_plus = as_vector((F1r, F2r))
        F_minus = as_vector((F1l, F2l))

        W_plus = as_vector((hr, mur, mvr))
        W_minus = as_vector((hl, mul, mvl))

        # conditional to prevent dividing by zero
        y = ((c_minus - vr) / (c_minus - c_s)) * (W_plus -
                                                  as_vector((0,
                                                            hr * (c_s - vl) * N[0],
                                                            hr * (c_s - vl) * N[1])))
        w_plus = conditional(eq(c_minus, c_s), zero(y.ufl_shape), y)

        # conditional to prevent dividing by zero
        y = ((c_plus - vl) / (c_plus - c_s)) * (W_minus -
                                                as_vector((0,
                                                          hl * (c_s - vr) * N[0],
                                                          hl * (c_s - vr) * N[1])))
        w_minus = conditional(eq(c_plus, c_s), zero(y.ufl_shape), y)

        Flux = ((0.5 * dot(N, F_plus + F_minus)) +
                (0.5 * (-((abs(c_minus) - abs(c_s)) * w_minus) +
                        ((abs(c_plus) - abs(c_s)) * w_plus) +
                        (abs(c_minus) * W_plus) -
                        (abs(c_plus) * W_minus))))

        return Flux


def Boundary_Flux(V, w):
    """ Calculates the boundary flux between the state vector and a solid reflective wall (zero velocity, same depth (-> to improve in the future add other boundary conditions options)

            :param w: State vector function
            :type w: :class:`Function:

            :param N: Normal
            :type N: :class:`FacetNormal'


    """

    d = V.mesh().geometric_dimension()

    gravity = parameters["flooddrake"]["gravity"]

    N = FacetNormal(V.mesh())

    # two different fluxes depending on dimension - both solid wall weak boundary conditions
    # has different normal (no restirction!)
    if d == 1:

        h, mu = split(w)
        mur = mu
        mul = Constant(0)  # reflecitve in wall - investigate why this it

        # Do HLLE flux
        vr = conditional(h <= 0, zero(mur.ufl_shape), (mur / h) * N[0])
        vl = conditional(h <= 0, zero(mul.ufl_shape), (mul / h) * N[0])

        h_star = h
        y = ((sqrt(h) * vr) + (sqrt(h) * vl)) / (sqrt(h) + sqrt(h))
        v_star = conditional(h <= 0, zero(y.ufl_shape), y)

        c_plus = Max(0, Max(vl + sqrt(gravity * h), v_star + sqrt(gravity * h_star)))
        c_minus = Min(0, Min(vr - sqrt(gravity * h), v_star - sqrt(gravity * h_star)))

        # Define F
        velocityr = conditional(h <= 0, zero(mur.ufl_shape), (mur * mur) / h)
        velocityl = conditional(h <= 0, zero(mul.ufl_shape), (mul * mul) / h)
        Fr = as_vector((mur, velocityr + (gravity / 2 * ((h * h)))))
        Fl = as_vector((mul, velocityl + (gravity / 2 * ((h * h)))))
        Wr = as_vector((h, mur))
        Wl = as_vector((h, mul))

        # Define Q
        I = as_matrix(((1, 0), (0, 1)))
        C = as_matrix(((0, 1), ((-v_star * v_star) + (gravity * h_star), (2 * v_star))))

        yl = ((c_plus + c_minus) / (c_plus - c_minus))
        yr = ((c_plus * c_minus) / (c_plus - c_minus))

        Q = (conditional(eq(c_plus - c_minus, 0), zero(yl.ufl_shape), yl) * C) - \
            (2 * conditional(eq(c_plus - c_minus, 0), zero(yr.ufl_shape), yr) * I)

        Flux = (0.5 * N[0] * (Fr + Fl)) + (0.5 * dot(Q, Wr - Wl))

        return Flux

    if d == 2:

        h, mu, mv = split(w)

        mul = Constant(0)
        mur = mu  # NO VELOCITY AT WALLS! BUT HEIGHT!

        mvr = mv
        mvl = Constant(0)

        # maybe fit a version where one can specify if one wants mv to be
        # opposite sign - would do this via adding two different variables
        # for mv

        # Do HLLC flux
        ul = conditional(h <= 0, zero(as_vector((mul / h, mvl / h)).ufl_shape),
                         as_vector((mul / h, mvl / h)))

        ur = conditional(h <= 0, zero(as_vector((mur / h, mvr / h)).ufl_shape),
                         as_vector((mur / h, mvr / h)))

        vr = dot(ur, N)
        vl = dot(ul, N)

        # set flux constants depending on wavelength
        c_minus = Min(vr - sqrt(gravity * h), vl - sqrt(gravity * h))
        c_plus = Min(vr + sqrt(gravity * h), vl + sqrt(gravity * h))

        # make sure we don't divide by 0 height
        y = (((0.5 * gravity * h * h) - (0.5 * gravity * h * h) +
             (h * vl * (c_plus - vl)) - (h * vr * (c_minus - vr))) /
             ((h * (c_plus - vl)) - (h * (c_minus - vr))))
        c_s = conditional(eq((h * (c_minus - vr)), (h * (c_plus - vl))), zero(y.ufl_shape), y)

        velocityl = conditional(h <= 0, zero(mul.ufl_shape), (mul * mvl) / h)
        velocity_ul = conditional(h <= 0, zero(mul.ufl_shape), (mul * mul) / h)
        velocity_ur = conditional(h <= 0, zero(mul.ufl_shape), (mur * mur) / h)
        velocityr = conditional(h <= 0, zero(mul.ufl_shape), (mur * mvr) / h)
        velocity_vr = conditional(h <= 0, zero(mvr.ufl_shape), (mvr * mvr) / h)
        velocity_vl = conditional(h <= 0, zero(mvl.ufl_shape), (mvl * mvl) / h)

        F1r = as_vector((mur,
                         velocity_ur + ((gravity / 2) * (h * h)),
                         velocityr))
        F2r = as_vector((mvr,
                         velocityr,
                         velocity_vr + ((gravity / 2) * (h * h))))

        F1l = as_vector((mul,
                         velocity_ul + ((gravity / 2) * (h * h)),
                         velocityl))
        F2l = as_vector((mvl,
                         velocityl,
                         velocity_vl + ((gravity / 2) * (h * h))))

        F_plus = as_vector((F1r, F2r))
        F_minus = as_vector((F1l, F2l))

        W_plus = as_vector((h, mur, mvr))
        W_minus = as_vector((h, mul, mvl))

        # conditional to prevent dividing by zero
        y = ((c_minus - vr) / (c_minus - c_s)) * (W_plus -
                                                  as_vector((0,
                                                            h * (c_s - vl) * N[0],
                                                            h * (c_s - vl) * N[1])))
        w_plus = conditional(eq(c_minus, c_s), zero(y.ufl_shape), y)

        # conditional to prevent dividing by zero
        y = ((c_plus - vl) / (c_plus - c_s)) * (W_minus -
                                                as_vector((0,
                                                          h * (c_s - vr) * N[0],
                                                          h * (c_s - vr) * N[1])))
        w_minus = conditional(eq(c_plus, c_s), zero(y.ufl_shape), y)

        Flux = ((0.5 * dot(N, F_plus + F_minus)) +
                (0.5 * (-((abs(c_minus) - abs(c_s)) * w_minus) +
                        ((abs(c_plus) - abs(c_s)) * w_plus) +
                        (abs(c_minus) * W_plus) -
                        (abs(c_plus) * W_minus))))

        return Flux
