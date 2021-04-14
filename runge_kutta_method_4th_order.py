from sympy import *


# returns an array of y(xs) points on the [x_start; x_end] interval
def runge_kutta_method_4th_order(dy_dx_func, y0, dy_dx_0, x_start, x_end, h):
    y = y0
    dy_dx = dy_dx_0
    x = x_start
    ys = []
    while x <= x_end:
        ys.append(y)

        # Ks are for y, Ls are for dy_dx
        K1 = h * dy_dx
        L1 = h * dy_dx_func(x, y, dy_dx)

        K2 = h * (dy_dx + L1 / 2)
        L2 = h * dy_dx_func(x + h / 2, y + K1 / 2, dy_dx + L1 / 2)

        K3 = h * (dy_dx + L2 / 2)
        L3 = h * dy_dx_func(x + h / 2, y + K2 / 2, dy_dx + L2 / 2)

        K4 = h * (dy_dx + L3 / 2)
        L4 = h * dy_dx_func(x + h / 2, y + K3 / 2, dy_dx + L3 / 2)

        Dy = (K1 + 2 * K2 + 2 * K3 + K4) / 6
        y += Dy

        Ddy_dx = (L1 + 2 * L2 + 2 * L3 + L4) / 6
        dy_dx += Ddy_dx

        x += h

    return ys
