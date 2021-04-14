from sympy import *


# returns an array of y(xs) points on the [x_start; x_end] interval
def euler_method_2nd_order(dy_dx_func, y0, dy_dx_0, x_start, x_end, h):
    y = y0
    dy_dx = dy_dx_0

    ys = []
    x = x_start
    while x <= x_end:
        ys.append(y)
        y += N(dy_dx * h)
        dy_dx += N(dy_dx_func(x, y, dy_dx) * h)

        x += h

    return ys
