from sympy import *
import matplotlib.pyplot as plt
from runge_kutta_method_4th_order import runge_kutta_method_4th_order
from euler_method_2nd_order import euler_method_2nd_order


##### - override this block for your task
# #11

# y'' + y'*tanx - y*cos^2 xs = 0

# function for y'' (y'' = y*cos^2 xs - y'*tanx)
def dy_dx_func(x, y, dy_dx):
    return N(y * cos(x) ** 2 - dy_dx * tan(x))


def true_y(x):
    return N(exp(sin(x)) + exp(-1 * sin(x)))


# true y'
def true_dy_dx(x):
    return 0


h = 0.1
x_start = 0
x_end = 1

# initial conditions y(x0), y'(x0)
y0 = true_y(x_start)
dy_dx_0 = true_dy_dx(x_start)

#####


xs = []
x = x_start
while x <= x_end:
    xs.append(x)
    x += h
x_end = xs[-1]

fig, ax = plt.subplots()
ax.set_axisbelow(True)
ax.grid(zorder=-1)

ys = euler_method_2nd_order(dy_dx_func, y0, dy_dx_0, x_start, x_end, h)
print("Euler calculated value: ", ys[-1])
ax.plot(xs, ys)

ys = runge_kutta_method_4th_order(dy_dx_func, y0, dy_dx_0, x_start, x_end, h)
print("R-K 4th order calculated value: ", ys[-1])
ax.plot(xs, ys)

print()
print("true value: ", true_y(x_end).evalf())

true_ys = []
for x in xs:
    true_ys.append(true_y(x))
ax.plot(xs, true_ys)
ax.legend(["Euler's method", "Runge-Kutta 4th ord. method", "True function"])
plt.show()
