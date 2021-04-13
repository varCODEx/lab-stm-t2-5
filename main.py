from sympy import *
import matplotlib.pyplot as plt

##### - override this block for your task
# #11

#y'' + y'*tanx - y*cos^2 x = 0

# function for y'' (y'' = y*cos^2 x - y'*tanx)
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

def euler_method_2nd_order(dy_dx_func, y0, dy_dx_0, x_start, x_end, h):
    y = y0
    dy_dx = dy_dx_0

    ys = [y]
    x_ = x_start
    while x_ < x_end:
        y += N(dy_dx * h)
        dy_dx += N(dy_dx_func(x_, y, dy_dx) * h)

        x_ += h

        ys.append(y)

    return ys


x = []
x_ = x_start
while x_ < x_end:
    x.append(x_)
    x_ += h
x_end = x[-1]


ys = euler_method_2nd_order(dy_dx_func, y0, dy_dx_0, x_start, x_end, h)
print("calculated value: ", ys[-1])

print()
print("true value: ", true_y(x_end).evalf())

fig, ax = plt.subplots()
ax.set_axisbelow(True)
ax.grid(zorder=-1)
ax.plot(x, ys)
true_ys = []
for x_ in x:
    true_ys.append(true_y(x_))
ax.plot(x, true_ys)
ax.legend(["Euler's method", "True function"])
plt.show()
