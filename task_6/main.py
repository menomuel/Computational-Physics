import numpy as np
import matplotlib.pyplot as plt


# dx/dt = -x, x(0) = 1, 0 < t < 3
def right_side(t, x):
    return -x


def ground_truth(t):
    return np.exp(-t)


def euler(t, f, init):
    x = np.zeros(t.shape[0])
    x[0] = init
    for i in range(x.shape[0] - 1):
        h = t[i + 1] - t[i]
        x[i + 1] = x[i] + h * f(t[i], x[i])
    return x


def runge_kutta_2(t, f, init, alpha=0.75):
    x = np.zeros(t.shape[0])
    x[0] = init
    for i in range(x.shape[0] - 1):
        h = t[i + 1] - t[i]
        x[i + 1] = x[i] + h * ((1 - alpha) * f(t[i], x[i]) + alpha * f(t + h/(2*alpha), x[i] + h/(2*alpha) * f(t[i], x[i])))
    return x


def runge_kutta_4(t, f, init):
    x = np.zeros(t.shape[0])
    x[0] = init
    for i in range(x.shape[0] - 1):
        h = t[i + 1] - t[i]
        k1 = f(t[i], x[i])
        k2 = f(t[i] + h/2, x[i] + h/2 * k1)
        k3 = f(t[i] + h / 2, x[i] + h / 2 * k2)
        k4 = f(t[i] + h, x[i] + h * k3)
        x[i + 1] = x[i] + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    return x


eps = 10e-3
t = np.linspace(0, 3, 10 ** 1)
x_ground = ground_truth(t)
x_euler = euler(t, right_side, init=1)  # x(0)=1
x_rk_2 = runge_kutta_2(t, right_side, init=1)
x_rk_4 = runge_kutta_4(t, right_side, init=1)

plt.plot(t, x_ground, label='ground')
plt.plot(t, x_euler, label='euler')
plt.plot(t, x_rk_2, label='runge_kutta (II)')
plt.plot(t, x_rk_4, label='runge_kutta (IV)')
plt.legend()
plt.show()
