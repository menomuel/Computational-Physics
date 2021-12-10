import numpy as np
import matplotlib.pyplot as plt


def d_J0(x: np.ndarray, t: np.ndarray) -> np.ndarray:
    return 1/np.pi * np.cos(x * np.sin(t))


def d_J1(x, t):
    return 1/np.pi * np.cos(t - x * np.sin(t))


def gomer_simpson(f, x, a, b, N):  # (1/N**4)
    h = (b - a) / N
    t_l = a
    t_r = a + h
    sum = 0
    for i in range(N):
        sum += (f(x, t_l) + 4 * f(x, (t_l + t_r) / 2) + f(x, t_r)) / 6
        t_l += h
        t_r += h
    return sum * h


N_x = 10 ** 2
delta = 10e-6
N_t = 2 ** 10 #empirical
a = 0
b = np.pi

err = []
x = np.linspace(0, 2 * np.pi, N_x)
for x_i in x[1:-2]:
    arg1 = gomer_simpson(d_J1, x_i, a, b, N_t)
    arg2 = (gomer_simpson(d_J0, x_i+delta, a, b, N_t) - gomer_simpson(d_J0, x_i-delta, a, b, N_t)) / (2*delta)
    err.append(arg1 + arg2)

count = sum(True for e in err if e > 10e-10)
print(count, count/len(err))


plt.plot(x[1:-2], err)
plt.show()
