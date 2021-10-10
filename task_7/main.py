import numpy as np
import matplotlib.pyplot as plt


def right_side(x, y, coef):
    f1 = coef[0][0] * x - coef[0][1] * x * y
    f2 = coef[1][0] * x * y - coef[1][1] * y
    return f1, f2


def runge_kutta_2(t, f, init, coef, alpha=0.75):
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    x[0], y[0] = init[0], init[1]
    for i in range(t.shape[0] - 1):
        h = t[i + 1] - t[i]
        f_i = f(x[i], y[i], coef)
        f_h = f(x[i] + h/(2*alpha) * f_i[0], y[i] + h/(2*alpha) * f_i[1], coef)
        x[i + 1] = x[i] + h * ((1 - alpha) * f_i[0] + alpha * f_h[0])
        y[i + 1] = y[i] + h * ((1 - alpha) * f_i[1] + alpha * f_h[1])
    return x, y


coef = [[10, 2], [2, 10]]
init = np.array([[10, 20], [2, 10], [5, 5], [1, 1], [20, 4]])  # x(0) = 30, y(0) = 80
t = np.linspace(0, 1, 10 ** 3)

fig, ax = plt.subplots(2, 1)
ax[0].title.set_text("Phase space")
ax[1].title.set_text("Time space")

choice = 0
for i in range(init.shape[0]):
    x, y = runge_kutta_2(t, right_side, init[i], coef)
    ax[0].plot(x, y, label=init[i])
    if i == choice:
        ax[1].plot(t, x, label=f'{init[i]}_x')
        ax[1].plot(t, y, label=f'{init[i]}_y')

ax[0].legend()
ax[1].legend()
plt.show()
