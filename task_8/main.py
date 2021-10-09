import numpy as np
import matplotlib.pyplot as plt


def ground_truth(t, u_0, v_0):
    alpha = u_0 + v_0
    beta = -2 * v_0 - u_0
    u = 2 * alpha * np.exp(-t) + beta * np.exp(-1000 * t)
    v = -alpha * np.exp(-t) - beta * np.exp(-1000 * t)
    return u, v


def f(u, v):
    return 998 * u + 1998 * v


def g(u, v):
    return -999 * u - 1999 * v


def F(u, v):
    return np.array([f(u, v), g(u, v)])


def euler(t, F, U_0):
    U = np.zeros((U_0.shape[0], t.shape[0]))
    U[:, 0] = U_0
    for i in range(t.shape[0] - 1):
        h = t[i + 1] - t[i]  # < 0.002
        U[:, i + 1] = U[:, i] + h * F(U[:, i][0], U[:, i][1])
    return U


def implicit_euler(t, F, U_0):
    J = np.array([[998, 1998], [-999, -1999]])
    E = np.eye(2)

    U = np.zeros((U_0.shape[0], t.shape[0]))
    U[:, 0] = U_0
    h = t[1] - t[0]  # < 0.002
    mat = np.linalg.inv(E / h - J)
    for i in range(t.shape[0] - 1):
        U[:, i+1] = U[:, i] + np.dot(mat, F(U[:, i][0], U[:, i][1]))
    return U


U_0 = np.array([1, 1])
t = np.linspace(0, 1, 10 ** 3)

gt_sol = ground_truth(t, U_0[0], U_0[1])
eu_sol = euler(t, F, U_0)
im_eu_sol = implicit_euler(t, F, U_0)

names = ['u(t)', 'v(t)']
fig, ax = plt.subplots(1, 2)

for i in range(2):
    ax[i].plot(t, gt_sol[i], label='ground truth')
    ax[i].plot(t, eu_sol[i], label='explicit euler')
    ax[i].plot(t, im_eu_sol[i], label='implicit euler')
    ax[i].title.set_text(names[i])
    ax[i].legend()

plt.show()