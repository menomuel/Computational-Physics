import matplotlib.pyplot as plt
import numpy as np


def phi(t):
    return 0


def psi(x, L=1):
    return x * (1 - x - L) ** 2


def crank_nicolson(t: np.array, x: np.array, _phi_1=phi, _phi_2=phi, _psi=psi, sigma=1):
    dt = t[1] - t[0]
    h = x[1] - x[0]

    a_j = - sigma * dt / (2 * h ** 2)
    c_j = a_j
    b_j = 1 + sigma * dt / h ** 2


    u = np.zeros(shape=(t.shape[0], x.shape[0]))
    u[0, :] = _psi(x)

    max_T = [np.amax(u[0, :])]

    alpha = np.zeros_like(x)
    beta = np.zeros_like(x)

    for n in range(t.shape[0] - 1):
        alpha[0] = 0
        beta[0] = _phi_1(t[n + 1])

        for j in range(1, x.shape[0] - 1):
            ksi_n_j = u[n, j] + sigma * dt /(2 * h ** 2) * (u[n, j + 1] - 2 * u[n, j] + u[n, j - 1])
            alpha[j] = - a_j / (b_j + c_j * alpha[j - 1])
            beta[j] = (ksi_n_j - c_j * beta[j - 1]) / (b_j + c_j * alpha[j - 1])
        
        u[n + 1, x.shape[0] - 1] = _phi_2(t[n + 1])

        for j in reversed(range(x.shape[0] - 1)):
            u[n + 1, j] = alpha[j] * u[n + 1, j + 1] + beta[j]

        max_T.append(np.amax(u[n + 1, :]))

    return u, max_T


L = 1
N = 10 ** 3

x = np.linspace(0, L, N // 10)
t = np.linspace(0, 10, N)

u, max_T = crank_nicolson(t, x)

plt.plot(t, max_T, label='crank-nicolson')
plt.plot(t, np.exp(-t), label='exp(-t)')
plt.legend()
plt.show()

