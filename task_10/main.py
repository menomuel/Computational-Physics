import matplotlib.pyplot as plt
import numpy as np

def phi(t):
    return 0


def psi(x, L=1):
    return x * (1 - x / L) ** 2


def tma(a, b, c, d):
    n = len(a)
    y = np.zeros(n)

    for i in range(1, n):
        w = a[i] / b[i - 1]
        b[i] = b[i] - w * c[i - 1]
        d[i] = d[i] - w * d[i - 1]

    y[n - 1] = d[n - 1] / b[n - 1]
    for i in reversed(range(n - 1)):
        y[i] = (d[i] - c[i] * y[i + 1]) / b[i]
    
    return y


def crank_nicolson(t: np.array, x: np.array, _phi_1=phi, _phi_2=phi, _psi=psi):
    dt = t[1] - t[0]
    h = x[1] - x[0]

    a_j = - dt / h ** 2
    c_j = a_j
    b_j = 1 + 2 * dt / h ** 2

    u = np.zeros(shape=(t.shape[0], x.shape[0]))
    u[0, :] = _psi(x)

    max_T = [np.max(u[0,:])]

    alpha = np.zeros(x.shape[0] - 1)
    beta = np.zeros(x.shape[0] - 1)

    for n in range(t.shape[0] - 1):
        alpha[0] = 0
        beta[0] = _phi_1(t[n + 1])

        for j in range(1, x.shape[0] - 1):
            ksi_n_j = u[n, j] + (dt / 2 / h ** 2) * (u[n, j + 1] - 2 * u[n, j] + u[n, j - 1])
            alpha[j] = -a_j / (b_j + c_j * alpha[j - 1])
            beta[j] = (ksi_n_j - c_j * beta[j - 1]) / (b_j + c_j * alpha[j - 1])

        u[n + 1, x.shape[0] - 1] = _phi_2(t[n + 1])

        for j in reversed(range(x.shape[0] - 1)):
            u[n + 1, j] = alpha[j] * u[n + 1, j + 1] + beta[j]
        
        max_T.append(np.max(u[n + 1, :]))

    return u, max_T


L = 1
N = 100

x = np.linspace(0, L, N)
t = np.linspace(0, 4, N)

u, max_T = crank_nicolson(t, x)

reduce_coef = 20
for i in range(len(t) // reduce_coef):
    plt.plot(x, u[i * reduce_coef, :], label=f'crank-nicolson {i * reduce_coef}')

#plt.plot(t, max_T, label='crank-nicolson')
#plt.plot(t, np.exp(-t), label='exp(-t)')
plt.legend()
plt.show()

