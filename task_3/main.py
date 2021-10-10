import numpy as np
import matplotlib.pyplot as plt


def f_a(x):
    return 1 / (1 + x**2)


def f_b(x):
    return x ** 1/3 * np.exp(np.sin(x))


def trapezoidal(f, a, b, N):
    x = np.linspace(a, b, N)
    h = x[1] - x[0]
    return h * (np.sum(f(x)) - (f(x[0]) + f(x[-1]))/2)


def simpson(f, a, b, N):  # (1/N**4)
    h = (b - a) / N
    x_l = a
    x_r = a + h
    sum = 0
    for i in range(N):
        sum += (f(x_l) + 4 * f((x_l + x_r) / 2) + f(x_r)) / 6
        x_l += h
        x_r += h

    return sum * h

a = -1
b = 1

ethalon = np.pi / 2
n_a = []
err_tr_a = []
err_simp_a = []
for i in range(1, 20):
    N = pow(2, i)
    n_a.append(N)
    err_tr_a.append(abs(ethalon - trapezoidal(f_a, a, b, N)))
    err_simp_a.append(abs(ethalon - simpson(f_a, a, b, N)))

print(f'k_tr_a={np.log(abs(err_tr_a[-1] / err_tr_a[0])) / np.log(abs(n_a[-1] / n_a[0]))}')
print(f'k_simp_a={np.log(abs(err_simp_a[-1] / err_simp_a[0])) / np.log(abs(n_a[-1] / n_a[0]))}')

n_b = []
tr_b = []
simp_b = []
for i in range(1, 20):
    N = pow(2, i)
    n_b.append(N)
    tr_b.append(trapezoidal(f_b, a, b, N))
    simp_b.append(simpson(f_b, a, b, N))

err_tr_b = []
err_simp_b = []
for i in range(len(tr_b)-1):
    err_tr_b.append(abs(tr_b[i+1]-tr_b[i]))
    err_simp_b.append(abs(simp_b[i + 1] - simp_b[i]))

print(f'k_tr_a={np.log(abs(err_tr_b[-1] / err_tr_b[0])) / np.log(abs(n_b[-1] / n_b[1]))}')
print(f'k_simp_a={np.log(abs(err_simp_b[-1] / err_simp_b[0])) / np.log(abs(n_b[-1] / n_b[1]))}')

fig, ax = plt.subplots(1, 2)

ax[0].set_xscale('log', base=2)
ax[0].set_yscale('log')
ax[0].title.set_text('Integral (a)')
ax[0].plot(n_a, err_tr_a, label='Trapezoidal')
ax[0].plot(n_a, err_simp_a, label='Simpsons')
ax[0].grid()
ax[0].legend()

ax[1].set_xscale('log', base=2)
ax[1].set_yscale('log')
ax[1].title.set_text('Integral (b)')
ax[1].plot(n_b[1:], err_tr_b, label='Trapezoidal')
ax[1].plot(n_b[1:], err_simp_b, label='Simpsons')
ax[1].grid()
ax[1].legend()

plt.show()