import numpy as np
import matplotlib.pyplot as plt


def l_i(x, roots, i):
    res = 1
    for j in range(len(roots)):
        if i != j:
            res *= (x - roots[j])
    return res


def P_n_legendre(x, n, y, roots):
    sum = 0
    for i in range(n + 1):
        sum += y(roots[i]) * l_i(x, roots, i) / l_i(roots[i], roots, i)
    return sum


def divided_differences(f, x):
    sum = 0
    for i in range(len(x)):
        denominator = 1
        for j in range(len(x)):
            if i != j:
                denominator *= (x[i] - x[j])
        sum += f(x[i]) / denominator
    return sum


def P_n_newton(x, n, y, roots):
    sum = y(roots[0])
    for k in range(1, n + 1):
        mult = divided_differences(y, roots[:k+1])
        for j in range(k):
            mult *= (x - roots[j])
        sum += mult
    return sum


N = np.linspace(1, 15, 15-4+1, dtype=int)
n = N[-5]
print(n)

k = list(range(n + 1))
roots = [1 + k_i/n for k_i in k]
x = np.linspace(1, 4, n * 100)

fig, ax = plt.subplots(1, 2)

P_n = [P_n_legendre(x, n, np.log, roots), P_n_newton(x, n, np.log, roots)]
names = ["Legendre P_n", "Newton P_n"]

for i in range(2):
    #ax[i].plot(x, P_n[i], color='b', label='P_n')
    #ax[i].plot(x, np.log(x), color='r', label='Func')
    #ax[i].scatter(roots, np.log(roots), color='g', label='Roots')
    ax[i].plot(x, P_n[i] - np.log(x), label='Diff')
    ax[i].legend()
    ax[i].title.set_text(names[i])


plt.show()


