import matplotlib.pyplot as plt
import numpy as np




x_dichotomy_history = []
y_dichotomy_history = []
iter_dichotomy = 0
def dichotomy(f, a, b, err):
    if f(a) * f(b) > 0:
        raise ValueError("Bag args")
    mid = (a + b) / 2
    x_dichotomy_history.append(mid)
    val = f(mid)
    y_dichotomy_history.append(val)
    #currErr = (b - a) / 2

    if np.abs(val) < err:
        return mid
    if val * f(a) <= 0:
        return dichotomy(f, a, mid, err)
    else:
        return dichotomy(f, mid, b, err)


def iteration_calc(f, x_0, err, recurrent_func):
    x_history = []
    y_history = []

    num_it = 0
    x = x_0
    x_history.append(x)
    val = f(x)
    y_history.append(val)
    while np.abs(val) >= err:
        if num_it > 10000:
            raise RuntimeError('No convergence')

        x = recurrent_func(x)
        x_history.append(x)
        val = f(x)
        y_history.append(val)
        num_it += 1
    return x_history, y_history


def simple_iteration(f, x_0, param, err):
    def simple_iteration_recurrent_func(x_n):
        return x_n - param * f(x_n)

    return iteration_calc(f, x_0, err, simple_iteration_recurrent_func)


def newton(f, df, x_0, err):
    def newton_recurrent_func(x_n):
        return x_n - f(x_n) / df(x_n)

    return iteration_calc(f, x_0, err, newton_recurrent_func)

width = 2
U_0 = 20


def func(ksi):
    return 1 / np.tan(np.sqrt(2 * width * width * U_0 * (1 - ksi))) - np.sqrt(1 / ksi - 1)


def d_func(ksi, delta=1e-9):
    #arg = np.sqrt(2 * width * width * U_0 * (1 - ksi))
    #return -width ** 2 / (np.sin(arg) ** 2 * arg) - 1 / (2 * ksi ** 2 * np.sqrt(1 / ksi - 1) * U_0)
    return (func(ksi + delta) - func(ksi)) / delta


eps = 1e-10
right_bound = 1 - eps
left_bound = 1 - (np.pi / width)**2 / (2*U_0) + eps
if left_bound <= 0:
    left_bound = eps

# Сalculation
dichotomy_sol = dichotomy(func, left_bound, right_bound, eps)
x_simple_sol, y_simple_sol = simple_iteration(func, (left_bound+right_bound)/2, 0.01, eps)
x_newton_sol, y_newton_sol = newton(func, d_func, (left_bound+right_bound)/2, eps)


print("Dichotomy:", dichotomy_sol, "(%d iter)" % len(x_dichotomy_history))
print("Simple iteration:", x_simple_sol[-1], "(%d iter)" % len(x_simple_sol))
print("Newton:", x_newton_sol[-1], "(%d iter)" % len(x_newton_sol))


x_sol = [x_dichotomy_history, x_simple_sol, x_newton_sol]
y_sol = [y_dichotomy_history, y_simple_sol, y_newton_sol]
names = ["Dichotomy", "Simple iteration", "Newton"]


fig = plt.figure()
x = np.linspace(eps, right_bound, 1000)
y = func(x)

for i in range(3):
    plt.subplot(1, 3, i+1)
    plt.title(names[i])
    plt.grid()
    plt.xlim(0, 1)
    plt.ylim(-5, 5)
    plt.plot(x_sol[i], y_sol[i], 'o')
    plt.plot(x_sol[i][-1], y_sol[i][-1], 'x')
    plt.plot(x, y)

plt.show()
