import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

'''
MAIN: y'' = sin(x)
Solution: y = c_1 * x + c_2 - sin(x)

y'' + p(x) * y' + q(x) * y = r(x)

c_t1 * y(a) + c_t2 * y'(a) = c_t
d_t1 * y(b) + d_t2 * y'(b) = d_t

Assume y_(j-1) = ksi_j * y_j + nu_j

'''

def groud_truth(x, init):
    return -np.sin(x) + (init[1] - init[0]) / np.pi * x + init[0]

def system(x, init):
    n = x.shape[0]
    h = x[1] - x[0]

    a = np.ones(n) / h ** 2# lower diagonal
    b = -2 * np.ones(n) / h ** 2 # main diagonal
    c = np.ones(n) / h ** 2# upper diagonal 
    d = np.sin(x)

    a[0] = 0
    c[n - 1] = 0
    # Boundaries
    c[0] = 0
    d[0] = init[0] * b[0]
    a[n - 1] = 0
    d[-1] = init[1] * b[n - 1]

    return a, b, c, d

#Tridiagonal matrix algorithm
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


def updateGraph():
    global slider_left
    global slider_right
    global graph_axes

    init[0] = slider_left.val 
    init[1] = slider_right.val

    N = 10 ** 4
    x = np.linspace(0, np.pi, N)
    a, b, c, d = system(x, init)
    
    y = tma(a, b, c, d)
    gt = groud_truth(x, init)

    graph_axes.clear()
    #graph_axes.plot(x, y, label='tma')
    #graph_axes.plot(x, gt, label='gt')
    graph_axes.plot(x, y - gt, label='tma')
    graph_axes.legend()

    plt.draw()
'''
#Tridiagonal matrix algorithm
def tma(x, init):
    h = x[1] - x[0]
    ksi = np.zeros_like(x)
    nu = np.zeros_like(x)
    y = np.zeros_like(x)

    c_t1, c_t2, c_t = init[0][0], init[0][1], init[0][2]
    d_t1, d_t2, d_t = init[1][0], init[1][1], init[1][2]

    ksi[0] = - c_t2 / (c_t1 * h - c_t2)
    nu[0] = c_t * h / (c_t1 * h - c_t2)
    for k in range(x.shape[0] - 1):
        ksi[k + 1] = 1 / (2 - ksi[k])
        nu[k + 1] = (nu[k] - h ** 2 * np.sin(x[k])) / (2 - ksi[k])

    y[-1] = (d_t2 * nu[-1] - d_t * h) / (d_t2 * (1 - ksi[-1]) + d_t1 * h)

    for j in reversed(range(1, len(x))):
        y[j - 1] = y[j] * ksi[j] + nu[j]

    return y


def updateGraph():
    global slider_left
    global slider_right
    global graph_axes

    init[0][2] = slider_left.val 
    init[1][2] = slider_right.val

    N = 10 ** 4
    x = np.linspace(0, np.pi, N)
    y = tma(x, init)

    graph_axes.clear()
    graph_axes.plot(x, y, label='tma')
    #graph_axes.plot(x, -np.sin(x), label='tma')
    #graph_axes.plot(x, y + np.sin(x), label='tma')
    graph_axes.legend()

    plt.draw()

# [c_t1, c_t2, c_t], [d_t1, d_t2, d_t]
#init = np.array([[1, 0, 0], [1, 0, 0]], dtype=float) # y(0) = 0, y(pi) = 0
'''

def onChangeSlider(value):
    updateGraph()


init = [0, 0] # # y(0) = 0, y(pi) = 0

fig, graph_axes = plt.subplots()
fig.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.2)

ax_slider_left = plt.axes([0.1, 0.1, 0.8, 0.04])
slider_left = Slider(ax_slider_left, label='y(0)', valmin=-5, valmax=5, valinit=0, valstep=0.5) 
ax_slider_right = plt.axes([0.1, 0.05, 0.8, 0.04])
slider_right = Slider(ax_slider_right, label='y(pi)', valmin=-5, valmax=5, valinit=0, valstep=0.5) 

slider_left.on_changed(onChangeSlider)
slider_right.on_changed(onChangeSlider)

updateGraph()

plt.show()
