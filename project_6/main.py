import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def phi(ax, L=1):
    return 1 - ax ** 2 / L ** 2


#Tridiagonal matrix algorithm
def tma(a, b, c, d):
    n = len(a)
    sol = np.zeros(n)

    for i in range(1, n):
        w = a[i] / b[i - 1]
        b[i] = b[i] - w * c[i - 1]
        d[i] = d[i] - w * d[i - 1]

    sol[n - 1] = d[n - 1] / b[n - 1]
    for i in reversed(range(n - 1)):
        sol[i] = (d[i] - c[i] * sol[i + 1]) / b[i]
        
    return sol

def local_unidimentional_method(t: np.array, x: np.array, y: np.array):
    dt = t[1] - t[0]
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    u = np.zeros(shape=(t.shape[0], x.shape[0], y.shape[0]))
    w = np.zeros(shape=(2 * t.shape[0], x.shape[0], y.shape[0]))

    for j in range(0, 2 * t.shape[0] - 2, 2):
        w[:, 0, :] = phi(y)
        w[:, -1, :] = phi(y)
        w[:, :, 0] = phi(x)
        w[:, :, -1] = phi(x)
        
        # Fill matrix (x - implicit, y - explicit)
        A_x = np.full(x.shape[0], - dt /2/ dx**2)
        B_x = np.full(x.shape[0], 1 + dt /2/ dx**2)
        C_x = np.full(x.shape[0], - dt /2/ dx**2)
        D_x = np.zeros(x.shape[0])
        # Bound conditions 
        B_x[0] = 1
        B_x[-1] = 1

        for m in range(1, y.shape[0] - 1):
            D_x[0] = phi(y[m])
            D_x[-1] = phi(y[m])
            for n in range(1, x.shape[0]-1):
                D_x[n] = w[j,n,m] + dt/2/dy**2 * (w[j,n,m+1] - 2*w[j,n,m] + w[j,n,m-1])
            w[j+1,m,:] = tma(A_x, B_x, C_x, D_x)
        
        # Fill matrix (x - explicit, y - implicit)
        A_y = np.full(y.shape[0], - dt /2/ dy**2)
        B_y = np.full(y.shape[0], 1 + dt / dy**2)
        C_y = np.full(y.shape[0], - dt /2/ dy**2)
        D_y = np.zeros(y.shape[0])
        # Bound conditions
        B_y[0] = 1
        B_y[-1] = 1

        for n in range(1, x.shape[0] - 1):
            D_y[0] = phi(x[n])
            D_y[-1] = phi(x[n])
            for m in range(1, y.shape[0]-1):
                D_y[m] = w[j+1,n,m] + dt/2/dx**2 * (w[j+1,n+1,m] - 2*w[j+1,n,m] + w[j+1,n-1,m])
            w[j+2,:,n] = tma(A_y, B_y, C_y, D_y)

    for i in range(0, t.shape[0] - 1):
        u[i] =  w[2*i]
    return u
        

def pseudoviscosity_method(t: np.array, x: np.array, y: np.array, eps=1000):
    dt = t[1] - t[0]
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    if(dt > 0.5 / (1/dx**2 + 1/dy**2)):
        print("Scheme not stable. Return to base. over.")
        return

    u = np.zeros(shape=(t.shape[0], x.shape[0], y.shape[0]))
    # Bound conditions
    u[:, 0, :] = phi(y)
    u[:, -1, :] = phi(y)
    u[:, :, 0] = phi(x)
    u[:, :, -1] = phi(x)

    for n in range(t.shape[0] - 1):
        diff = u[n,:,:]
        #while np.max(np.abs(diff)) > eps:
        for k in range(eps):
            for i, j in zip(range(1, x.shape[0] - 1), range(1, y.shape[0] - 1)):
                u[n+1,i,j] = u[n,i,j] + dt * ((u[n,i+1,j] - 2*u[n,i,j] + u[n,i-1,j])/dx**2 + \
                                            +(u[n,i,j+1] - 2*u[n,i,j] + u[n,i,j-1])/dy**2)
                #print(f'n = {n}; i = {i}; j = {j}')
                #print(u[n,i+1,j], u[n,i,j], u[n,i-1,j])
                #print(u[n,i,j+1],u[n,i,j], u[n,i,j-1])
            diff = u[n+1,:,:] - u[n,:,:]

    #print(f'Converge in {n} iterations')
    return u

def update():
    global slider_time
    global graph_axes
    global u
    global clb
    
    graph_axes.clear()
    clb.remove()
    im = graph_axes.imshow(u[slider_time.val, :, :], cmap='hot')
    clb=plt.colorbar(im, ax=graph_axes)

    plt.draw()

def onChangeSlider(value):
    update()

L = 1
N = 256
t = np.linspace(0, 1, N)
x = np.linspace(-L, L, 16)
y = np.linspace(-L, L, 16)

u = pseudoviscosity_method(t, x, y)
#u = local_unidimentional_method(t, x, y)

fig, graph_axes = plt.subplots()
fig.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.2)

ax_slider_time = plt.axes([0.05, 0.05, 0.8, 0.04])
slider_time = Slider(ax_slider_time, label='t', valmin=0, valmax=N-2, valinit=0, valstep=1) 

slider_time.on_changed(onChangeSlider)

im = graph_axes.imshow(u[0, :, :], cmap='hot')
clb=plt.colorbar(im, ax=graph_axes)

plt.show()