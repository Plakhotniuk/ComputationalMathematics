import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def threeDiagProg(a_, b_, c_, d_):
    n_ = len(b_)
    p_, q_ = np.zeros(n_), np.zeros(n_)
    c_ = c_ + [0]
    p_[0] = c_[0] / b_[0]
    q_[0] = d_[0] / b_[0]
    for i in range(1, n_):
        p_[i] = c_[i] / (b_[i] - (a_[i - 1] * p_[i - 1]))
        q_[i] = (d_[i] - (a_[i - 1] * q_[i - 1])) / (b_[i] - (a_[i - 1]*p_[i - 1]))

    x_ = np.zeros(n_)
    m_ = n_ - 1

    x_[m_] = q_[-1]
    for i in range(n_ - 2, -1, -1):
        x_[i] = q_[i] - (x_[i + 1] * p_[i])
    return x_


lm = 1e-4
u_an = lambda t, x, y: cos(pi*x)*sin(5*pi*y)*exp(-50*(pi**2)*lm*t)
u0x = lambda t, y: sin(5*pi*y)*exp(-50*(pi**2)*lm*t)
u1x = lambda t, y: -sin(5*pi*y)*exp(-50*(pi**2)*lm*t)


u = np.zeros([])

def calc(N, Nx, Ny):
    global u
    err_max = 0
    tt = np.linspace(0, 1, N+1)
    xx = np.linspace(0, 1, Nx+1)
    yy = np.linspace(0, 1, Ny+1)

    tau = np.diff(tt)[0]
    hx = np.diff(xx)[0]
    hy = np.diff(yy)[0]

    u = np.zeros([N+1, Nx+1, Ny+1])

    for i in range(Nx + 1):
        for j in range(Ny + 1):
            u[0, i, j] = cos(pi*xx[i]) * sin(5*pi*yy[j])

    u[:, 0, :] = [[u0x(t, y) for y in yy] for t in tt]
    u[:, Nx, :] = [[u1x(t, y) for y in yy] for t in tt]

    a = (-50*lm/(hx**2)) - (2*lm/(hy**2)) - (1/tau)
    b = lm/(hy**2)
    c = 25*lm/(hx**2)

    a2 = [0] + [c for m in range(1, Nx)]
    a1 = [0] + [b for m in range(1, Nx)]
    a0 = [1] + [a for l in range(1, Nx)] + [1]
    a_1 = [b for m in range(1, Nx)] + [0]
    a_2 = [c for m in range(1, Nx)] + [0]

    for n in range(N):
        if n % 2 == 0:
            for j in range(1, Ny):
                d = [u0x(tt[n+1], yy[j])] + [-(u[n, i, j]/tau)-(b*u[n+1, i, j-1])-(b*u[n, i, j+1]) for i in range(1, Nx)] + [u1x(tt[n+1], yy[j])]
                u[n+1, :, j] = threeDiagProg(a_2, a0, a2, d)
        else:
            for i in range(1, Nx):
                d = [0] + [-(u[n, i, j]/tau)-(c*u[n+1, i-1, j])-(c*u[n, i+1, j]) for j in range(1, Ny)] + [0]
                u[n+1, i, :] = threeDiagProg(a_1, a0, a1, d)

        real = np.cos(np.pi * xx) * np.sin(5 * np.pi * yy) * np.exp(- np.pi * np.pi * 50 * 1e-4 * (tt[n+1] + 1) * tau)
        err = np.linalg.norm(real - u[n+1])
        if err > err_max:
            err_max = err

    return err

errors = []
for i in [10, 20, 40, 80, 120, 150]:
    errors.append(calc(100, i, i))

# print(errors)
x = np.log([10, 20, 40, 80, 120, 150])
y = np.log(errors)
p = np.polyfit(x, y, 1)
print(p)
plt.plot(x, y, label=f'k = {p[0]}')
plt.xlabel('log(x)')
plt.ylabel('log(err)')
plt.legend()
plt.grid()
plt.show()
# sns.heatmap(u[-1, :, :])
