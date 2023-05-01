import numpy as np
import matplotlib.pyplot as plt
import math

# constants and initial conditions
gamma = 5/3
L = 10 # m, x=[-L, L]

# left side
v_left = 0 # m/s
rho_left = 13 # kg/m3
p_left = 10 # atm

# right side
v_right = 0 # m/s
rho_right = 1.3 # kg/m3
p_right = 1 # atm

# v = (rho, u, e)  w = (rho, rho*u, rho*e)
init_conds = [rho_left, v_left, p_left/((gamma-1)*rho_left), rho_right, v_right, p_right/((gamma-1)*rho_right)]

time_end = 2
cfl_max = 0.01
tau = 1e-3
h = 0.2

x_limits = [-L, L]
t_limits = [0, time_end]


def sound_speed(e):
    return math.sqrt(gamma * (gamma - 1) * e)


def A(u, e):
    return np.array([[0, 1, 0], [-u * u, 2 * u, gamma - 1], [-gamma * u * e, gamma * e, u]])


def OmegaT(u, e):
    c = sound_speed(e)
    return np.array([[-u * c, c, gamma - 1], [-c * c, 0, gamma - 1], [u * c, -c, gamma - 1]])


def LambdaAbs(u, e):
    c = sound_speed(e)
    return np.array([[math.fabs(u + c), 0, 0], [0, math.fabs(u), 0], [0, 0, math.fabs(u - c)]])


def OmegaTinv(u, e):
    c = sound_speed(e)
    return np.array([[1, -2, 1], [u + c, -2 * u, u - c], [c * c / (gamma - 1), 0, c * c / (gamma - 1)]]) / (2 * c * c)


def solve(init_conds_, x_limits_, t_limits_, gamma_, tau_, h_, cfl_max_):
    x_n_steps = math.floor((x_limits_[1] - x_limits_[0]) / h_)
    result = np.empty((x_n_steps + 1, 1, 4))
    x_coord = lambda i: x_limits_[0] + i * h_
    to_v = lambda w: np.array([w[0], w[1] / w[0], w[2] / w[0], (gamma - 1) * w[2]])
    to_w = lambda v: np.array([v[0], v[0] * v[1], v[0] * v[2]])

    # applying initial conditions
    for i in range(x_n_steps + 1):
        if x_coord(i) < 0:
            result[i, 0] = np.array(init_conds_[:3] + [(gamma - 1) * init_conds_[0] * init_conds_[2]])
        else:
            result[i, 0] = np.array(init_conds_[3:] + [(gamma - 1) * init_conds_[3] * init_conds_[5]])

    curr_time = t_limits_[0]
    while curr_time < t_limits_[1]:
        t_layer = np.empty((x_n_steps + 1, 1, 4))
        for i in range(1, x_n_steps):
            matA = A(result[i, -1, 1], result[i, -1, 2])
            matOmegaT = OmegaT(result[i, -1, 1], result[i, -1, 2])
            matLambdaAbs = LambdaAbs(result[i, -1, 1], result[i, -1, 2])
            matOmegaTinv = OmegaTinv(result[i, -1, 1], result[i, -1, 2])
            lambda_max = matLambdaAbs.max()

            # calc time step
            curr_tau = tau_
            curr_cfl = curr_tau * lambda_max / h_
            while curr_cfl > cfl_max_:
                curr_tau /= 2
                curr_cfl = curr_tau * lambda_max / h_

            w = to_w(result[i, -1, :3]) - tau_ / (2 * h_) * matA @ (
                        to_w(result[i + 1, -1, :3]) - to_w(result[i - 1, -1, :3])) + tau_ / (2 * h) * matOmegaTinv @ (
                            matLambdaAbs @ (matOmegaT @ (
                                to_w(result[i + 1, -1, :3]) - 2 * to_w(result[i, -1, :3]) + to_w(
                            result[i - 1, -1, :3]))))
            t_layer[i, 0] = to_v(w + [0])

        # applying border conditions
        t_layer[0, 0] = t_layer[1, 0]
        t_layer[-1, 0] = t_layer[-2, 0]
        result = np.concatenate((result, t_layer), axis=1)
        curr_time += curr_tau
    return result


sol = solve(init_conds, x_limits, t_limits, gamma, tau, h, cfl_max) # pho, u, e, p
print(sol.shape)

nx = math.floor((x_limits[1] - x_limits[0])/h)
x_arr = np.array([x_limits[0] + i * h for i in range(nx + 1)])

plt.scatter(x_arr, sol[:, 2000, 3], label=r'Численное решение')
plt.legend()
plt.title('Задача Римана о распаде произвольного разрыва')
plt.grid()
plt.xlabel('x, m')
plt.ylabel(r'P, atm')
plt.savefig('data/P.png')
plt.show()

