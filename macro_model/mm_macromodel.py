import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


def A(t):
    """
    transition rate matrix
    :param t: main timeline (from integrator)
    :return: normalized transition rate matrix with stimulus
    """

    test_matrix = [[-4 * imp(t, 2, 5, 10),    2*2 * imp(t, 2, 5, 10),    0,     0,    0],
                   [0.5,  -0.5 - 2 * imp(t, 2, 5, 10),    2 * imp(t, 2, 5, 10),     0,    0],
                   [0,  2*0.5,  -2.6,   0.1,  1.5],
                   [0,      0,   0.3,  -0.3,    0],
                   [0,      0,     1,     0,   -1]]

    return np.array(test_matrix)


def imp(t, a, b, v):
    """
    time dependant transition rate step modification
    :param t: main timeline (from integrator)
    :param a: step start time
    :param b: step end time
    :param v: step height (rate multiplication)
    :return: multiplier scalar for given time
    """
    if a < t < b:
        return v
    else:
        return 0


def dpdt(t, P, A):
    """
    Differential equation: dP/dt = A * P
    :param t: blank parameter for integrator and stimulus
    :param P: zero-time + blank to fill matrix
    :param A: transition rate matrix row normalized
    :return: differential equation for integrator
    """
    return np.dot(P, A(t))


def solveKFW(t, dpdt, A, P0, t0):
    """
    Discretize and integrate the Kolmogorov Forward Equation.
    :param t: main timeline (from integrator)
    :param dpdt: differential equation for integrator
    :param A: normalized transition rate matrix with stimulus
    :param P0: starting conditions
    :param t0: starting time
    :return:
    """

    rk45 = ode(dpdt).set_integrator('dopri5')
    rk45.set_initial_value(P0, t0).set_f_params(A)
    samples = 1000
    dt = t/samples

    P = np.zeros((samples + 1, 5))
    T = np.zeros(samples + 1)
    P[0, :] = P0
    T[0] = t0

    idx = 0
    while rk45.successful() and rk45.t < t:
        rk45.integrate(rk45.t+dt)
        P[idx, :] = np.array(rk45.y)
        T[idx] = rk45.t
        idx += 1

    return P, T

P, T = solveKFW(20, dpdt, A, np.array([1, 0, 0, 0, 0]), 0)

plt.plot(T, P[:, 0], 'r--', linewidth=2.0)
plt.plot(T, P[:, 1], 'b--', linewidth=2.0)
plt.plot(T, P[:, 2], 'g--', linewidth=2.0)
plt.plot(T, P[:, 3], 'y--', linewidth=2.0)
plt.plot(T, P[:, 4], 'c-', linewidth=4.0)
plt.legend(['R', 'AR', 'A2R', 'A2D', 'A2O'])
plt.xlabel('time []')
plt.ylabel('state probability')
plt.show()