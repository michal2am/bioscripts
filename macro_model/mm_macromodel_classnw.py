import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


def A(t):
    '''
    Transition probability matrix
    '''
    test_matrix = [[-4 * imp(t),    2*2 * imp(t),    0,     0,    0],
                   [0.5,  -0.5 - 2 * imp(t),    2 * imp(t),     0,    0],
                   [0,  2*0.5,  -2.6,   0.1,  1.5],
                   [0,      0,   0.3,  -0.3,    0],
                   [0,      0,     1,     0,   -1]]

    return np.array(test_matrix)


def imp(t):
    if 4 < t < 7:
        return 10
    else:
        return 0


class NumericalSolver:

    def __init__(self, P0, t0, t):
        self.P0 = P0
        self.t0 = t0
        self.t = t
        self.P, self.T = self.solveKFW()

    def dpdt(self, P, A):
        """
        Differential equation: dP/dt = A * P
        :param t: blank parameter for integrator
        :param P: zero-time condition matrix
        :param A: transition rate matrix row normalized
        :return: differential equation for integrator
        """
        print(np.dot(P, A(self.t)))
        return np.dot(P, A(self.t))

    def solveKFW(self):
        """
        Discretize and integrate the Kolmogorov Forward Equation.
        :param Tf:
        :param dpdt:
        :param A:
        :param P0:
        :param t0:
        :return:
        """

        rk45 = ode(self.dpdt).set_integrator('dopri5')
        rk45.set_initial_value(self.P0, self.t0).set_f_params(A)
        samples = 1000
        dt = self.t/samples

        P = np.zeros((samples + 1, 5))
        T = np.zeros(samples + 1)
        P[0, :] = self.P0
        T[0] = self.t0

        idx = 0
        while rk45.successful() and rk45.t < self.t:
            rk45.integrate(rk45.t+dt)
            P[idx, :] = np.array(rk45.y)
            T[idx] = rk45.t
            idx += 1

        return P, T

solver = NumericalSolver(np.array([1, 0, 0, 0, 0]), 0, 20)

plt.plot(solver.T, solver.P[:, 0], 'r--', linewidth=2.0)
plt.plot(solver.T, solver.P[:, 1], 'b--', linewidth=2.0)
plt.plot(solver.T, solver.P[:, 2], 'g--', linewidth=2.0)
plt.plot(solver.T, solver.P[:, 3], 'y--', linewidth=2.0)
plt.plot(solver.T, solver.P[:, 4], 'c-', linewidth=4.0)
plt.legend(['R', 'AR', 'A2R', 'A2D', 'A2O'])
plt.xlabel('time []')
plt.ylabel('state probability')
plt.show()
