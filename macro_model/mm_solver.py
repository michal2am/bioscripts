import numpy as np
from scipy.integrate import ode


class Solver:

    def __init__(self, A, P0, t0, te):
        """
        differential equation solver
        :param A: normalized transition rate matrix with stimulus
        :param P0: starting conditions
        :param t0: starting time
        :param te: ending time
        """
        self.A, self.P0, self.t0, self.te = A, P0, t0, te
        self.P, self.T = self.solve_kfw()

    def solve_kfw(self):
        """
        differential equation solver
        :return: array of probabilities, array of times
        """

        def dpdt(t, P, A):
            """
            differential equation: dP/dt = A * P
            :param t: blank parameter for integrator and stimulus
            :param P: zero-time + blank to fill matrix
            :param A: transition rate matrix row normalized
            :return: differential equation for integrator
            """
            return np.dot(P, A(t))

        rk45 = ode(dpdt).set_integrator('dopri5')
        rk45.set_initial_value(self.P0, self.t0).set_f_params(self.A)
        samples = 1000
        dt = self.te / samples

        P = np.zeros((samples + 1, 5))
        T = np.zeros(samples + 1)
        P[0, :] = self.P0
        T[0] = self.t0

        idx = 0
        while rk45.successful() and rk45.t < self.te:
            rk45.integrate(rk45.t+dt)
            P[idx, :] = np.array(rk45.y)
            T[idx] = rk45.t
            idx += 1

        return P, T

    def get_results(self):
        """
        brings numeric results
        :return: time array, probability array
        """
        return self.T, self.P
