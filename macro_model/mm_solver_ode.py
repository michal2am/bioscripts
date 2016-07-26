import numpy as np
import pandas as pd
import scipy.integrate as itg
import scipy.optimize as opt
import logging as log


class SolverOde:

    def __init__(self, A, P0, names, t0, te, steady):
        """
        differential equation solver
        :param A: normalized transition rate matrix with stimulus
        :param P0: starting conditions (taking first row, if multiple specified)
        :param names: names of states
        :param t0: starting time
        :param te: ending time
        :param steady:
        """
        self.A, self.P0, self.names, self.t0, self.te, self.steady = A, P0, names, t0, te, steady
        self.TP = self.solve_kfw()

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
            if self.steady: t = 0
            return np.dot(P, A(t))

        rk45 = itg.ode(dpdt).set_integrator('lsoda', nsteps=1e4, atol=0.0001)
        rk45.set_initial_value(self.P0, self.t0).set_f_params(self.A)
        samples = 1e4
        dt = self.te / samples

        P = np.zeros((samples + 1, len(self.P0)))
        T = np.zeros(samples + 1)
        P[0, :] = self.P0
        T[0] = self.t0

        idx = 0
        log.info("### Initiating RK ODE Solver")
        while rk45.successful() and rk45.t < self.te:
            rk45.integrate(rk45.t+dt)
            result = np.array(rk45.y)
            log.info("Step: {} Occupancies: {}".format(idx, result))
            P[idx, :] = result
            T[idx] = rk45.t
            idx += 1
        log.info("### Done RK ODE Solver")

        P = P[0:-1]
        T = T[0:-1]

        T = np.array([T]).transpose()
        TP = np.concatenate((T, P), axis=1)
        TP = pd.DataFrame(data=TP, columns=['time'] + self.names)

        return TP

    def get_results(self):
        """
        brings numeric results
        :return: pandas data frame (columns: time and state names)
        """
        return self.TP
