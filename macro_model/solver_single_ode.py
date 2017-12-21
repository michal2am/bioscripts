import numpy as np
import pandas as pd
import scipy.integrate as itg
import logging as log


class SolverOde:

    def __init__(self, a, p0, names, t0, te):
        """
        differential equation solver
        :param a: normalized transition rate matrix with stimulus
        :param p0: starting conditions (taking first row, if multiple specified)
        :param names: names of states
        :param t0: starting time
        :param te: ending time
        """
        self.a, self.p0, self.names, self.t0, self.te = a, p0, names, t0, te
        self.tp = self.solve_kfw()

    def solve_kfw(self):
        """
        differential equation solver
        :return: array of probabilities, array of times
        """

        def dpdt(t, p, a):
            """
            differential equation: dP/dt = A * P
            :param t: blank parameter for integrator and stimulus
            :param p: zero-time + blank to fill matrix
            :param a: transition rate matrix row normalized
            :return: differential equation for integrator
            """
            return np.dot(p, a(t))

        rk45 = itg.ode(dpdt).set_integrator('dopri5', nsteps=1e3, atol=0.001)
        rk45.set_initial_value(self.p0, self.t0).set_f_params(self.a)
        samples = int(1e4)
        dt = self.te / samples

        p = np.zeros((samples + 1, len(self.p0)))
        t = np.zeros(samples + 1)
        p[0, :] = self.p0
        t[0] = self.t0

        idx = 0
        log.info("### Initiating RK ODE Solver")
        while rk45.successful() and rk45.t < self.te:
            rk45.integrate(rk45.t+dt)
            result = np.array(rk45.y)
            log.info("Step: {} Occupancies: {}".format(idx, result))
            p[idx, :] = result
            t[idx] = rk45.t
            idx += 1
        log.info("### Done RK ODE Solver")

        p = p[0:-1]
        t = t[0:-1]

        t = np.array([t]).transpose()
        tp = np.concatenate((t, p), axis=1)
        tp = pd.DataFrame(data=tp, columns=['time'] + self.names)
        tp = tp.set_index('time', drop=True)

        return tp

    def get_results(self):
        """
        brings numeric results
        :return: pandas data frame (time + states)
        """
        return self.tp
