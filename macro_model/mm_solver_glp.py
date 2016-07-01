import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sta


class SolverGlp:

    def __init__(self,  A, P0, t0, te):
        """
        Monte Carlo Gillespie solver
        :param A: normalized transition rate matrix with stimulus
        :param P0: starting conditions
        :param t0: starting time
        :param te: ending time
        """
        self.A, self.P0, self.t0, self.te = A, P0, t0, te
        self.parno = np.sum(P0)
        self.stano = len(P0)
        print("### Initiating Gillespie Monte Carlo Solver")
        print("Particles: {}, States: {}".format(self.parno, self.stano))
        print("Initial concentrations: \n {}".format(self.P0))
        print("Initial transition matrix: \n{}".format(self.A(0)))

        self.P, self.T = self.solve_glp()

    def solve_glp(self):
        """
        Monte Carlo Gillespie solver
        :return: array of probabilities, array of times
        """
        P = np.zeros((self.te + 1, self.stano))
        T = np.zeros(self.te + 1)
        P[0, :] = self.P0
        T[0] = self.t0

        print("### Gillespie iteration begins!")

        for step in range(self.t0, self.te):
            print("Step: {}, time: {}".format(step, T[step]))
            print(" {}".format(P[step]))
            print(self.A(step))

            current = np.where(P[step] > 0)

            current_rate = np.sum(self.A(T[step])[current])
            if current_rate == 0:
                print("Zero transition rate, staying in previous state")
                P[step + 1, :] = P[step]
            else:
                P[step + 1, :] = [0, 0, 0, 0, 0]
                new_state = self.get_uni((self.A(T[step])[current][0]))
                print("Non zero transition rate, changing to state {}".format(new_state))
                P[step + 1, new_state] = 1

            T[step + 1] = T[step] + 0.1

        return P, T

    def get_results(self):
        """
        brings numeric results
        :return: time array, probability array
        """
        return self.T, self.P

    def get_exp(self, rates, show=False):
        """
        generates a random number from exponential distribution of lambda equal to sum of give transition rates
        :param rates: vector (transition rates)
        :param show: plot exponential pdf with random value
        :return: scalar (random from exponential distribution)
        """
        all_rates = np.sum(rates)
        random = np.random.exponential(1.0/all_rates)
        if show:
            x = np.linspace(0.0, 0.01 * all_rates, num=100 * all_rates)            # magic: sampling depending on lambda
            pdf = np.array([all_rates * np.exp(-all_rates * xi) for xi in x])
            limit = np.where(pdf > 0.00001 * all_rates * np.max(pdf))              # magic: omitting the long tail
            plt.plot(x[limit], pdf[limit])
            plt.plot(random, 0.5 * np.max(pdf), '*')
            plt.show()
        return random

    def get_uni(self, rates):
        """
        chooses o transition on the basis on weighted probability by rate value
        :param rates: vector (transition rates)
        :return: scalar (index of random rate from weighted distribution)
        """
        idxs = [i for i in range(len(rates))]
        probs = [rate/float(np.sum(rates)) for rate in rates]
        random = np.random.choice(idxs, p=probs)
        return random


