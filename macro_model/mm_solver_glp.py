import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sta
# from mm_kinetic import Kinetic



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
        print("### Initiating Gillspie Monte Carlo Solver")
        print("Particles: {}, States: {}, Initial states: {}".format(self.parno, self.stano, self.P0))
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

        print(self.A(0))

        for step in range(self.t0, self.te):
            print(self.A(step))

        print(P)
        print(T)

        return 0, 0

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


