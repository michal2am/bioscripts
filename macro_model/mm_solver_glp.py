import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sta


class SolverGlp:

    def __init__(self,  model, P0, t0, te):
        """
        Monte Carlo Gillespie solver
        :param A: normalized transition rate matrix with stimulus
        :param P0: starting conditions
        :param t0: starting time
        :param te: ending time
        """
        self.model, self.A, self.P0, self.t0, self.te = model, model.trm, P0, t0, te
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
        P = np.zeros((1, self.stano))
        T = np.zeros(1)
        P[0, :] = self.P0
        T[0] = self.t0

        print("### Gillespie iteration begins!")
        initial = np.where(self.P0 > 0)[0]
        print("Initial state: {}".format(self.model.states[initial].name))

        # TODO: add array to store occupancy times of each state
        # TODO: add array to store cumulative times of each shut and open period
        # TODO: plot histograms of all

        step = 0
        while T[step] < self.te:
            current = np.where(P[step] > 0)
            current_rate = self.A(T[step])[current]
            current_rate_sum = np.sum(current_rate)
            dt = self.get_exp(current_rate) if current_rate_sum > 0 else 0.01

            print("Step: {}, time: {}".format(step, T[step]))
            print(" {}".format(P[step]))
            print("Occupancy time: {}".format(dt))
            print(self.A(T[step]))

            if current_rate_sum == 0:
                print("Zero transition rate, staying in previous state")
                P = np.append(P, [P[step]], axis=0)
            else:
                new = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])
                new_state = self.get_uni(current_rate[0])
                print("Non zero transition rate, changing to state {}".format(new_state))
                new[0][new_state] = 1.0
                P = np.append(P, new, axis=0)

            T = np.append(T, np.array([T[-1] + dt]))
            step += 1

        else:
            print("### Gillespie iteration ends!")
            print("Final time {} achieved after {} steps".format(self.te, step))

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


