import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sta
import time as tm
import numexpr as ne


class SolverGlp:

    def __init__(self,  model, P0, p, t0, te, suspend):
        """
        Monte Carlo Gillespie solver
        :param A: normalized transition rate matrix with stimulus
        :param P0: starting conditions
        :param t0: starting time
        :param te: ending time
        """
        self.model, self.A, self.P0, self.t0, self.te, self.suspend = model, model.trm, P0, t0, te, suspend
        self.parno = p
        self.stano = len(P0[0])
        print("### Initiating Gillespie Monte Carlo Solver")
        print("Particles: {}, States: {}".format(self.parno, self.stano))
        print("Initial concentrations: \n {}".format(self.P0))
        print("Initial transition matrix: \n{}".format(self.A(0)))

        start_time = tm.time()
        self.allP, self.allT = self.solve_glp()
        print("--- %s seconds ---" % (tm.time() - start_time))
        start_time = tm.time()
        self.mcP, self.mcT = self.get_cumulative()
        print("--- %s seconds ---" % (tm.time() - start_time))


    def solve_glp(self):
        """
        Monte Carlo Gillespie solver
        :return: array of probabilities, array of times
        """

        allP = []
        allT = []

        for particle in range(self.parno):

            P = np.array(self.P0)
            T = np.array([self.t0])

            # P = np.zeros((1, self.stano))
            # T = np.zeros(1)
            # P[0, :] = self.P0[0]
            # T[0] = self.t0

            print("### Gillespie iteration begins!")
            initial = np.where(self.P0 > 0)[0]
            print("Initial state: {}".format(self.model.states[initial].name))

            step = 0
            start_time = tm.time()

            while T[step] < self.te:
                current = np.where(P[step] > 0)
                current_rate = self.A(T[step])[current][0]
                current_rate_sum = ne.evaluate('sum(current_rate)')

                print("Step: {}, time: {}".format(step, T[step]))
                print(" {}".format(P[step]))
                print(" {}".format(current_rate))
                print(self.A(T[step]))

                if current_rate_sum == 0:
                    print("Zero transition rate, staying in previous state")
                    P = np.append(P, [P[step]], axis=0)
                    # dt = 0.01
                    for period in self.suspend:
                        if period[0] <= T[step] < period[1]:
                            dt = period[1] - T[step]
                else:
                    new = np.zeros(self.stano)
                    new_state = self.get_uni(current_rate)
                    print("Non zero transition rate, changing to state {}".format(new_state))
                    new[new_state] = 1.0
                    P = np.append(P, [new], axis=0)
                    dt = self.get_exp(current_rate)

                print("Occupancy time: {}".format(dt))

                T = np.append(T, np.array([T[-1] + dt]))
                step += 1

            else:
                print("### Gillespie iteration ends!")
                print("Final time {} achieved after {} steps".format(self.te, step))

            allP.append(P)
            allT.append(T)

            print("--- %s seconds ---" % (tm.time() - start_time))


        return allP, allT

    def get_cumulative(self):
        """

        :return:
        """
        print("###Cumulative analysis begins!")

        samples = 1000

        sampleT = np.linspace(self.t0, self.te, samples)
        sampleP = np.array([np.zeros(self.stano) for i in range(samples)])
        previousP = np.array([(np.zeros(self.stano)) for i in range(self.parno)])

        for frame, f0 in enumerate(sampleT):

            # abort at last frame
            if frame == len(sampleT) - 1:
                break
            f1 = sampleT[frame+1]

            print("Frame:{} time {} ms to {} ms".format(frame, f0, f1))

            # for each particle
            for particle in range(self.parno):
                single_sampleP = np.zeros(self.stano)
                indexes = np.where(np.logical_and(self.allT[particle] >= f0, self.allT[particle] < f1))
                times = self.allT[particle][indexes]
                states = self.allP[particle][indexes]
                print("Particle: {} Timestep: {} States: {}".format(particle, times, states))

                for state in states:
                    single_sampleP += state

                if np.sum(single_sampleP) == 0.:
                    single_sampleP = previousP[particle]
                previousP[particle] = single_sampleP

                sampleP[frame] += single_sampleP

        print("###Cumulative analysis ends!")
        sum_sampleP = sampleP.sum(axis=1)
        norm_sampleP = sampleP / sum_sampleP[:, np.newaxis]
        return norm_sampleP, sampleT

    def get_distributions(self):
        # TODO: add array to store occupancy times of each state
        # TODO: add array to store cumulative times of each shut and open period
        # TODO: plot histograms of all
        pass

    def get_results(self):
        """
        brings numeric results
        :return: time array, probability array
        """
        return self.allT, self.allP, self.mcT, self.mcP

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
        chooses a transition on the basis on weighted probability by rate value
        :param rates: vector (transition rates)
        :return: scalar (index of random rate from weighted distribution)
        """
        idxs = np.arange(self.stano)
        probs = rates/ne.evaluate('sum(rates)')
        random = np.random.choice(idxs, p=probs)
        return random


