import numpy as np
import matplotlib.pyplot as plt



class SolverGlp:

    def __init__(self):
        pass

    def get_exp(self, rates, show=False):
        all_rates = np.sum(rates)
        random = np.random.exponential(1.0/all_rates)
        if show:
            x = np.arange(0.0, 1.5 * random, 0.01)
            pdf = [all_rates * np.exp(-all_rates * xi) for xi in x]
            plt.plot(x, pdf)
            plt.plot(random, 0.5 * np.max(pdf), '*')
            print
        return random

    def get_uni(self, rates):
        pass

test = SolverGlp()
exp = test.get_exp(np.array([3.7, 4.4, 1.4]), show=True)
plt.show()


