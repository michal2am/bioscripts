import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class RateMatrix:

    def __init__(self, rates, borders):
        self.rates = np.array(rates)
        self.num_rates = len(rates)
        self.borders = borders

        self.diagonal()

    def diagonal(self):
        """
        puts row sums on the column
        """
        sums = np.sum(self.rates, axis=0)
        sums *= -1
        self.rates[np.diag_indices_from(self.rates)] = sums

    def derivatives(self):
        def derivative(border, t):
            rates = self.rates

            return [for row_rate in rates]

            return [border[0] * rates[0][0] + border[1] * rates[0][1] + border[2] * rates[0][2] + border[3] * rates[0][3] + border[4] * rates[0][4],
                    border[0] * rates[1][0] + border[1] * rates[1][1] + border[2] * rates[1][2] + border[3] * rates[1][3] + border[4] * rates[1][4],
                    border[0] * rates[2][0] + border[1] * rates[2][1] + border[2] * rates[2][2] + border[3] * rates[2][3] + border[4] * rates[2][4],
                    border[0] * rates[3][0] + border[1] * rates[3][1] + border[2] * rates[3][2] + border[3] * rates[3][3] + border[4] * rates[3][4],
                    border[0] * rates[4][0] + border[1] * rates[4][1] + border[2] * rates[4][2] + border[3] * rates[4][3] + border[4] * rates[4][4]]


        t = np.linspace(0, 10, 1000)
        border = [1, 0, 0, 0, 0]
        C = odeint(derivative, border, t)
        plt.plot(t, C[:, 0], 'r--', linewidth=2.0)
        plt.plot(t, C[:, 1], 'b--', linewidth=2.0)
        plt.plot(t, C[:, 2], 'g--', linewidth=2.0)
        plt.plot(t, C[:, 3], 'y--', linewidth=2.0)
        plt.plot(t, C[:, 4], 'c-', linewidth=4.0)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend(['R', 'AR', 'A2R', 'A2D', 'A2O'])
        plt.show()


test_matrix = [[0.0,   0.5,   0.0,   0.0,   0.0],
               [2*2.0, 0.0,   2*0.5, 0.0,   0.0],
               [0.0,   2.0,   0.0,   0.3,   1.0],
               [0.0,   0.0,   0.1,   0.0,   0.0],
               [0.0,   0.0,   1.5,   0.0,   0.0]]

test_border = [1.0, 0.0, 0.0, 0.0, 0.0]

test_rates = RateMatrix(test_matrix, test_border)
test_rates.derivatives()

def A(t):
    """
    transition rate matrix
    :param t: main timeline (from integrator)
    :return: normalized transition rate matrix with stimulus
    """

    test_matrix = [[-4 * Stimulus.square(2, 5, 10)(t), 2 * 2 * Stimulus.square( 2, 5, 10)(t), 0, 0, 0],
                   [0.5, -0.5 - 2 * Stimulus.square(2, 5, 10)(t), 2 * Stimulus.square(2, 5, 10)(t), 0, 0],
                   [0,  2*0.5,  -2.6,   0.1,  1.5],
                   [0,      0,   0.3,  -0.3,    0],
                   [0,      0,     1,     0,   -1]]

    return np.array(test_matrix)