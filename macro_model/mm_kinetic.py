import numpy as np
import matplotlib.pyplot as plt
from mm_solver_ode import SolverOde


class Kinetic:
    """
    kinetic Q-matrix based model
    """
    def __init__(self, states):
        """
        :param states: vector of model states
        """
        self.states = states
        self.trm = self.trm_create()

    class State:
        def __init__(self, no, name, rates, border):
            """
            :param no: state number (0 start convention)
            :param name: state name
            :param rates: vector of transition rates to other states
            :param border: initial probability of the state
            """
            self.no, self.name, self.rates, self.border = no, name, rates, border

        class Rate:
            def __init__(self, name, value, stimulus=False):
                """
                :param name: rate name
                :param value: rate value
                :param stimulus: stimulus if rate is time dependant
                """
                self.name, self.value, self.stimulus = name, value, stimulus

    def trm_create(self):
        """
        creates transition rate matrix (Q-matrix)
        :return: function returning time dependant transition rate matrix
        """
        trm = [[rate for rate in state.rates] for state in self.states]       # rate object matrix
        trm_r = [[rate.value for rate in row]for row in trm]                  # rate value no time matrix

        print("Zero time transition rates (i(row) -> j(column)")
        print([state.name for state in self.states])
        print(np.array(trm_r))

        def trm_f(t):
            trm_t = np.array([[rate.value * rate.stimulus(t) if rate.stimulus else rate.value for rate in row]
                              for row in trm])                                # rate value with stimulus
            trm_t[np.diag_indices_from(trm_t)] = -1 * np.sum(trm_t, axis=1)   # row normalization
            return trm_t

        return trm_f


class Stimulus:
    """
    stimulus shall be defined as static method returning time dependant function returning single scalar for
    given time, notice, that time value is passed by integrator
    """

    @staticmethod
    def square(a, b, v):
        """
        square stimulus
        :param a: step start time
        :param b: step end time
        :param v: step height (rate multiplication)
        :return: time dependant function
        """
        def square_t(t):
            if a < t < b:
                return v
            else:
                return 0
        return square_t


# model definition #

stimulus = Stimulus.square(2, 12, 10)

r_kon   = Kinetic.State.Rate('kon',   2.0, stimulus)
r_2kon  = Kinetic.State.Rate('2kon',  4.0, stimulus)
r_koff  = Kinetic.State.Rate('koff',  0.5)
r_2koff = Kinetic.State.Rate('2koff', 1.0)
r_d     = Kinetic.State.Rate('d',     0.5)
r_r     = Kinetic.State.Rate('r',     0.3)
r_b     = Kinetic.State.Rate('b',     2.5)
r_a     = Kinetic.State.Rate('a',     1.0)

r_0     = Kinetic.State.Rate('block', 0.0)

st_r    = Kinetic.State(0, 'R',   [r_0,     r_2kon,   r_0,    r_0,  r_0],   1)
st_ar   = Kinetic.State(1, 'AR',  [r_koff,  r_0,      r_kon,  r_0,  r_0],   0)
st_a2r  = Kinetic.State(2, 'A2R', [r_0,     r_2koff,  r_0,    r_d,  r_b],   0)
st_a2d  = Kinetic.State(3, 'A2D', [r_0,     r_0,      r_r,    r_0,  r_0],   0)
st_a2o  = Kinetic.State(4, 'A2O', [r_0,     r_0,      r_a,    r_0,  r_0],   0)

jwm = Kinetic([st_r, st_ar, st_a2r, st_a2d, st_a2o])

# ode solving & plot

solver = SolverOde(jwm.trm, np.array([1, 0, 0, 0, 0]), 0, 25)
T, P = solver.get_results()

plt.plot(T, P[:, 0], 'r--', linewidth=2.0)
plt.plot(T, P[:, 1], 'b--', linewidth=2.0)
plt.plot(T, P[:, 2], 'g--', linewidth=2.0)
plt.plot(T, P[:, 3], 'y--', linewidth=2.0)
plt.plot(T, P[:, 4], 'c-', linewidth=4.0)
plt.plot(T, [0.01 * stimulus(ti) - 0.25 for ti in T], 'k-', linewidth=2.0)
plt.legend(['R', 'AR', 'A2R', 'A2D', 'A2O'])
plt.xlabel('time []')
plt.ylabel('state probability')
plt.show()
