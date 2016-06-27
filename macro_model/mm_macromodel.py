import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


def A():
    '''
    Transition probability matrix
    '''
    test_matrix = [[-4,    2*2,    0,     0,    0],
                   [0.5,  -2.5,    2,     0,    0],
                   [0,  2*0.5,  -2.6,   0.1,  1.5],
                   [0,      0,   0.3,  -0.3,    0],
                   [0,      0,     1,     0,   -1]]

    return np.array(test_matrix)

def dPdt(t,P,A):
    '''
    Differential equation: dP/dt = A * P
    '''
    return np.dot(P,A())

def solveKFW(Tf):
    '''
    Discretize and integrate the Kolmogorov
    Forward Equation.
    '''

    #Define initial conditions
    P0 = np.array([1, 0, 0, 0, 0])
    t0 = 0
    dt = Tf/1000

    #set the integrator
    rk45 = ode(dPdt).set_integrator('dopri5')
    rk45.set_initial_value(P0,t0).set_f_params(A)

    #integrate and store each solution in
    #a numpy matrix, P and T.
    P = np.zeros((1000000,5))
    T = np.zeros(1000000)
    P[0,:] = P0
    T[0]   = t0

    idx    = 0
    while rk45.successful() and rk45.t < Tf:
        rk45.integrate(rk45.t+dt)
        P[idx,:] = np.array(rk45.y)
        T[idx]   = rk45.t
        idx     += 1

    return P,T

if __name__=='__main__':
    P,T = solveKFW(5)

    plt.plot(T,P[:,0],'-r')
    plt.plot(T,P[:,1],'-g')
    plt.plot(T,P[:,2],'-k')
    plt.plot(T,P[:,3],'-k')
    plt.plot(T,P[:,4],'-k')
    plt.legend(['R', 'AR', 'A2R', 'A2D', 'A2O'])
    plt.xlabel('Time (m)')
    plt.ylabel('Pr of state, s')
    plt.show()
