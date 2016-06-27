import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

def A():
    '''
    Transition probability matrix
    '''
    return np.array([[-.1,.1,0], [.1,-.2,.1], [0,.1,-.1]])

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
    P0 = np.array([1,0,0])
    t0 = 0
    dt = Tf/100

    #set the integrator
    rk45 = ode(dPdt).set_integrator('dopri45')
    rk45.set_initial_value(P0,t0).set_f_params(A)

    #integrate and store each solution in
    #a numpy matrix, P and T.
    P = np.zeros((100,3))
    T = np.zeros(100)
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
    P,T = solveKFW(120)

    plt.plot(T,P[:,0],'s-r')
    plt.plot(T,P[:,1],'^-g')
    plt.plot(T,P[:,2],'*-k')
    plt.legend(['Inactive A','Active, B', 'Inactive, C'])
    plt.xlabel('Time (m)')
    plt.ylabel('Pr of state, s')
    plt.show()
