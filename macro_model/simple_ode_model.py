import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def rxn(C, t):
    R = C[0]
    AR = C[1]
    A2R = C[2]
    A2D = C[3]
    A2O = C[4]

    kon = 2.0 * imp(t)
    koff = 0.5
    r = 0.3
    d = 0.1
    a = 1.0
    b = 1.5

    dRdt = -2*kon * R  + koff * AR
    dARdt = 2*kon * R  - koff * AR - kon * AR + 2*koff * A2R
    dA2Rdt = kon * AR - 2*koff * A2R - d * A2R + r * A2D - b * A2R + a * A2O
    dA2Ddt = d * A2R - r * A2D
    dA2Odt = b * A2R - a * A2O

    return [dRdt, dARdt, dA2Rdt, dA2Ddt, dA2Odt]


def imp(t):
    if 2 < t < 5:
        return 10
    else:
        return 0

t =np.linspace(0, 10, 1000)
C0 = [1, 0, 0, 0, 0]
C = odeint(rxn, C0, t)

plt.plot(t, C[:, 0], 'r--', linewidth=2.0)
plt.plot(t, C[:, 1], 'b--', linewidth=2.0)
plt.plot(t, C[:, 2], 'g--', linewidth=2.0)
plt.plot(t, C[:, 3], 'y--', linewidth=2.0)
plt.plot(t, C[:, 4], 'c-', linewidth=4.0)
plt.plot(t, [0.01 * imp(ti) - 0.25 for ti in t], 'k-', linewidth=2.0)
plt.xlabel('x')
plt.ylabel('y')
plt.legend(['R', 'AR', 'A2R', 'A2D', 'A2O'])
plt.show()

