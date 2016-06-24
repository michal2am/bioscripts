import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as pl

def getIp(time):
        if time > 300 and time < 600:
                x = 5
        else:
                x = 0
        return x
def func1(Rates,time):
        # 0 = x, 1 = y
        # dx/dt = -x+ ax+by, dy/dt =-y+ cx+dy
        leak = -0.3
        a = d = 0.1
        b = c = -0.2
        dxbydt = leak*Rates[0]+ a*Rates[0]+b/2.*Rates[1]+getIp(time)
        dybydt = leak*Rates[1]+ c*Rates[0]+d*Rates[1]+getIp(time)
        return [dxbydt, dybydt]

timeGrid = np.arange(0,1000,0.01)
ip = np.zeros((len(timeGrid)))

#ip[300:600] = 5.0
initR = np.ones((2))*10
fR = odeint(func1,initR,timeGrid)

pl.figure()
pl.plot(timeGrid,ip,'k-',label='ip')
pl.plot(timeGrid,fR[:,0],'b-',label='x')
pl.plot(timeGrid,fR[:,1],'r-',label='y')
pl.legend()
pl.show()