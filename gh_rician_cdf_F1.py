import numpy as np
import math
from mpmath import *

def F1(x,xi,wi,K,omega,m_z,v_z,N):
    # Gauss-Hermite f(x) function for hybrid approx CDF with lognormal
    # x is a single value
    mp.dps = 32 #Set precision
    eta = math.log(10)/10
    b = (K+1)/omega
    temp = mpf(0)
    

    np_exp = np.frompyfunc(mp.exp,1,1)
    z = np_exp(math.sqrt(2)*eta*v_z*np.array(xi) + eta*m_z)
    #print(z)
    for i in range(0,int(N)+1):
        for j in range(0,i+1):
            E_ln = mpf(0)
            for k in range(0,len(wi)):
                E_ln = E_ln + (wi[k]/(mpf(z[k])**mpf(j)) * mp.exp(-b*x/z[k]))
                
            E_ln = E_ln / mpf(math.sqrt(math.pi))
            #print(E_ln)
            temp = temp + (((mpf(K)**mpf(i))*(mpf(b)**mpf(j)))/(mp.factorial(i)*mp.factorial(j))) * mp.exp(-K) * mpf(x)**mpf(j) * E_ln
            
    return float(1.0-temp)
                
