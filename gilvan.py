import math
import scipy
from scipy.special import expi
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import itertools  


def Whitt(x):
#Returns the value of Whittaker´s function W(1/2,1/2,z)
#Author: G.Feitosa, July 10, 1991, PhD, U.Tulsa (originally in FORTRAN)
#VBA Code: L.Kubota, Feb 2015
#Python Code: L.Kubota, Sept 2019
#Reference: SPE 19845 by Oliver
    if x< 2.0:
        soma=0
        for k in np.arange(0,9,1):
            PSIK1 = PSIN(k + 1)
            PSIK2 = PSIN(k + 2)
            PSIK3 = PSINHALF(k)
            TBRAC = PSIK1 + PSIK2 - PSIK3 - np.log(x)
            GNHALF = GAMMHALF(k)
            FACTK = math.factorial(k)
            FACTK1 = math.factorial(k + 1)
            DEN = FACTK * FACTK1
            soma = soma + GNHALF * (x**k) * TBRAC / DEN
        TLBRAC = soma + 2 * np.sqrt(np.pi) / x
        XMULT = x * np.exp((-x / 2.0)) / (2.0 * np.pi)
        whitt = XMULT * TLBRAC    
    else:
        whitt = np.sqrt(x) * np.exp(-x / 2.0) * (1.0 + (1.0 / (4.0 * x)) - (3.0 / (32.0 * (x**2))) + (15.0 / (128.0 * (x**3))) - (525.0 / (2048.0 * (x**4))))

    return whitt

def GAMMHALF(n):
#Returns the value of gamma(n+1/2)where n is an integer such that n<=8
#Author: G.Feitosa, July 10, 1991, PhD, U.Tulsa (originally in FORTRAN)
#VBA Code: L.Kubota, Feb 2015
#Python Code: L.Kubota, Sept 2019
    g=np.ones(9)
    g[0] = 1.7724538509
    g[1] = 0.88622693
    g[2] = 1.3293404
    g[3] = 3.323351
    g[4] = 11.631728
    g[5] = 52.342778
    g[6] = 287.88528
    g[7] = 1871.2543
    g[8] = 14034.407  
    
    return g[n]
    
def PSIN(n):
#Returns the value of psi(n)
#Author: G.Feitosa, July 10, 1991, PhD, U.Tulsa (originally in FORTRAN)
#VBA Code: L.Kubota, Feb 2015
#Python Code: L.Kubota, Sept 2019

    Euler = 0.577215664901533
    if n==1:
        psin=-Euler
    else:
        soma=0
        for i in np.arange(1,n,1):
            soma=soma+1.0/i
        psin=-Euler+soma    
    return psin

def PSINHALF(n):
#Returns the value of psi(n+1/2)
#Author: G.Feitosa, July 10, 1991, PhD, U.Tulsa (originally in FORTRAN)
#VBA Code: L.Kubota, Feb 2015
#Python Code: L.Kubota, Sept 2019
    phalf = -1.96351002602142
    if n==0:
        psinhalf=phalf
    else:
        soma=0
        for i in np.arange(1,n+1,1):
            soma=soma+1.0/(2.0*i-1.0)
        psinhalf=phalf+2*soma
    
    return psinhalf    

def Kernel(rd,td):
#Returns the value of the Kernel function K1(rd,td)=2*K(rd,td) where K=K(rd,td) is Oliver´s Kernel.
#Author: G.Feitosa, July 10, 1991, PhD, U.Tulsa (originally in FORTRAN)
#VBA Code: L.Kubota, Feb 2015
#Python Code: L.Kubota, Sept 2019
#Reference: G.Feitosa, JPT, July 1994 "Determination of Permeability Distribution From Well-Test Pressure Data"

    kernel = np.sqrt(np.pi) * (rd / td) * np.exp(-(rd * rd) / (2 * td)) * Whitt((rd * rd) / td)

    return kernel

def Integral(rd , td , kr ):
    m = len(rd)
    soma = 0
    for j in np.arange(1,m,1):
        ds = ((Kernel(rd[j], td) * (1 / kr[j])) + (Kernel(rd[j - 1], td) * (1 / kr[j - 1]))) * (rd[j] - rd[j - 1]) * 0.5
        soma = soma + ds
    integral = soma
    return integral









def main():
	
	
 if __name__ == "__main__":
	 main()