import numpy as np
import matplotlib.pyplot as plt

def Equation13(K0,T,D1,D2,PHI,flag):
    # break-up of new fragments

    DS=np.minimum(D1,D2)
    D=DS
    a = np.pi*DS**2
    
    if (flag == 1):
        # collisions of graupel with graupel / hail - rime mass frac >0.5 and <0.9
        a0=3.78e4*(1.+0.0079/D**1.5)
        T0=-15.
        A=a0/3. + np.maximum(2*a0/3.-a0/9.*(T-ttr-T0),0.)
        phi=3.5e-3
        C=6.30e6*phi
        Nmax=100
        gamma=0.3
        zeta=0.001
    elif (flag == 2):
        # collisions between hail - maybe rime mass frac > 0.9
        a0=4.35e5
        T0=-15.
        A=a0/3. + np.maximum(2*a0/3.-a0/9.*(T-ttr-T0),0.)
        phi=3.5e-3
        C=3.31e5
        Nmax=1000
        gamma=0.54    
        zeta=1e-6
    elif (flag == 3):
        # collisions between ice (low rime fraction - dendrites)
        A=1.41e6*(1.+100.*PHI**2)*(1.+3.98e-5/D**1.5)
        phi=3.5e-3
        C=3.09e6*phi
        Nmax=100
        gamma=0.50-0.25*PHI    
        zeta=0.001
    elif (flag == 4):
        # collisions between ice and spatial planar (low rime fraction)
        A=1.58e7*(1.+100.*PHI**2)*(1.+1.33e-4/D**1.5)
        phi=3.5e-3
        C=7.08e6*phi
        Nmax=100
        gamma=0.50-0.25*PHI    
        zeta=0.001
    
    
    N=np.minimum(a*A*(1.-np.exp(-(C*K0/(a*A))**gamma)), Nmax)
    
    return(N)
    
    
    
# calculation of rain-drop freezing fragmentation / mode 1
# use for heterogenous nucleation and collisions with less massive ice
ttr=273.15


T=258.0
rhoi1=900.
rhoi2=800.
DS=50e-6 # size of small - used to calculate surface area
DL=5e-3 # size of large

MS=rhoi1*np.pi/6.*DS**3
ML=rhoi2*np.pi/6.*DL**3

VS=0.1
VL=5.
PHI=0.1 # rime fraction

K0=0.5*(MS*ML)/(MS+ML)*(VS-VL)**2 

K0=np.logspace(-12,-1,100) # CKE

N=Equation13(K0,T,DS,DL,PHI,4)


plt.ion()
plt.plot(K0,N)
# 
# plt.ion()
# fig=plt.figure()

