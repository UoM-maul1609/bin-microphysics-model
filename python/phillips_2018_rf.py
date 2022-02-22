import numpy as np
import matplotlib.pyplot as plt

def CalcFrag(D,t):

    Dthresh=np.minimum(D,1.6e-3)
    X=np.log10(Dthresh*1000.)

    # table 3, phillips et al.
    beta1=0. #-0.1839*X*X-0.2017*X-0.0512
    log10zeta = 2.4268*X*X*X+3.3274*X*X+2.0783*X+1.2927
    log10nabla = 0.1242*X*X*X-0.2316*X*X-0.9874*X-0.0827
    T0 = -1.3999*X*X*X - 5.3285*X*X -3.9847*X - 15.0332

    # table 4, phillips et al. 
    zetaB = -0.4651*X*X*X-1.1072*X*X-0.4539*X+0.5137
    nablaB = 28.5888*X*X*X + 49.8504*X*X + 22.4873*X + 8.0481
    TB0 = 13.3588*X*X*X+15.7432*X*X-2.6545*X-18.4875

    Sigma = np.minimum(np.maximum((D-50.e-6)/10.e-6,0.),1.)
    Omega = np.minimum(np.maximum((-3-t)/3.,0.),1.)


    N = Sigma*Omega*(10**log10zeta * (10**log10nabla)**2 / ((t-T0)**2+(10**log10nabla)**2)+beta1*t)

    N=N*D/Dthresh
    NB = np.minimum(Sigma*Omega*(zetaB*nablaB**2/((t-TB0)**2+nablaB**2)),N)
    return(N,NB)
    
    
    
# calculation of rain-drop freezing fragmentation / mode 1
# use for heterogenous nucleation and collisions with less massive ice
ttr=273.15
rhoi=920.
rhow=1000.
DT=10e-6

# change diameter to that of interest
D=0.150e-3 # diameter of freezing drop
#D=2.6e-3 # diameter of freezing drop

m=rhow*np.pi*D**3/6.

t=np.linspace(273.15,243.15,100)-ttr

(N,NB)=CalcFrag(D,t)



plt.ion()
fig=plt.figure()
plt.plot(N,t)
plt.plot(NB,t)
plt.xscale('log')
plt.xlabel('Number of fragments')
plt.ylabel('T ($^\circ$C)')
ax=fig.gca()
ax.invert_yaxis()

mB=1./2.5*m
mT=rhoi*np.pi*DT**3/6.

plt.legend(['Small splinters , D$_S$=' + str(round(1e6*(6.*mT/(np.pi*1000))**(1./3.),0)) + ' $\mu$m', \
    'Large splinters, D$_L$=' + str(round(1e6*(6.*mB/(np.pi*1000))**(1./3.),0)) + ' $\mu$m'])
plt.title('Number of fragments from a D=' +str(round(1000*D,3)) + 'mm freezing drop')