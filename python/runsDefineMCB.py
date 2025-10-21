import numpy as np


NaClMR= np.logspace(-14,-4,50)

winit=0.3
logSig=0.25
Dm=162.e-9

N_aer1=np.array([7246000,3132000,49800])*1e6
#N_aer1=np.array([7246000,3132000,0])*1e6
N_aer1=N_aer1/np.sum(N_aer1)
logSig=np.array([0.63345,0.3607466,0.88649])
Dm=np.array([0.0395,0.1162,0.4477])*1e-6




# see moments of lognormal distribution
a=NaClMR/np.sum((N_aer1*np.pi/6.*2165.*np.exp(3.*np.log(Dm)+4.5*logSig**2)))
N_aer=a*np.ones([len(N_aer1),1])*(N_aer1*np.ones([50,1])).transpose()
    

outputDir='/tmp'
