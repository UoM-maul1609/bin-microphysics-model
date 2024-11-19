import numpy as np


NaClMR= np.logspace(-14,-4,50)

winit=0.3
logSig=0.25
Dm=162.e-9
# see moments of lognormal distribution
N_aer=NaClMR / (np.pi/6.*2165.*np.exp(3.*np.log(Dm)+4.5*logSig**2))
    

outputDir='/tmp'
