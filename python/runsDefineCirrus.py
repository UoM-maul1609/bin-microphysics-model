import numpy as np



winit=np.logspace(-2,1,50)
logSig=0.47
Dm=110.e-9

tinit=235

# see moments of lognormal distribution
massLoading=1.2e-9

# Number per cc
N_aer=massLoading / (np.pi/6.*1830.0*np.exp(3.*np.log(Dm)+4.5*logSig**2)) 

outputDir='/tmp'
