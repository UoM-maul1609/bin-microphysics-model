import numpy as np

nbins=50
NaClMR= np.logspace(-14,-4,nbins)

winit=0.3
logSig=0.25
Dm=162.e-9

COOPER=1
HARRISON_EFFERVESCENCE=2
HARRISON_DELAVAL=3
EDMUND_SCF=4
RAYLEIGH=5

spray_method=RAYLEIGH

# N_aer1=np.array([7246000,3132000,49800])*1e6
# #N_aer1=np.array([7246000,3132000,0])*1e6
# logSig=np.array([0.63345,0.3607466,0.88649])
# Dm=np.array([0.0395,0.1162,0.4477])*1e-6

if spray_method==COOPER:
	N_aer1		=np.array([1.55,194,17.7])*1e6
	logSig		=np.array([0.129,0.625,0.666])
	Dm			=np.array([0.0263,0.0588,0.269])*1e-6
elif spray_method==HARRISON_EFFERVESCENCE:
	N_aer1		=np.array([94800,245000,137000])*1e6
	logSig		=np.array([0.252,0.497,0.883])
	Dm			=np.array([0.0264,0.0434,0.0574])*1e-6
elif spray_method==HARRISON_DELAVAL:
	N_aer1		=np.array([214000,66800,4800])*1e6
	logSig		=np.array([0.703,0.534,0.9])
	Dm			=np.array([0.0156,0.123,0.398])*1e-6
elif spray_method==EDMUND_SCF:
	N_aer1		=np.array([6460000,3830000,32700])*1e6
	logSig		=np.array([0.571,0.391,0.739])
	Dm			=np.array([0.0366,0.109,0.696])*1e-6
elif spray_method==RAYLEIGH:
	N_aer1		=np.array([1.0])*1e6
	logSig		=np.array([0.25])
	Dm			=np.array([0.162])*1e-6

N_aer1=N_aer1/np.sum(N_aer1)




# see moments of lognormal distribution
a=NaClMR/np.sum((N_aer1*np.pi/6.*2165.*np.exp(3.*np.log(Dm)+4.5*logSig**2)))
N_aer=a*np.ones([len(N_aer1),1])*(N_aer1*np.ones([nbins,1])).transpose()
    

outputDir='/tmp'
