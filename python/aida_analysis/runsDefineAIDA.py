import numpy as np

bmm_run=True
bam_run=True


bam_location='../../../bam/'


namelist_fn='/namelist-aida.in'
NaClMR= np.logspace(-14,-4,50)

"""
	AIDA - totM1=3000e6; winit=1.3
	psd_type1=2
	psd_type=1
	
	more natural - totM1=800e6; w=0.3
"""

w_flag=False
totM1=800.e6
totM1=500.e6
# totM1=5000.e6
winit=0.2
logSig=0.5
Dm=100.e-9
# see moments of lognormal distribution
N_aer=NaClMR / (np.pi/6.*2165.*np.exp(3.*np.log(Dm)+4.5*logSig**2))
    
rhoa=100000./280./287.0
kappa_back=0.61
kappa_add=1.28
#kappa_add=0.5
density_back=1770.
density_add=2165.
#density_add=1770.
density_w=1000.

# equivalent van't hoff factors
mole_w=18e-3
mole_back=132.14e-3
mole_add=58.44e-3
#mole_add=132.14e-3

nu_back = kappa_back*density_w/density_back*mole_back/mole_w
nu_add = kappa_add*density_w/density_add*mole_add/mole_w

N_aer=np.linspace(0.0,1000e6,50)
if w_flag:
	winit=np.linspace(0.01,10,50)

psd_type1=2 # psd of the background aerosol
psd_type=6 # psd of the added aerosol

if psd_type==1:
	N2=np.array([0.49,0.38,1e-8])
	logSig2=np.array([0.25,0.84,0.25])
	Dm2=np.array([0.247e-6,0.205e-6,100e-9])
elif psd_type==2:
	N2=np.array([0.18,0.74,1e-8])
	logSig2=np.array([0.19,0.45,0.25])
	Dm2=np.array([0.122e-6,0.140e-6,100e-9])
elif psd_type==3:
	N2=np.array([0.16,0.91,1e-8])
	logSig2=np.array([0.19,0.43,0.25])
	Dm2=np.array([0.084e-6,0.115e-6,100e-9])
elif psd_type==4:
	N2=np.array([0.2,1.06,1e-8])
	logSig2=np.array([0.23,0.47,0.25])
	Dm2=np.array([0.061e-6,0.102e-6,100e-9])
elif psd_type==5:
	N2=np.array([0.6,1.37,1e-8])
	logSig2=np.array([0.49,0.76,0.25])
	Dm2=np.array([0.038e-6,0.08e-6,100e-9])
elif psd_type==6:
	N2=np.array([0.56,1.118,1e-8])
	logSig2=np.array([0.46,0.68,0.25])
	Dm2=np.array([0.029e-6,0.053e-6,100e-9])
elif psd_type==7:
	N2=np.array([1.0,1.e-8,1e-8])
	logSig2=np.array([0.25,0.68,0.25])
	Dm2=np.array([0.162e-6,0.053e-6,100e-9])
elif psd_type==8:
	N2=np.array([1.0,1.e-8,1e-8])
	logSig2=np.array([0.79,0.54,0.25])
	Dm2=np.array([0.072e-6,0.019e-6,100e-9])

# psd of the background aerosol	
if psd_type1==1:
	N1=np.array([0.49,0.38,1e-8])
	logSig1=np.array([0.25,0.84,0.25])
	Dm1=np.array([0.247e-6,0.205e-6,100e-9])
elif psd_type1==2:
	N1=np.array([0.18,0.74,1e-8])
	logSig1=np.array([0.19,0.45,0.25])
	Dm1=np.array([0.122e-6,0.140e-6,100e-9])
elif psd_type1==3:
	N1=np.array([0.16,0.91,1e-8])
	logSig1=np.array([0.19,0.43,0.25])
	Dm1=np.array([0.084e-6,0.115e-6,100e-9])
elif psd_type1==4:
	N1=np.array([0.2,1.06,1e-8])
	logSig1=np.array([0.23,0.47,0.25])
	Dm1=np.array([0.061e-6,0.102e-6,100e-9])
elif psd_type1==5:
	N1=np.array([0.6,1.37,1e-8])
	logSig1=np.array([0.49,0.76,0.25])
	Dm1=np.array([0.038e-6,0.08e-6,100e-9])
elif psd_type1==6:
	N1=np.array([0.56,1.118,1e-8])
	logSig1=np.array([0.46,0.68,0.25])
	Dm1=np.array([0.029e-6,0.053e-6,100e-9])
elif psd_type1==1006:
	# dust aerosol
	N1=np.array([197.63,196.55,1e-8])
	logSig1=np.array([0.67,0.31,0.25])
	Dm1=np.array([0.207e-6,0.374e-6,100e-9])
elif psd_type1==7:
	# dust aerosol
	N1=np.array([46.64,153.42,166.77])
	logSig1=np.array([0.348,0.354,0.465])
	Dm1=np.array([0.018e-6,0.039e-6,0.154e-6])

N2a=np.maximum(N_aer*N2[0]/np.sum(N2)/rhoa,1e-8)
N2b=np.maximum(N_aer*N2[1]/np.sum(N2)/rhoa,1e-8)
N2c=np.maximum(N_aer*N2[0]/np.sum(N2)/rhoa,1e-8)

if w_flag:
	N2a[:]=1.e-8
	N2b[:]=1.e-8
	N2c[:]=1.e-8

# N1=np.array([0.18,0.74])
N11=totM1*N1/np.sum(N1)/rhoa


outputDir='/tmp'
