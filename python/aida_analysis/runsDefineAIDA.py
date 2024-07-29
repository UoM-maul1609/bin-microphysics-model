import numpy as np

namelist_fn='/namelist-aida-ice.in'
NaClMR= np.logspace(-14,-4,50)

w_flag=False
totM1=1000.e6
# totM1=5000.e6
winit=1.3
logSig=0.5
Dm=100.e-9
# see moments of lognormal distribution
N_aer=NaClMR / (np.pi/6.*2165.*np.exp(3.*np.log(Dm)+4.5*logSig**2))
    

N_aer=np.linspace(0.0,10000e6,50)
if w_flag:
	winit=np.linspace(0.01,10,50)

psd_type1=1006 # psd of the background aerosol
psd_type=1 # psd of the added aerosol

if psd_type==1:
	N2=np.array([0.49,0.38])
	logSig2=np.array([0.25,0.84])
	Dm2=np.array([0.247e-6,0.205e-6])
elif psd_type==2:
	N2=np.array([0.18,0.74])
	logSig2=np.array([0.19,0.45])
	Dm2=np.array([0.122e-6,0.140e-6])
elif psd_type==3:
	N2=np.array([0.16,0.91])
	logSig2=np.array([0.19,0.43])
	Dm2=np.array([0.084e-6,0.115e-6])
elif psd_type==4:
	N2=np.array([0.2,1.06])
	logSig2=np.array([0.23,0.47])
	Dm2=np.array([0.061e-6,0.102e-6])
elif psd_type==5:
	N2=np.array([0.6,1.37])
	logSig2=np.array([0.49,0.76])
	Dm2=np.array([0.038e-6,0.08e-6])
elif psd_type==6:
	N2=np.array([0.56,1.118])
	logSig2=np.array([0.46,0.68])
	Dm2=np.array([0.029e-6,0.053e-6])

# psd of the background aerosol	
if psd_type1==1:
	N1=np.array([0.49,0.38])
	logSig1=np.array([0.25,0.84])
	Dm1=np.array([0.247e-6,0.205e-6])
elif psd_type1==2:
	N1=np.array([0.18,0.74])
	logSig1=np.array([0.19,0.45])
	Dm1=np.array([0.122e-6,0.140e-6])
elif psd_type1==3:
	N1=np.array([0.16,0.91])
	logSig1=np.array([0.19,0.43])
	Dm1=np.array([0.084e-6,0.115e-6])
elif psd_type1==4:
	N1=np.array([0.2,1.06])
	logSig1=np.array([0.23,0.47])
	Dm1=np.array([0.061e-6,0.102e-6])
elif psd_type1==5:
	N1=np.array([0.6,1.37])
	logSig1=np.array([0.49,0.76])
	Dm1=np.array([0.038e-6,0.08e-6])
elif psd_type1==6:
	N1=np.array([0.56,1.118])
	logSig1=np.array([0.46,0.68])
	Dm1=np.array([0.029e-6,0.053e-6])
elif psd_type1==1006:
	# dust aerosol
	N1=np.array([197.63,196.55])
	logSig1=np.array([0.67,0.31])
	Dm1=np.array([0.207e-6,0.374e-6])

N2a=N_aer*N2[0]/np.sum(N2)
N2b=N_aer*N2[1]/np.sum(N2)

if w_flag:
	N2a[:]=1.e-8
	N2b[:]=1.e-8

# N1=np.array([0.18,0.74])
N11=totM1*N1/np.sum(N1)


outputDir='/tmp'
