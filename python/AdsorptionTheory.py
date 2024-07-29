# adsorption theory from https://acp.copernicus.org/articles/11/3527/2011/acp-11-3527-2011.pdf

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def mean_free_path(T,P,dk):
	return kb*T/(np.sqrt(2.0)*np.pi*dk**2*P)

def calc_X(frac,Xroot,Ddry,lam1):

	
	Dve=frac[0]*Ddry
	C1=Cslip(Ddry,lam1)
	C2=Cslip(Dve,lam1)
	X=Ddry*C2/(Dve*C1)
	print (X,Xroot,frac[0])
	return X-Xroot

def Cslip(Di,lam1):
	return 1.0+2.0*lam1/Di*(1.142+0.558*np.exp(-0.999*Di*0.5/lam1))
if __name__=="__main__":

	P=100000.0
	kb=1.380649e-23
	T=270.0
	dk=364e-12
	
	AFHH=2.96
	BFHH=1.28
# 	AFHH=1.36
# 	BFHH=1.12
	Mw=18e-3
	Ms=5000e-3
	rhos=1500.0
	sigma=0.072
	R_gas=8.314
	Dh2o=2.75e-10
	rhow=1000.0
	Ddry=400e-9
	X_exp=-0.82
	XFHH =-0.83
	X=1.3
	
	lam1=mean_free_path(T,P,dk)
	ns=np.pi/6.0*rhos*Ddry**3/Ms
	nw=np.logspace(-16,-10,10000)
	rhosol=(nw*Mw+ns*Ms)/(nw*Mw/rhow+ns*Ms/rhos)
	d=(6.0*(nw*Mw+ns*Ms)/np.pi/rhosol)**(1.0/3.0)
	
	Dve=fsolve(calc_X,[1.1],args=(X,Ddry,lam1))
	Dve=Dve*Ddry
	Dse=(3.0*X*Dve-Ddry)*0.5
	
	plt.plot(d, \
		np.exp(4.0*sigma*Mw/(R_gas*T*rhow*d) - \
		AFHH*((d-Ddry)/(2.*Dh2o))**-BFHH) )
	plt.plot(d, \
		np.exp(4.0*sigma*Mw/(R_gas*T*rhow*d) - \
		AFHH*((d-Dse)/(2.*Dh2o))**-BFHH) )
		
	kappa=1.28	
	Ddry=100.e-9
	plt.plot(d,np.exp(4*sigma*Mw/(R_gas*T*rhow*d))*(d**3-Ddry**3)/(d**3-Ddry**3*(1.-kappa)) )
		