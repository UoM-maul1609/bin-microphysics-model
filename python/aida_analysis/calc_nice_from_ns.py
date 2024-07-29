import numpy as np
from scipy.integrate import quad


def ns_integrand(x,N,lnsig,dm,ns):
	dNdD=0.0
	for i in range(len(dm)):
		dNdD=dNdD+ \
			N[i]/(x*np.sqrt(2.0*np.pi)*lnsig[i])* \
			np.exp(-(np.log(x/dm[i])**2.0)/(2*lnsig[i]**2))
	return dNdD*(1.0-np.exp(-ns*x**2*np.pi)) #*np.pi*x**2
	
	
if __name__=="__main__":
	N=np.array([197.63289436422826, 196.5507641307173])
	
	N=1000e6*N/np.sum(N)
	lnsig=np.array([0.6683776379438333, 0.3143097539903099])
	dm=np.array([0.20704688461788912e-6, 0.3736013880519923e-6])
	
	t=273.15-24.4
	ns=np.exp(-1.6*(t-273.15)-16.0067)
	
	#ns=1.e30
	result=quad(ns_integrand,1e-9,100e-6,args=(N,lnsig,dm,ns))
	print(result[0]/1000.0)

