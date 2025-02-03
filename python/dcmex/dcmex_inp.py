import numpy as np
import matplotlib.pyplot as plt


if __name__=="__main__":
	print('Calculating the INP parameterisation from Daly et al')
	
	T=np.linspace(273.15,273.15-40,1000)
	
	
	NINP=np.zeros(len(T))
	Tc=T-273.15
	log_ns = -3.25 + (-0.793*Tc) + (-6.91e-2*Tc**2) + (-4.17e-3*Tc**3) + \
		(-1.05e-4*Tc**4) + (-9.08e-7*Tc**5)
	ns = (10.0**log_ns)*1e4
	
	Adust=7.5 # micrometers^2 per cm^-3
	Adust = Adust*1e6 # micrometers ^2 m^-3
	Adust = Adust/(1e12) # m^2 m^-3
	NINP1 = Adust*0.05*ns
	NINP1[:] = 0.0
	ax=plt.plot(Tc,NINP1/1000.)
	plt.ylabel('N$_{INP}$ (L$^{-1}$)')
	plt.xlabel('T ($\circ$C)')
	#ax[0].axes.invert_xaxis()
	
	NINP2=1000.*np.exp(-50.0+45.25*(-4.-Tc)**0.046)
	plt.plot(Tc,NINP2/1000.)
	
	# temperatures less than -23 degC
	#inds,=np.where(T<(273.15-23.0))
	#NINP[inds]=NINP1[inds]
	# temperatures greater than or equal to -23 degC
	#inds,=np.where(T>=(273.15-23.0))
	#NINP[inds]=NINP2[inds]
	NINP=NINP1+NINP2
	# first ice -4
	inds,=np.where(T>=(273.15-4))
	NINP[inds]=0.0
	
	
	plt.plot(Tc,NINP/1000.)
	#NINP[inds]=	
	print('done')
	