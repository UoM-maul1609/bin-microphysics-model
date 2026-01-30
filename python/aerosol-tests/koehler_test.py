import numpy as np
import matplotlib.pyplot as plt






if __name__=="__main__":
	print('Testing')
	
	Mw=18e-3;
	Ms=58.44e-3
	nu=2.
	rho_h2o=1000.
	rho_nacl=2160.
	sig1=0.072
	Rgas=8.314
	T=290.0
	dDry=30e-9
	mDry=np.pi/6*rho_nacl*dDry**3
	
	
	dWetMax=1e-5;
	mWetMax=np.pi/6*rho_h2o*dWetMax**3
	mWat=np.logspace(np.log10(mDry/10),np.log10(mWetMax),1000)
	
	
	nw=mWat/Mw
	ns=mDry/Ms
	
	# total volume adds
	# mwet/rho = msol/rhos+mwat/rhow
	
	rho_tot=(mDry+mWat)/(mDry/rho_nacl+mWat/rho_h2o)
	
	dWet=((mWat+mDry)*6./(np.pi*rho_tot))**(1/3)
	
	RHeq1=nw/(nw+nu*ns)
	RHeq = RHeq1*np.exp(4*Mw*sig1/(Rgas*T*rho_h2o*dWet))
	
	plt.ion()
	plt.plot(dWet,RHeq1)
	plt.plot(dWet,RHeq)
	
	plt.xlabel('$d_w$')
	plt.ylabel('$RH_{eq}$')
	plt.legend(['no Kelvin term','Kelvin term'])