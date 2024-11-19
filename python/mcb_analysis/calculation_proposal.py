import numpy as np
import matplotlib.pyplot as plt

Re=6371e3
# table 3 from Rasch et al. 2024
regions = {'R1': {'lat': [0,30], 'lon': [-150,-110], \
	'Name': 'NE Pacific', 'Acronym': 'NEP'}, \
	'R2': {'lat': [-30,0], 'lon': [-110,-70], \
	'Name': 'SE Pacific', 'Acronym': 'SEP'}, \
	'R3': {'lat': [-30,0], 'lon': [-25,15], \
	'Name': 'SE Atlantic', 'Acronym': 'SEA'}, \
	'R4': {'lat': [30,50], 'lon': [170,-120], \
	'Name': 'N Pacific', 'Acronym': 'NP'}, \
	'R5': {'lat': [-50,-30], 'lon': [-170,-90], \
	'Name': 'S Pacific', 'Acronym': 'SP'}, \
	'R6': {'lat': [60,90], 'lon': [-180,180], \
	'Name': 'Northern Oceans', 'Acronym': 'NO'}}
	
def CalcArea(Re,lons,lat1,lat2):
	dlon=np.diff(lons)
	if dlon<0:
		dlon=360.0+dlon
	dlon=dlon*np.pi/180.0
	return Re**2*dlon*(np.sin(lat2*np.pi/180.)-np.sin(lat1*np.pi/180.))
	
if __name__=="__main__":
	A1=CalcArea(Re,regions['R1']['lon'], \
		regions['R1']['lat'][0],regions['R1']['lat'][1])
		
	A2=CalcArea(Re,regions['R2']['lon'], \
		regions['R2']['lat'][0],regions['R2']['lat'][1])
		
	A3=CalcArea(Re,regions['R3']['lon'], \
		regions['R3']['lat'][0],regions['R3']['lat'][1])
		
	A4=CalcArea(Re,regions['R4']['lon'], \
		regions['R4']['lat'][0],regions['R4']['lat'][1])
		
	A5=CalcArea(Re,regions['R5']['lon'], \
		regions['R5']['lat'][0],regions['R5']['lat'][1])
		
	A6=CalcArea(Re,regions['R6']['lon'], \
		regions['R6']['lat'][0],regions['R6']['lat'][1])
	
	AT=4.0 *np.pi *Re**2 * 0.71
	print(A1/AT,A2/AT,A3/AT,A4/AT,A5/AT,A6/AT)
	
	# PSD parameters:
	lnSigma=0.5
	Dm=80e-9
	rhoAer=2165
	
	M3=np.exp(3*np.log(Dm)+4.5*lnSigma**2)
	
	burden=50e12*1e-3
	# number of aerosol generated per year
	numAer = burden / (rhoAer*M3*np.pi/6.0)
	# per day
	numAer = numAer / (365)
	# if the lifetime is one day, then divide by the volume of air
	# per cc
	print(numAer / (A2*1000.0) / 1e6)
	
	
	"""
		difference between Koehler and K-Koehler for NaCl
	"""
	kappa = 0.91 # 1.28
	T=293.15
	nu=2.0
	Mw=18e-3
	Ms=58.44e-3
	rhos=2165.0
	rhow=1000.0
	ddry=60.e-9
	mdry=np.pi/6.0*rhos*ddry**3
	ns=mdry/Ms
	nw=np.logspace(np.log10(0.01*ns),np.log10(1e5*ns),1000)
	sigma_t=0.072
	R_gas=8.314
	# mixing rule
	mwat = nw*Mw
	mtot = mdry + mwat
	vtot = mdry/rhos+mwat/rhow
	rhosol = (mtot)/vtot 
	dw = (6.0/np.pi/rhosol*mtot)**(1.0/3.0)
	rh_eq=nw/(ns*nu+nw)*np.exp(4.0*Mw*sigma_t/(R_gas*T*rhow*dw))
	# kappa koehler theory
	aw_inv = 1.0 + kappa*(mdry/rhos) / (mwat/rhow)
	rh_eq_kappa=(1.0/aw_inv) * \
		np.exp(4.0*Mw*sigma_t/(R_gas*T*rhow*dw))
# 	rh_eq_kappa=(dw**3-ddry**3)/(dw**3-(ddry**3)*(1.0-kappa)) * \
# 		np.exp(4.0*Mw*sigma_t/(R_gas*T*rhow*dw))

	plt.ion()
	ax = [plt.subplot(2,1,i+1) for i in range(2)]
	lims = [[1.0,1.01],[0,1]]
	for i in range(2):
		ax[i].plot(dw,rh_eq)
		ax[i].plot(dw,rh_eq_kappa,'--')
		ax[i].set_ylim(lims[i])
		arr=['{:1.4f}'.format(ax[i].get_yticks()[j]) \
			for j in range(len(ax[i].get_yticks()))]
		ax[i].set_yticks(ax[i].get_yticks())
		ax[i].set_yticklabels(arr)
	
	plt.subplots_adjust(wspace=0, hspace=0)