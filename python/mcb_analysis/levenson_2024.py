import numpy as np
import scipy.optimize as sco
from scipy.optimize import curve_fit

"""
	parameter constants
	https://doi.org/10.1016/j.pss.2024.105866
"""
sigma_sb=5.67e-8
S0=1362.9
ecc=0.0167 # eccentricity of earth's orbit
a=1.0 # AU distance of earth from sun
S=S0/(a**2*np.sqrt(1.0-ecc**2)) # equation 8
Asurf=0.295 # surface shortwave albedo
epsilon = 0.98 # surface longwave emissivity
alpha=epsilon
Ps=101325. # surface pressure
Patm=101325. # 1 atmosphere in Pa
P=Ps/Patm
#PCH4 = 0.061 # partial pressures
PCH4 = 1.9e-6*Ps
#PCO2 = 28.4
PCO2 = 420e-6*Ps
PH2O = 366
PN2O = 0.028
PO3  = 0.004 #0.039
Fgeo = -0.087 
# feedbacks
feedback_flag=True
Acloud=0.448


def cloud_fit(Ts,a,b):
	return 1.0/(1.0+10**(a*(Ts-b)))

def do_calc(Tguess):
	global PH2O
	
	# mean cloud cover
	n=0.53
	n=0.54
	if feedback_flag:
		PH2O = 402*np.exp(0.0698*(Tguess-288))
		#n=1.0/(1.0+10**(0.0079223*(Tguess-294.69)))
		n=1.0/(1.0+10**(7.43733582e-03*(Tguess-2.97363021e+02)))
	# equation 11 (12-17) for longwave partial optical depths
	# note that total pressure is in atmospheres
	tau_lw_CH4 = 0.0656*PCH4**0.336*P**0.461
	tau_lw_CO2 = 0.105*PCO2**0.357*P**0.268
	tau_lw_H2O = 0.143*PH2O**0.315*P**0.653
	tau_lw_N2O = 0.0973*PN2O**0.336*P**0.461
	tau_lw_O3  = 0.105*PO3**0.336*P**0.221
	tau_lw_cld = 0.491*n**0.112*P**0.957
	tau_lw=tau_lw_CH4+tau_lw_CO2+tau_lw_H2O+tau_lw_N2O+tau_lw_O3+tau_lw_cld
	
	# equation 18-21 - note, seems to be a ratio 1/0.95 out
	frac=1.0 #0.95
	tau_sw_CO2 = 0.00508/frac*PCO2**0.404
	tau_sw_H2O = 0.0551/frac*PH2O**0.315
	tau_sw_O3 = 0.351/frac*PO3**0.336
	tau_sw_cld = 0.0791/frac*n**0.112
	
	tau_sw=tau_sw_CO2+tau_sw_H2O+tau_sw_O3+tau_sw_cld
	
	# calculate consistent with OD, surface illumination
	Fsi = 0.25*S*np.exp(-tau_sw) # or set to 188
	#tau_sw = -np.log(Fsi/(0.25*S))
	
	# now we need to calculate the albedo, see Fig 2.
	# first step is the calculate:
	# 1. fland, Aland, Asea => Abase
	fland=0.292
	fsea=0.708
	Aland=0.2
	Asea=0.045
	Abase = fland*Aland+(1.0-fland)*Asea # equation 22
	# 2. fice, Aice, Abase => Asurf 
	fice=0.061 # fraction of ice under the clouds
	Aice=0.63
	if feedback_flag:
		#pass
		fice=((328.0-Tguess)/95.4)**3.22 # - equation 44
		#fice=((328.0-Tguess)/70)**5 # - equation 44
		
	Asurf = fice*Aice+(1.0-fice)*Abase # equation 23
	# Fsolar - absorbed sunlight - equation 2, basically
	Fsolar=(1.0-Asurf)*Fsi
	Acover=n*Acloud+(1.0-n)*Asurf # equation 24
	C=1.0+1.5*PCO2/Ps # eqaution 30
	
	# rayleigh scattering
	tau_ray=0.062*P**0.858*C # slight edit
	Aray=1.0-np.exp(-tau_ray) # equation 29
	
	Tbase=1.831-3.557*np.exp(-tau_sw)+3.028*np.exp(-tau_sw)**2
	Tbase=np.minimum(1.0,np.maximum(Tbase,0))
	Acc=Acover*Tbase # equation 34
	Asolar=Acc+Aray
	
	D=0.1
	kspec=1.243*D**0.09407
	A=kspec*Asolar
	# Flux density absorbed by the planet's climate system
	F=0.25*S*(1.0-A) # equation 4
	
	Te=(F/sigma_sb)**(0.25)
	Fgreen=3/4*alpha*(sigma_sb*Te**4)*tau_lw # absorbed atmospheric back radiation
	
	k1=0.0999
	k2=1.23
	Fconv=k1*(Fsolar*tau_lw)**k2
	
	
	Frad=Fgeo+Fsolar+Fgreen-Fconv
	
	if feedback_flag:
		Frad=Frad-0.084*(Tguess-288.15) # - equation 46
	
	Ts=(Frad/sigma_sb/epsilon)**0.25
	
	return Ts-Tguess
	
if __name__=="__main__":
	if feedback_flag:
		roots=sco.fsolve(do_calc, 289.)
		Ts=roots[0]
	else:
		Ts=do_calc(0.)
	print(Ts-273.15)
	
# 	popt,pcov = curve_fit(cloud_fit,np.array([288.,252.]),\
# 		np.array([0.54,0.685]),p0=[0.0079223,294.69],method='trf')
	
