def tau_vis(tau):
   if (tau < 0.723):
      return 0.
   else:
      return 0.36*(tau-0.723)**0.411

import numpy as np

pres=.75e5    # pressure to convert mixing ratios to partial pressures

mr_co2=800.  # mixing ratio of co2 (ppmv)
#mr_co2=0.1e6
mr_ch3=1.75  # mixing ratio of methane (ppmv)
mr_h2o=2700. # mixing ratio of water vapour (ppmv)

Ap=0.30        # planetary albedo
sigma=5.67e-8  # stefan-boltzmann constant
Sflux=1368.    # solar flux (W m-2)
epsil_s=0.95   # infrared emissivity of the surface

e_co2=mr_co2/1.e6*pres # conversion from total pressure to partial pressure
e_ch3=mr_ch3/1.e6*pres
e_h2o=mr_h2o/1.e6*pres

# Power law to convert to optical depth
(tau_co2,tau_ch3,tau_h2o)=(0.029,0.725,0.097)*np.sqrt((e_co2,e_ch3,e_h2o))

tau=tau_co2+tau_ch3+tau_h2o # total infrared optical depth

epsil=1./(1.+0.75*tau)      # IR emissivity of atmosphere

T_0=( ((1.-Ap)*Sflux/4.) /(sigma*epsil))**0.25 # raw greenhouse temperature

F=((1.-Ap)*Sflux/4.)
F0=epsil_s*sigma*T_0**4

tvis=tau_vis(tau)

Labs=F-F*np.exp(-tvis)

Fsi=F-Labs
Fabs=(1.-0.15)*Fsi+epsil_s*(F0-F)
Fc=0.369*Fabs*tau/(-0.6+2.*tau)

Fs=F0-Labs-Fc

T_s=(Fs/(epsil_s*sigma))**0.25  # surface temperature

print ('The surface temperature is %.2f K' % T_s)



