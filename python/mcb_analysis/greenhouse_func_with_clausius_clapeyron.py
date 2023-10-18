import numpy as np
import scipy.optimize as sco


def tau_vis(tau):
   if (tau < 0.723):
      return 0.
   else:
      return 0.36*(tau-0.723)**0.411


def surface_temp(Tguess,mr_co2,Sflux,cc_flag):
   pres=0.8e5    # pressure to convert mixing ratios to partial pressures

   #mr_co2=400.  # mixing ratio of co2 (ppmv)
   #mr_co2=0.1e6
   mr_ch3=1.75  # mixing ratio of methane (ppmv)
   #mr_h2o=2700. # mixing ratio of water vapour (ppmv)

   if cc_flag:
      mr_h2o=0.2*610.*np.exp(2.5e6/461.*(-1./Tguess[0]+1./273.15))/pres*1e6
   else:
      mr_h2o=0.2*610.*np.exp(2.5e6/461.*(-1./289.511485670835+1./273.15))/pres*1e6
      
   Ap=0.306        # planetary albedo
   sigma=5.67e-8  # stefan-boltzmann constant
   #Sflux=1368.    # solar flux (W m-2)
   epsil_s=0.95   # infrared emissivity of the surface

   e_co2=mr_co2/1.e6*pres # conversion from total pressure to partial pressure
   e_ch3=mr_ch3/1.e6*pres
   e_h2o=mr_h2o/1.e6*pres

   # Power law to convert to optical depth
   (tau_co2,tau_ch3,tau_h2o)=(0.029,0.725,0.087)*np.sqrt((e_co2,e_ch3,e_h2o))

   #print("%e %e %e" % (tau_co2,tau_ch3,tau_h2o,))
   tau=tau_co2+tau_ch3+tau_h2o # total infrared optical depth

   epsil=1./(1.+0.75*tau)      # IR emissivity of atmosphere

   T_0=( ((1.-Ap)*Sflux/4.) /(sigma*epsil/epsil_s))**0.25 # raw greenhouse temperature

   F=((1.-Ap)*Sflux/4.)
   F0=epsil_s*sigma*T_0**4

   print("tau: %e" % tau)
   tvis=tau_vis(tau)

   

   Labs=F-F*np.exp(-tvis)

   #Labs=F*tvis

   Fsi=F-Labs
   Fabs=(1.-0.15)*Fsi+epsil_s*(F0-F)
   Fc=0.369*Fabs*tau/(-0.6+2.*tau)

#   Fc=93.*1.1 # replace with what the actual value is

   Fs=F0-Labs-Fc
   T_s=(Fs/(epsil_s*sigma))**0.25  # surface temperature

   print("%e %e %e %e" % (Labs,mr_h2o,Fabs,Fc))
   return T_s-Tguess

def run_model(mr_co2=400.,Sflux=1368.,cc_flag=False):
   roots=sco.fsolve(surface_temp, 289., args=(mr_co2,Sflux,cc_flag))
   root=roots[0]
   print("root is %15.12f and value is %e" % (root,surface_temp([root], mr_co2, Sflux,cc_flag)))



