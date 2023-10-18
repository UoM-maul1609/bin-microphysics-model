import numpy as np

def tau_vis(tau):
    ret1=np.zeros_like(tau)
    ind1,ind2=np.where(tau>=0.723)
    ret1[ind1,ind2]=0.36*(tau[ind1,ind2]-0.723)**0.411

    return ret1

def surface_temp(albedo_geo,frac_tot,co2ppm):
   Sflux = 1365.
   cloud_fraction=0.70
   albedo_clear=0.15
   albedo_c = 0.356
   frac_tot1=0.175*frac_tot
   albedo_p=(cloud_fraction-frac_tot1)*albedo_c+ \
    (1.-(cloud_fraction-frac_tot1))*albedo_clear+frac_tot1*albedo_geo;
    
   pres=0.8e5    # pressure to convert mixing ratios to partial pressures

   #mr_co2=400.  # mixing ratio of co2 (ppmv)
   #mr_co2=0.1e6
   mr_ch3=1.75  # mixing ratio of methane (ppmv)
   #mr_h2o=2700. # mixing ratio of water vapour (ppmv)

   mr_h2o=0.12*610.*np.exp(2.5e6/461.*(-1./289.+1./273.15))/0.8e5*1e6

   Ap=albedo_p        # planetary albedo
   sigma=5.67e-8  # stefan-boltzmann constant
   #Sflux=1368.    # solar flux (W m-2)
   epsil_s=0.95   # infrared emissivity of the surface

   e_co2=co2ppm/1.e6*pres # conversion from total pressure to partial pressure
   e_ch3=mr_ch3/1.e6*pres
   e_h2o=mr_h2o/1.e6*pres

   # Power law to convert to optical depth
   tau_co2=0.029*np.sqrt(e_co2)
   tau_ch3=0.725*np.sqrt(e_ch3)
   tau_h2o=0.097*np.sqrt(e_h2o)

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

   return T_s
   #print 'The surface temperature is %.2f K' % T_s



