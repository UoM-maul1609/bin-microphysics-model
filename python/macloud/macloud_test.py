import numpy as np

import numpy as np
from scipy.stats import norm
from scipy.optimize import root_scalar


def lognormal_number_between(D1, D2, Ntot, Dg, sigmag):
    """
    Number of particles between diameters D1 and D2 for a lognormal distribution.

    Parameters
    ----------
    D1, D2 : float
        Lower/upper diameter bounds (same units as Dg). Use np.inf for open upper bound.
    Ntot : float
        Total number concentration / total number (same units as desired output).
    Dg : float
        Geometric mean (median) diameter.
    sigmag : float
        Geometric standard deviation (> 1).

    Returns
    -------
    float
        Number between D1 and D2.
    """
    if D1 <= 0 or D2 <= 0 or Dg <= 0:
        raise ValueError("Diameters must be positive.")
    if not (sigmag > 1):
        raise ValueError("sigmag must be > 1.")

    s = np.log(sigmag)

    z2 = (np.log(D2) - np.log(Dg)) / s if np.isfinite(D2) else np.inf
    z1 = (np.log(D1) - np.log(Dg)) / s if np.isfinite(D1) else -np.inf

    frac = norm.cdf(z2) - norm.cdf(z1)
    return Ntot * frac
    
    
def calculate_temp(rh_guess):
	"""
		basically guess the RH
		calculate how much water on the aerosol assuming equilibrium RH
		calculate how much water needs to have evaporated and calculate the RH
	"""
	

"""
	store total water, calculate the temperature
"""
rgas=8.314
mair=29e-3
ra=rgas/mair
tinit=290.15
pinit=1e5
rhinit=0.65
rhoa_init=pinit/tinit/ra
lv=2.5e6
cp=1005.

"""
	aerosol properties
"""
NaClMR=1e-9
N_aer1=np.array([7246000,3132000,49800])*1e6
N_aer1=np.array([7246000,3132000,0])*1e6
N_aer1=N_aer1/np.sum(N_aer1)
logSig=np.array([0.63345,0.3607466,0.88649])
Dm=np.array([0.0395,0.1162,0.4477])*1e-6
rho_a=2165.0




# see moments of lognormal distribution
a=NaClMR/np.sum((N_aer1*np.pi/6.*2165.*np.exp(3.*np.log(Dm)+4.5*logSig**2)))
N_aer=a*np.ones([len(N_aer1),1])*(N_aer1*np.ones([1,1])).transpose()

# create bins from 5nm to 1micron, edges
Dpe=np.logspace(np.log10(5e-9),np.log10(1e-6),1001)
Dp=np.array([0.5*(Dpe[i]+Dpe[i+1]) for i in range(1000)])

N_between1 = [lognormal_number_between(
    D1=Dpe[i], D2=Dpe[i+1],
    Ntot=N_aer[0],        # e.g. cm^-3
    Dg=Dm[0],
    sigmag=np.exp(logSig[0])
) for i in range(1000)]
N_between2 = [lognormal_number_between(
    D1=Dpe[i], D2=Dpe[i+1],
    Ntot=N_aer[1],        # e.g. cm^-3
    Dg=Dm[1],
    sigmag=np.exp(logSig[1])
) for i in range(1000)]
N_between3 = [lognormal_number_between(
    D1=Dpe[i], D2=Dpe[i+1],
    Ntot=N_aer[2],        # e.g. cm^-3
    Dg=Dm[2],
    sigmag=np.exp(logSig[2])
) for i in range(1000)]
N_between=np.array(N_between1)+np.array(N_between2)+np.array(N_between3)

mp=np.pi/6*rho_a*Dp**3

msw=100/3.5*mp
S=3.5/100.
rho_sw=1025.
Dsw=Dp*(rho_a/(rho_sw*S))**(1/3)

total_mr_sw=np.sum(N_between*np.pi/6*rho_sw*Dsw**3)
total_mr_h2o=total_mr_sw*(1-S)

deltaT=total_mr_h2o*lv/(cp*rhoa_init)




