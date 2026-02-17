"""
	MACLOUD: Effect of water in the spray on the sub-cloud region
	P. Connolly
	Assumptions:
	1. Base case, conserve q, and theta_dry, calculate profile until saturation - done
	2. Water evaporates instantaneously and finds equilibrium - do full method first, so bin micro
	3. Then aerosols stay in equilibrium during ascent until saturation
"""

import numpy as np
from scipy.stats import norm
from scipy.optimize import root_scalar
from scipy.optimize import minimize_scalar
from scipy.integrate import odeint
import svp
import matplotlib.pyplot as plt
import sys

import numpy as np
from scipy.stats import norm

argdict=dict()
"""
	Constants and initial conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
rgas=8.314
mair=29e-3
mh2o=18e-3
ra=rgas/mair
rv=rgas/mh2o
tinit=290.15
pinit=1e5
rhinit=0.65
rhoa_init=pinit/tinit/ra
lv=2.5e6
cp=1005.
cpv=1870.
cpl=4.27e3
dz=10.;
grav=9.81
rhow=1000.
winit=0.3
eps1=mh2o/mair
N_flux=2e15 # spray particles per second
area_SEP=1.42e7*1e6 # m^2
depth_SEP=1000.0 # m
aerosol_lifetime=2*86400.0
rho_a=2165.0
rh_thresh=0.99
COOPER=1
HARRISON_EFFERVESCENCE=2
HARRISON_DELAVAL=3
EDMUND_SCF=4
RAYLEIGH=5
tau_cond=0.1
tau_cond2=10.0
nbins=60
long_integration=True
rlinit=0.0

def lognormal_volume_between(D1, D2, Ntot, Dg, sigmag):
    """
    Total particle VOLUME between diameters D1 and D2 for a single lognormal mode.

    Parameters
    ----------
    D1, D2 : float
        Lower/upper diameter bounds [m] (same units as Dg). Use np.inf for open upper bound.
    Ntot : float
        Total number concentration [#/m^3] for this mode.
    Dg : float
        Geometric mean (median) diameter [m].
    sigmag : float
        Geometric standard deviation (> 1), dimensionless.

    Returns
    -------
    float
        Total particle volume concentration between D1 and D2 [m^3 of particles per m^3 of air] (dimensionless).
    """
    if D1 <= 0 or D2 <= 0 or Dg <= 0:
        raise ValueError("Diameters must be positive.")
    if not (sigmag > 1):
        raise ValueError("sigmag must be > 1.")

    s = np.log(sigmag)

    # For a lognormal in D, the k-th moment truncated to [D1,D2] is:
    # E[D^k; D in [D1,D2]] = exp(k*mu + 0.5*k^2*s^2) * [Phi(z2-k*s) - Phi(z1-k*s)]
    mu = np.log(Dg)
    k = 3.0  # because volume ~ D^3

    z1 = (np.log(D1) - mu) / s if np.isfinite(D1) else -np.inf
    z2 = (np.log(D2) - mu) / s if np.isfinite(D2) else np.inf

    trunc_factor = norm.cdf(z2 - k * s) - norm.cdf(z1 - k * s)

    # Total volume = N * (pi/6) * E[D^3; truncated]
    vol = Ntot * (np.pi / 6.0) * np.exp(k * mu + 0.5 * (k ** 2) * (s ** 2)) * trunc_factor
    return vol


def sum_lognormal_volume_between(D1, D2, N, Dg, sigmag):
    """
    Total particle VOLUME between diameters D1 and D2 for a SUM of lognormal modes.

    Parameters
    ----------
    D1, D2 : float
        Lower/upper diameter bounds [m].
    N : array-like
        Total number concentration per mode [#/m^3].
    Dg : array-like
        Geometric mean (median) diameter per mode [m].
    sigmag : array-like
        Geometric standard deviation per mode (dimensionless, >1).

    Returns
    -------
    float
        Total particle volume concentration between D1 and D2 [m^3/m^3].
    """
    N = np.asarray(N, dtype=float)
    Dg = np.asarray(Dg, dtype=float)
    sigmag = np.asarray(sigmag, dtype=float)

    if not (N.shape == Dg.shape == sigmag.shape):
        raise ValueError("N, Dg, sigmag must have the same shape.")

    total = 0.0
    for i in range(N.size):
        total += lognormal_volume_between(D1, D2, N[i], Dg[i], sigmag[i])
    return total


def iterate_qv_bulk(qv,argdict):

	cp=argdict['cp']
	cpv=argdict['cpv']
	cpl=argdict['cpl']
	ra=argdict['ra']
	rv=argdict['rv']

	lv=argdict['lv']
	
	wv=qv
	wl=argdict['rtot']-qv # set
	
	rm=ra+rv*wv
	cpm=cp+wv*cpv+wl*cpl
	
	# now, calculate the temperature if there is more vapour now, the temp will be lower
	told=argdict['t']
	t=argdict['t']-(qv-argdict['wv'])*lv/cpm
	
	# calculate the RH
	es=svp.svp([t],'buck2','liq')[0]
	smr=argdict['eps1']*es/(argdict['p']-es)
	argdict['rh']=qv/smr
	
# 	GF1=(1.+argdict['rh']/(1.-argdict['rh'])* \
# 		argdict['mh2o']/argdict['molw_s_back']* \
# 		argdict['rhos_back']/argdict['rhow'] ) **(1./3.)
# 
# 	GF2=(1.+argdict['rh']/(1.-argdict['rh'])* \
# 		argdict['mh2o']/argdict['molw_s_spray']* \
# 		argdict['rhos_spray']/argdict['rhow'] ) **(1./3.)
	
		
	(q_back,q_spray)=calc_masses_bulk2(t)

	
	qtot=q_back+q_spray+qv  #+argdict['wl2']
	argdict['t']=told
	
	

	return qtot-argdict['rtot']

def calc_masses_bulk1(t):
	# lognormal 3rd moment of the dry aerosol
	qa1=argdict['rhos_back']*np.pi/6.0*np.sum(argdict['N_back']* \
		np.exp(3.*np.log(argdict['dm_back'])+4.5*argdict['lnsig_back']**2))

	qa2=argdict['rhos_spray']*np.pi/6.0*np.sum(argdict['N_spray']* \
		np.exp(3.*np.log(argdict['dm_spray'])+4.5*argdict['lnsig_spray']**2))

	# just the water, from raout term
	q1=argdict['rh']/(1.-argdict['rh'])*argdict['nu_back']* \
		qa1/argdict['molw_s_back']*argdict['mh2o']
	
	argdict['mass_s']=qa1
	argdict['molw_s']=argdict['molw_s_back']
	argdict['nu']=argdict['nu_back']
	argdict['rhos']=argdict['rhos_back']
	argdict['t']=t
	argdict['rh_mult']=1.0
	argdict['N_this']=argdict['N_back']
	# here, you need to solve a "bulk" koehler curve, changing the q1 until the RH is 
	# calculated properly
	# you can calc the average wet diameter for the kelvin term?
	sol2=root_scalar(koehler_curve_bulk,bracket=[argdict['mass_s']/1000,q1],\
			method='bisect',rtol=1e-10,xtol=1e-30,args=(argdict))
	
	q1=sol2.root
	# salinity S for the spray
	
	q2=qa2/argdict['S']-qa2
	return (q1,q2)

def calc_masses_bulk2(t):
	# lognormal 3rd moment of the dry aerosol
	qa1=argdict['rhos_back']*np.pi/6.0*np.sum(argdict['N_back']* \
		np.exp(3.*np.log(argdict['dm_back'])+4.5*argdict['lnsig_back']**2))

	qa2=argdict['rhos_spray']*np.pi/6.0*np.sum(argdict['N_spray']* \
		np.exp(3.*np.log(argdict['dm_spray'])+4.5*argdict['lnsig_spray']**2))

	# just the water, from raout term
	q1=argdict['rh']/(1.-argdict['rh'])*argdict['nu_back']* \
		qa1/argdict['molw_s_back']*argdict['mh2o']
		
	q2=argdict['rh']/(1.-argdict['rh'])*argdict['nu_spray']* \
		qa2/argdict['molw_s_spray']*argdict['mh2o']

	argdict['mass_s']=qa1
	argdict['molw_s']=argdict['molw_s_back']
	argdict['nu']=argdict['nu_back']
	argdict['rhos']=argdict['rhos_back']
	argdict['t']=t
	argdict['rh_mult']=1.0
	argdict['N_this']=argdict['N_back']
	# here, you need to solve a "bulk" koehler curve, changing the q1 until the RH is 
	# calculated properly
	# you can calc the average wet diameter for the kelvin term?
	sol2=root_scalar(koehler_curve_bulk,bracket=[argdict['mass_s']/1000,q1],\
			method='bisect',rtol=1e-10,xtol=1e-30,args=(argdict))
	q1=sol2.root

	argdict['mass_s']=qa2
	argdict['molw_s']=argdict['molw_s_spray']
	argdict['nu']=argdict['nu_spray']
	argdict['rhos']=argdict['rhos_spray']
	argdict['N_this']=argdict['N_spray']
	# here, you need to solve a "bulk" koehler curve, changing the q1 until the RH is 
	# calculated properly
	# you can calc the average wet diameter for the kelvin term?
	sol2=root_scalar(koehler_curve_bulk,bracket=[argdict['mass_s']/1000,q2],\
			method='bisect',rtol=1e-10,xtol=1e-30,args=(argdict))
	
	q2=sol2.root

	return (q1,q2)

def koehler_curve_bulk(x,argdict):
	# x is the water mass for this bin
	
	x=np.maximum(x,0)
	
	nw=x/argdict['mh2o']
	ms=argdict['mass_s']
	ns=ms/argdict['molw_s']
	nu=argdict['nu']
	rhos=argdict['rhos']
	rhow=argdict['rhow']
	t=argdict['t']
	rgas=argdict['rgas']
	mh2o=argdict['mh2o']
	
	rhoat=x/rhow+ms/rhos # volume or material per kg of air
	rhoat=(x+ms)/rhoat # density kg/m^3
	
	dw=((x+ms)/np.sum(argdict['N_this'])*6.0/(np.pi*rhoat))**(1.0/3.0)
	sigma=surface_tension(t)
	
	rh=np.exp(4.0*mh2o*sigma/rgas/t/rhoat/dw)*nw/(nw+nu*ns)*argdict['rh_mult']
	
	return rh-argdict['rh']

	
	
def iterate_qv(qv,argdict):
	# use this temperature to calculate the equilibrium
	cp=argdict['cp']
	cpv=argdict['cpv']
	cpl=argdict['cpl']
	ra=argdict['ra']
	rv=argdict['rv']

	lv=argdict['lv']
	
	wv=qv
	wl=argdict['rtot']-qv # set
	
	rm=ra+rv*wv
	cpm=cp+wv*cpv+wl*cpl
	
	# now, calculate the temperature if there is more vapour now, the temp will be lower
	told=argdict['t']
	t=argdict['t']-(qv-argdict['wv'])*lv/cpm
	
	# calculate the RH
	es=svp.svp([t],'buck2','liq')[0]
	smr=argdict['eps1']*es/(argdict['p']-es)
	argdict['rh']=qv/smr
	

	(mwat_back,mwat_spray)=calc_masses(t)

	qtot=(np.sum(mwat_back*argdict['num_bins_back'])+ \
		np.sum(mwat_spray*argdict['num_bins_spray']))+qv
	
	#print(qtot,argdict['rtot'],t,told,qtot-argdict['rtot'])
	argdict['t']=told

	return qtot-argdict['rtot']


def calc_masses(t):
	argdict['nu']=argdict['nu_back']
	argdict['rhos']=argdict['rhos_back']
	argdict['molw_s']=argdict['molw_s_back']
	# t=t
	argdict['t']=t
	argdict['rh_mult']=1.0
	mwat_back=np.zeros(argdict['nbins'])
	for i in range(argdict['nbins']):
		# argdict['rh']=set
		argdict['mass_s']=argdict['maer_back'][i]
		sol2=root_scalar(koehler_curve,bracket=[argdict['mass_s']/1000,argdict['mmax']],method='bisect',\
			rtol=1e-10,xtol=1e-30,args=(argdict))
		mwat_back[i]=sol2.root			

	argdict['nu']=argdict['nu_spray']
	argdict['rhos']=argdict['rhos_spray']
	argdict['molw_s']=argdict['molw_s_spray']
	mwat_spray=np.zeros(argdict['nbins'])
	for i in range(argdict['nbins']):
		#argdict['rh']=set
		argdict['mass_s']=argdict['maer_spray'][i]
		sol2=root_scalar(koehler_curve,bracket=[argdict['mass_s']/1000,argdict['mmax']],method='bisect',\
			rtol=1e-10,xtol=1e-30,args=(argdict))
		mwat_spray[i]=sol2.root		
		
	return (mwat_back,mwat_spray)

def surface_tension(t):
	tc=t-argdict['ttr']
	tc=np.maximum(tc,-40.0)
	surface_tension1 = 75.93 + 0.115 * tc + 6.818e-2 * tc**2 + 6.511e-3 * tc**3 + 2.933e-4 * tc**4 +  6.283e-6 * tc**5 + 5.285e-8 * tc**6
	
	if(tc > 0.0):
		surface_tension1 = 76.1 - 0.155*tc
		
	surface_tension1 = surface_tension1*1.e-7 # convert to j/cm2 
	surface_tension1 = surface_tension1*1.e4 # convert to j/m2 
	
	return surface_tension1

def koehler_curve(x,argdict):
	# x is the water mass for this bin
	
	x=np.maximum(x,0)
	
	nw=x/argdict['mh2o']
	ms=argdict['mass_s']
	ns=ms/argdict['molw_s']
	nu=argdict['nu']
	rhos=argdict['rhos']
	rhow=argdict['rhow']
	t=argdict['t']
	rgas=argdict['rgas']
	mh2o=argdict['mh2o']
	
	rhoat=x/rhow+ms/rhos # volume
	rhoat=(x+ms)/rhoat # density
	
	dw=((x+ms)*6.0/(np.pi*rhoat))**(1.0/3.0)
	
	sigma=surface_tension(t)
	
	rh=np.exp(4.0*mh2o*sigma/rgas/t/rhoat/dw)*nw/(nw+nu*ns)*argdict['rh_mult']
	
	return rh-argdict['rh']

def f_min_num(x,argdict):
	num1=lognormal_between_tot(argdict['minD'],x,\
		argdict['N_this'],argdict['dm_this'],np.exp(argdict['lnsig_this']))
	
	return num1-argdict['num_per_bin_this']

def lognormal_between_tot(D1,D2, N,Dg, sigmag):
	Ntot=0.0
	for i in range(len(N)):
		Ntot=Ntot+lognormal_number_between(D1,D2,N[i],Dg[i],sigmag[i])

	return Ntot

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
    
    	

def parcelwarm_simple(y,t,argdict):
	"""
		y=z,t,p,rv,rl
	"""
	cp=argdict['cp']
	cpv=argdict['cpv']
	cpl=argdict['cpl']
	ra=argdict['ra']
	rv=argdict['rv']
	eps1=argdict['eps1']
	tau_cond=argdict['tau_cond']
	lv=argdict['lv']
	
	wv=y[3]
	wl=y[4]+argdict['wl']
	
	rm=ra+rv*wv
	cpm=cp+wv*cpv+wl*cpl
	
	dydt=np.zeros_like(y)
	
	dydt[0]=argdict['winit']
	"""
		rate of change of pressure dp/dz=-rho*g
	"""
	dydt[2]=-y[2]/(rm*y[1])*argdict['grav']*argdict['winit']
	
	"""
		rate of change of temperature
	"""
	dydt[1]=rm/y[2]*dydt[2]*y[1]/cpm
	
	"""
		RH
	"""
	es=svp.svp([y[1]],'buck2','liq')[0]
	smr=eps1*es/(y[2]-es)*argdict['rh_thresh']
	#rh=y[3]/smr
	
	
	
	dydt[3]=1.0/tau_cond * (np.minimum(y[4],np.maximum(0,smr-y[3])) - \
		np.maximum(0.0,y[3]-smr))
	dydt[4]=-dydt[3]
	
	"""
		add latent heat
	"""
	dydt[1]=dydt[1]+lv*dydt[4]/cpm
	
	return dydt
	
def parcelwarm_equilibration(y,t,argdict):
	"""
		y=z,t,p,rv
	"""
	cp=argdict['cp']
	cpv=argdict['cpv']
	cpl=argdict['cpl']
	ra=argdict['ra']
	rv=argdict['rv']
	eps1=argdict['eps1']
	tau_cond=argdict['tau_cond']
	tau_cond2=argdict['tau_cond2']
	lv=argdict['lv']
	nbins=argdict['nbins']
	
	wv=y[3]
	num_bins=np.array(list(argdict['num_bins_back'])+list(argdict['num_bins_spray']))
	mbins=np.array(list(argdict['mwat_back'])+list(argdict['mwat_spray']))
	wl=np.sum(num_bins*mbins)+y[4]
	
	rm=ra+rv*wv
	cpm=cp+wv*cpv+wl*cpl
	
	dydt=np.zeros_like(y)
	
	dydt[0]=argdict['winit']
	"""
		rate of change of pressure dp/dz=-rho*g
	"""
	dydt[2]=-y[2]/(rm*y[1])*argdict['grav']*argdict['winit']
	
	"""
		rate of change of temperature
	"""
	dydt[1]=rm/y[2]*dydt[2]*y[1]/cpm
	
	"""
		RH
	"""
	es=svp.svp([y[1]],'buck2','liq')[0]
	smr=eps1*es/(y[2]-es)*argdict['rh_thresh']
	#rh=y[3]/smr
	
	
	"""
		calcaulate the equilibrium mass for each aerosol at this humidity ++++++++++++++++
		
	"""
	dydt[3]=1.0/tau_cond * (np.minimum(y[4],np.maximum(0,smr-y[3])) - \
		np.maximum(0.0,y[3]-smr)) 
	dydt[4]=-dydt[3]
	"""
		----------------------------------------------------------------------------------
	"""
	
	"""
		add latent heat
	"""
	dydt[1]=dydt[1]+lv*(dydt[4])/cpm
	
	return dydt

def do_it(dx_grid=100.,mass_flag=False,spray_method=1,M_flux1=50e9):

	M_flux=M_flux1/(365.25*86400) # 50 Tg per year in kg/s
	z=np.mgrid[0:2000+dz:dz]
	t=np.zeros_like(z)
	p=np.zeros_like(z)
	rv1=np.zeros_like(z)
	rl1=np.zeros_like(z)
	maer=np.zeros((len(z),nbins))
	
	esinit=svp.svp([tinit],'buck2','liq')[0]
	rvinit=eps1*esinit/(pinit-esinit)*rhinit
	#dx_grid=100.0
	dy_grid=dx_grid
	dz_grid=dz
	area_grid=dx_grid*dy_grid
	vol_grid=dx_grid*dy_grid*dz_grid
	#mass_flag=False # if false do it by number flux
	
	
	"""
		--------------------------------------------------------------------------------------
	"""
	
	"""
		aerosol parameters +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	"""
# 	spray_method=HARRISON_EFFERVESCENCE
# 	spray_method=5
	
	N_back		=np.array([46.642, 153.421, 166.774])
	lnsig_back	=np.array([0.348, 0.354, 0.465])
	dm_back		=np.array([0.018, 0.039, 0.154])
	
	if spray_method==COOPER:
		N_spray		=np.array([1.55,194,17.7])
		lnsig_spray	=np.array([0.129,0.625,0.666])
		dm_spray	=np.array([0.0263,0.0588,0.269])	
	elif spray_method==HARRISON_EFFERVESCENCE:
		N_spray		=np.array([94800,245000,137000])
		lnsig_spray	=np.array([0.252,0.497,0.883])
		dm_spray	=np.array([0.0264,0.0434,0.0574])	
	elif spray_method==HARRISON_DELAVAL:
		N_spray		=np.array([214000,66800,4800])
		lnsig_spray	=np.array([0.703,0.534,0.9])
		dm_spray	=np.array([0.0156,0.123,0.398])	
	elif spray_method==EDMUND_SCF:
		N_spray		=np.array([6460000,3830000,32700])
		lnsig_spray	=np.array([0.571,0.391,0.739])
		dm_spray	=np.array([0.0366,0.109,0.696])	
	elif spray_method==RAYLEIGH:
		N_spray		=np.array([1.0])
		lnsig_spray	=np.array([0.25])
		dm_spray	=np.array([0.162])	
	
	minD=5e-9 
	maxD=10e-6
	maxD1=1.0
	#print(lognormal_between_tot(1e-9,1e-6,N_spray*1e6,dm_spray*1e-6,np.exp(lnsig_spray)))
	
	num_tot_back=lognormal_between_tot(minD,maxD,N_back*1e6,dm_back*1e-6,np.exp(lnsig_back))
	num_per_bin_back=num_tot_back/nbins
	dedges_back=np.zeros(nbins+1)
	dedges_back[0]=minD
	num_tot_spray=lognormal_between_tot(minD,maxD,N_spray*1e6,dm_spray*1e-6,np.exp(lnsig_spray))
	num_per_bin_spray=num_tot_spray/nbins
	dedges_spray=np.zeros(nbins+1)
	dedges_spray[0]=minD
	argdict['N_back']=N_back*1e6
	argdict['dm_back']=dm_back*1e-6
	argdict['lnsig_back']=lnsig_back
	argdict['num_per_bin_back']=num_per_bin_back
	
	argdict['N_spray']=N_spray*1e6
	argdict['dm_spray']=dm_spray*1e-6
	argdict['lnsig_spray']=lnsig_spray
	argdict['num_per_bin_spray']=num_per_bin_spray
	
	argdict['N_this']=N_back*1e6
	argdict['dm_this']=dm_back*1e-6
	argdict['lnsig_this']=lnsig_back
	argdict['num_per_bin_this']=num_per_bin_back
	vol_back=np.zeros(nbins)
	vol_spray=np.zeros(nbins)
	
	for i in range(nbins):
		argdict['minD']=dedges_back[i]
		sol1=root_scalar(f_min_num,bracket=[dedges_back[i],maxD1],method='brentq',\
			rtol=1e-10,args=(argdict))
		dedges_back[i+1]=sol1.root
		
		vol_back[i]=sum_lognormal_volume_between(\
			dedges_back[i], dedges_back[i+1], argdict['N_back'], \
			argdict['dm_back'], np.exp(argdict['lnsig_back']))/num_per_bin_back
		
		
	# 	argdict['num_per_bin_this']=0.0
	# 	print(f_min_num(dedges_back[i+1],argdict)/num_per_bin_back)
	# 	argdict['num_per_bin_this']=num_per_bin_back
	
	argdict['N_this']=N_spray*1e6
	argdict['dm_this']=dm_spray*1e-6
	argdict['lnsig_this']=lnsig_spray
	argdict['num_per_bin_this']=num_per_bin_spray
	for i in range(nbins):
		argdict['minD']=dedges_spray[i]
		try:
			sol1=root_scalar(f_min_num,bracket=[dedges_spray[i],maxD1],method='brentq',\
				rtol=1e-10,args=(argdict))
			dedges_spray[i+1]=sol1.root
		except:
			dedges_spray[i+1]=dedges_spray[i]*100.0
			
	
		vol_spray[i]=sum_lognormal_volume_between(\
			dedges_spray[i], dedges_spray[i+1], argdict['N_spray'], \
			argdict['dm_spray'], np.exp(argdict['lnsig_spray']))/num_per_bin_spray
	
	# 	argdict['num_per_bin_this']=0.0
	# 	print(f_min_num(dedges_spray[i+1],argdict)/num_per_bin_spray)
	# 	argdict['num_per_bin_this']=num_per_bin_spray
			
		
	
	
	# bin mid-points
	dmid_back=np.array([0.5*(dedges_back[i]+dedges_back[i+1]) for i in range(nbins)])
	dmid_spray=np.array([0.5*(dedges_spray[i]+dedges_spray[i+1]) for i in range(nbins)])
	
	dmid_back=(vol_back*6/np.pi)**(1/3)
	dmid_spray=(vol_spray*6/np.pi)**(1/3)
	
	argdict['num_bins_back']=argdict['num_per_bin_back']*np.ones(nbins)
	argdict['num_bins_spray']=argdict['num_per_bin_spray']*np.ones(nbins)
	
	# for the background aerosol, equilibriate at this RH, T, P, etc
	argdict['winit']=winit
	argdict['ra']=ra
	argdict['rv']=rv
	argdict['grav']=grav
	argdict['cp']=cp
	argdict['cpv']=cpv
	argdict['cpl']=cpl
	argdict['eps1']=eps1
	argdict['tau_cond']=tau_cond
	argdict['tau_cond2']=tau_cond2
	argdict['lv']=lv
	argdict['rh']=rhinit
	argdict['t']=tinit
	argdict['p']=pinit
	argdict['ttr']=273.15
	argdict['rgas']=rgas
	argdict['mh2o']=mh2o
	argdict['rhow']=rhow
	argdict['rhos_back']=1770.0
	argdict['nu_back']=2.5
	argdict['molw_s_back']=132.0e-3
	argdict['rhos_spray']=2165.0
	argdict['nu_spray']=2.0
	argdict['molw_s_spray']=58.0e-3
	argdict['S']=35e-3 # kg salt per kg of sea water
	argdict['rh_mult']=1.0
	argdict['rh_thresh']=rh_thresh
	
	"""
		nw=x/argdict['mh2o']
		ms=argdict['mass_s']
		ns=ms/argdict['molw_s']
		nu=argdict['nu']
		rhos=argdict['rhos']
		rhow=argdict['rhow']
		t=argdict['t']
		rgas=argdict['rgas']
	
	"""
	
	"""
		background aerosol
	"""
	mwat_back=np.zeros(nbins)
	maer_back=np.zeros(nbins)
	rh_back=np.zeros(nbins)
	
	mmax=10000e-6**3*np.pi/6.0*argdict['rhow']
	argdict['nu']=argdict['nu_back']
	argdict['rhos']=argdict['rhos_back']
	argdict['molw_s']=argdict['molw_s_back']
	argdict['t']=tinit
	for i in range(nbins):
		argdict['rh']=rhinit
		argdict['mass_s']=dmid_back[i]**3*np.pi/6.0*argdict['rhos']
		sol2=root_scalar(koehler_curve,bracket=[argdict['mass_s']/1000.,mmax],method='bisect',\
			rtol=1e-10,xtol=1e-30,args=(argdict))
		mwat_back[i]=sol2.root
		maer_back[i]=argdict['mass_s']
		
		argdict['rh']=0.0
		rh_back[i]=koehler_curve(mwat_back[i],argdict)
	
	"""
		spray aerosol
	"""
	mwat_spray=np.zeros(nbins)
	mwat_spray2=np.zeros(nbins)
	maer_spray=np.zeros(nbins)
	rh_spray=np.zeros(nbins)
	
	argdict['mmax']=mmax
	argdict['nu']=argdict['nu_spray']
	argdict['rhos']=argdict['rhos_spray']
	argdict['molw_s']=argdict['molw_s_spray']
	argdict['t']=tinit
	for i in range(nbins):
		argdict['rh']=rhinit
		argdict['mass_s']=dmid_spray[i]**3*np.pi/6.0*argdict['rhos']
		sol2=root_scalar(koehler_curve,bracket=[argdict['mass_s']/1000.,mmax],method='bisect',\
			rtol=1e-10,xtol=1e-30,args=(argdict))
		mwat_spray[i]=sol2.root
		maer_spray[i]=argdict['mass_s']
		
		argdict['rh']=0.0
		rh_spray[i]=koehler_curve(mwat_spray[i],argdict)
		
		mwat_spray2[i]=maer_spray[i]/argdict['S']-maer_spray[i]
	
	"""
		RHcrit and maer_crit of spray aerosol
	"""
	mwat_spray3=np.zeros(nbins)
	rh_spray3=np.zeros(nbins)
	
	argdict['rh']=0.0
	argdict['rh_mult']=-1.0
	for i in range(nbins):
		argdict['mass_s']=dmid_spray[i]**3*np.pi/6.0*argdict['rhos']
		sol2=minimize_scalar(koehler_curve,bracket=[mwat_spray[i],mmax],method='golden',\
			tol=1e-30,args=(argdict))
		mwat_spray3[i]=sol2.x
		rh_spray3[i]=-sol2.fun
			
	argdict['rh_mult']=1.0
	
	
	
	"""
		--------------------------------------------------------------------------------------
	"""
	
	
	
	"""
		calculate the mass in the spray system++++++++++++++++++++++++++++++++++++++++++++++++
	"""
	# kg/m^3
	mass_spray=np.sum((N_spray*1e6*np.pi/6.*rho_a*np.exp(3.*np.log(dm_spray*1e-6)+4.5*lnsig_spray**2)))
	mass_content_from_mass_flux=M_flux/(area_SEP*depth_SEP)*aerosol_lifetime # kg/m^3
	# ratio by which we need to scale the number concs to match this bass flux
	ratio_of_mass_contents=mass_content_from_mass_flux / mass_spray 
	N_conc_for_this_mass_flux=np.sum(N_spray*ratio_of_mass_contents) # /cc
	
	# now what number flux would this be?over the whole SEP?
	N_flux_calc=N_conc_for_this_mass_flux*1e6*area_SEP*depth_SEP/(aerosol_lifetime)
	
	# next, calculate the number of sulphate aerosol in each bin, and number of nacl in each bin
	# adjust the water mass in each bin
	
	if mass_flag:
		argdict['num_bins_spray']=argdict['num_bins_spray']*N_conc_for_this_mass_flux/ \
			np.sum(N_spray)
		N_spray=N_spray*N_conc_for_this_mass_flux/ \
			np.sum(N_spray)
		argdict['N_spray']=N_spray*1e6
	else:
		# Nflux/area/winit / N_spray *
		argdict['num_bins_spray']=argdict['num_bins_spray']*N_flux/(winit*area_grid*1e6)/ \
			np.sum(N_spray)
		N_spray=N_spray*N_flux/(winit*area_grid*1e6)/ \
			np.sum(N_spray)
		argdict['N_spray']=N_spray*1e6
	
		
	
	"""
		--------------------------------------------------------------------------------------
	"""
	
	
	
	
	"""
	 1. calculate mr. theta, etc, just use
	"""
	
	
	"""
	  1.a Simple ODE solver+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	"""
	y0=[z[0],tinit,pinit,rvinit,rlinit]
	t[0]=tinit
	p[0]=pinit
	rv1[0]=rvinit
	rl1[0]=rlinit
	t_init=0.0
	argdict['wl']=0.0
	for i in range(len(z)-1):
		tsolve=t_init+dz/winit
		sol = odeint(parcelwarm_simple, y0, [t_init,tsolve], args=(argdict,))
		
		t[i+1]=sol[1][1]
		p[i+1]=sol[1][2]
		rv1[i+1]=sol[1][3]
		rl1[i+1]=sol[1][4]
		y0=[z[i+1],t[i+1],p[i+1],rv1[i+1],rl1[i+1]]
		t_init=tsolve
	
	"""
	 1b. Calculate the RH
	"""
	es=np.array(svp.svp(t,'buck2','liq'))
	smr=argdict['eps1']*es/(p-es)
	rh=rv1/smr
	"""
	   ---------------------------------------------------------------------------
	"""
	
	"""
	  2a. with aerosol equilibration
	"""
	if long_integration:
		t_a			=np.zeros_like(t)
		p_a			=np.zeros_like(p)
		rv1_a		=np.zeros_like(rv1)
		rl1_a		=np.zeros_like(rv1)
		rl2_a		=np.zeros_like(rv1)
		mwat_back_a	=np.zeros((len(t_a),nbins))
		mwat_spray_a=np.zeros((len(t_a),nbins))
		argdict['nbins']=nbins
		argdict['maer_back']=maer_back
		argdict['maer_spray']=maer_spray
		argdict['mwat_back']=mwat_back
		argdict['mwat_spray']=mwat_spray2
		rl2_a[0]=np.sum(mwat_back*argdict['num_bins_back'])+ \
			np.sum(mwat_spray2*argdict['num_bins_spray'])
		t_a[0]=tinit
		p_a[0]=pinit
		rv1_a[0]=rvinit
		rl1_a[0]=rlinit
		RH=np.zeros(nbins)
		#argdict['num_bins_spray'][:]=0.0
		
		y0=[z[0],tinit,pinit,rvinit,rlinit]
		#y0=y0+list(mwat_back)+list(mwat_spray)
		#y0=y0+list(mwat_spray2)
		t_init=0.0
		for i in range(len(z)-1):
			print(t_init)
			tsolve=t_init+dz/winit
			sol = odeint(parcelwarm_equilibration, y0, [t_init,tsolve], args=(argdict,), \
				atol=[1e-2,1e-5,1e-1,1e-10,1e-10])
			t_a[i+1]=sol[1][1]
			p_a[i+1]=sol[1][2]
			rv1_a[i+1]=sol[1][3]
			rl1_a[i+1]=sol[1][4]
			
			# do equilibration, first calculate the RH
			es=svp.svp([t_a[i+1]],'buck2','liq')[0]
			smr=argdict['eps1']*es/(p_a[i+1]-es)
			rh1=rv1_a[i+1]/smr
			
			# only if RH is less than RH thresh do we calculate the sub-saturated water content
			if rh1<argdict['rh_thresh']:
				argdict['t']=t_a[i+1]
				argdict['mwat_back']=mwat_back # the old water content
				argdict['mwat_spray']=mwat_spray2 # the old water content
				argdict['wv']=rv1_a[i+1]
				argdict['wl']=rl1_a[i+1]
				argdict['wl2']=np.sum(argdict['mwat_back']*argdict['num_bins_back'])+ \
					np.sum(argdict['mwat_spray']*argdict['num_bins_spray'])
			
				# total water - this is what I will root find
				# i.e. vapour plus liquid has to match this
				argdict['rtot']=argdict['wv']+argdict['wl']+argdict['wl2']
				
				argdict['p']=p_a[i+1]
				sol3=root_scalar(iterate_qv,x0=rv1_a[i+1],method='newton',\
					rtol=1e-6,xtol=1e-10,args=(argdict))
				rm=ra+rv*sol3.root
				cpm=cp+sol3.root*cpv+(argdict['rtot']-sol3.root)*cpl
			
				argdict['t']=t_a[i+1]
				t1=argdict['t']-(sol3.root-argdict['wv'])*lv/cpm
				argdict['t']=t1
				
				t_a[i+1]=t1
				rv1_a[i+1]=sol3.root
				# now calc the new masses
				es=svp.svp([t_a[i+1]],'buck2','liq')[0]
				smr=argdict['eps1']*es/(p_a[i+1]-es)
				argdict['rh']=rv1_a[i+1]/smr
				(mwat_back,mwat_spray2)=calc_masses(t_a[i+1])
				
				# calculate equilibrium rh
	# 			argdict['t']=t_a[i+1]
	# 			argdict['rhos']=argdict['rhos_spray']
	# 			argdict['nu']=argdict['nu_spray']
	# 			argdict['molw_s']=argdict['molw_s_spray']
	# 			argdict['rh']=0.0
	# 			for j in range(nbins):
	# 				# x is the water mass for this bin
	# 				argdict['mass_s']=maer_spray[j]
	# 				RH[j]=koehler_curve(mwat_spray2[j],argdict)
	# 			
	# 			argdict['rh']=rv1_a[i+1]/smr
	# 			print(RH,argdict['rh'])
	
	# 			test2=rv1_a[i+1]+rl1_a[i+1]+np.sum(mwat_back*argdict['num_bins_back'])+ \
	# 				np.sum(mwat_spray2*argdict['num_bins_spray'])
	# 			print(argdict['rtot']/test2)
			
			rl2_a[i+1]=np.sum(mwat_back*argdict['num_bins_back'])+ \
				np.sum(mwat_spray2*argdict['num_bins_spray'])
				
			y0=[z[i+1],t_a[i+1],p_a[i+1],rv1_a[i+1],rl1_a[i+1]] #+ list(sol[1][5:])
			t_init=tsolve
		
		"""
		 2b. Calculate the RH
		"""
		es_a=np.array(svp.svp(t_a,'buck2','liq'))
		smr_a=argdict['eps1']*es_a/(p_a-es_a)
		rh_a=rv1_a/smr_a
	
	# sys.exit()
	"""
	  3. do the simple scheme based on moments
	"""
	t_b			=np.zeros_like(t)
	p_b			=np.zeros_like(p)
	rv1_b		=np.zeros_like(rv1)
	rl1_b		=np.zeros_like(rv1)
	rl2_b		=np.zeros_like(rv1)
	t_b[0]=tinit
	p_b[0]=pinit
	rv1_b[0]=rvinit
	rl1_b[0]=rlinit
	
	y0=[z[0],tinit,pinit,rvinit,rlinit]
	t_init=0.0
	
	# at this humidity, etc calculate the new diameter, etc
	argdict['rh']=rhinit
	(q_back,q_spray)=calc_masses_bulk1(t_b[0])
	rl2_b[0]=q_back+q_spray
	argdict['wl']=q_back+q_spray
	for i in range(len(z)-1):
		tsolve=t_init+dz/winit
		argdict['wl']=q_back+q_spray
		sol = odeint(parcelwarm_simple, y0, [t_init,tsolve], args=(argdict,))
		t_b[i+1]=sol[1][1]
		p_b[i+1]=sol[1][2]
		rv1_b[i+1]=sol[1][3]
		rl1_b[i+1]=sol[1][4]
		
		# do equilibration
		es=svp.svp([t_b[i+1]],'buck2','liq')[0]
		smr=argdict['eps1']*es/(p_b[i+1]-es)
		rh1=rv1_b[i+1]/smr
		
		
		if rh1<argdict['rh_thresh']:
			argdict['t']=t_b[i+1]
	
			argdict['wv']=rv1_b[i+1]
			argdict['wl']=rl1_b[i+1]+q_back+q_spray
		
			# total water - this is what I will root find
			argdict['rtot']=argdict['wv']+argdict['wl']
			argdict['wl2']=rl1_b[i+1] # this is the liquid water not on aerosol
			argdict['p']=p_b[i+1]
			sol3=root_scalar(iterate_qv_bulk,x0=rv1_b[i+1],method='newton',\
				rtol=1e-6,xtol=1e-10,args=(argdict))
			rm=ra+rv*sol3.root
			cpm=cp+sol3.root*cpv+(argdict['rtot']-sol3.root)*cpl
		
			t1=argdict['t']-(sol3.root-argdict['wv'])*lv/cpm
			argdict['t']=t1
			
			t_b[i+1]=t1
			rv1_b[i+1]=sol3.root
			# now calc the new masses
			es=svp.svp([t_b[i+1]],'buck2','liq')[0]
			smr=argdict['eps1']*es/(p_b[i+1]-es)
			argdict['rh']=rv1_b[i+1]/smr
			(q_back,q_spray)=calc_masses_bulk2(t_b[i+1])
			
		rl2_b[i+1]=q_back+q_spray
		
		y0=[z[i+1],t_b[i+1],p_b[i+1],rv1_b[i+1],rl1_b[i+1]] #+ list(sol[1][5:])
		t_init=tsolve
	
	"""
	 3b. Calculate the RH
	"""
	es_b=np.array(svp.svp(t_b,'buck2','liq'))
	smr_b=argdict['eps1']*es_b/(p_b-es_b)
	rh_b=rv1_b/smr_b
	
	
	outDict=dict()
	outDict['z']=z
	outDict['t']=t
	outDict['p']=p
	outDict['rh']=rh
	outDict['rv1']=rv1
	outDict['rl1']=rl1
	outDict['t_a']=t_a
	outDict['p_a']=p_a
	outDict['rh_a']=rh_a
	outDict['rv1_a']=rv1_a
	outDict['rl1_a']=rl1_a
	outDict['rl2_a']=rl2_a
	outDict['t_b']=t_b
	outDict['p_b']=p_b
	outDict['rh_b']=rh_b
	outDict['rv1_b']=rv1_b
	outDict['rl1_b']=rl1_b
	outDict['rl2_b']=rl2_b

	return outDict

# sys.exit()
if __name__=="__main__":
	outDict=do_it()
	



