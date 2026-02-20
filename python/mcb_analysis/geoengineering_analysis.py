import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import getpass
from scipy.interpolate import interp1d
import greenhouse_func_geo
import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, '../')
import lutMCB
from scipy.optimize import fsolve
import multi_root
import levenson_2024
import scipy.optimize as sco
from scipy.signal import savgol_filter, find_peaks

username=getpass.getuser()

# stern review
def scenario(temp_C, kind="central"):
    base = loss_pct(temp_C)
    if kind == "central":
        return base
    if kind == "upper":
        return loss_upper(temp_C)
    if kind == "lower":
        return loss_lower(temp_C)
    if kind == "ineq_125":
        return 1.25 * base
    if kind == "nonmarket_add12":
        return base + 12.0
    if kind == "amp_add6":
        return base + 6.0
    raise ValueError("unknown scenario")

def loss_lower(temp_C):
    temp_C = np.asarray(temp_C)
    out = np.zeros_like(temp_C, dtype=float)
    # below 2.5°C: ~0%
    # above: ramp with a power law anchored at 5.5°C -> 5%
    mask = temp_C > 2.5
    # one-point anchored curve: choose exponent ~3 like the central fit
    b_l = b
    a_l = 5.0 / (5.5**b_l)
    temp_C = np.maximum(temp_C,0.0)
    out[mask] = a_l * temp_C[mask]**b_l
    return out


def loss_upper(temp_C):
	T2 = np.array([2.5, 5.5])
	L_upper = np.array([3.0, 10.0])
	
	b_u = np.log(L_upper[1]/L_upper[0]) / np.log(T2[1]/T2[0])
	a_u = L_upper[0] / (T2[0]**b_u)
	temp_C = np.asarray(temp_C)
	temp_C = np.maximum(temp_C,0.0)
	return a_u * temp_C**b_u


def loss_pct(temp_C):
	# Midpoint “targets”
	T = np.array([2.5, 5.5])
	L = np.array([1.5, 7.5])  # percent GDP loss
	
	# Fit L = a*T^b using log regression
	b = np.log(L[1]/L[0]) / np.log(T[1]/T[0])
	a = L[0] / (T[0]**b)
	temp_C = np.asarray(temp_C)
	temp_C = np.maximum(temp_C,0.0)
	return a * temp_C**b



def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def maxAmaxMR(mr,A):
	
	diff2=np.flip(np.gradient(np.gradient(A)))
	diff1=np.flip(np.gradient(A))
	for i in range(len(diff1)-1):
		if ((diff1[i]<0) and (diff1[i+1]>=0) and (diff2[i]<0) ):
			break
	
	mr2=np.flip(mr)
	A2=np.flip(A)
	mr_peak=0.5*(mr2[i+1]+mr2[i+1])
	A_peak=0.5*(A2[i+1]+A2[i+1])
	
		
	
	return (mr_peak,A_peak)
	

""""++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% this function outputs the cloud albedo given the mixing ratio of NaCl and
% the diameter of the aerosol.
% input D in nm, mr is kg/kg
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
def cloud_albedo(mr,D,beta):
    a=beta[1]+beta[2]*D;
    b=beta[3]+beta[4]*D;
    yhat=(1.-beta[0])/(1.+np.exp(-a*(b+np.log10(mr))))+beta[0];
    
    return yhat

def solve1(x):
	return scint1(10**x)-albedo_geog
    
if __name__=="__main__":
    """
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % function to do the cost-benefit analysis of geoengineering
    %
    % 1. First of all it plots out the dependence of cloud albedo on the mass
    %       and size of NaCl particles added to a typical background marine
    %       stratocumulus cloud.
    % 2. We estimate the number of ships, the amount of sea water that needs to
    %       be pumped into the clouds, and use a climate model 'albedo_geo.m'
    %       to estimate the rise in temperature of different amounts of co2 and
    %       geoengineering. This is converted into a cost.
    % 3. We estimate the cost to society of this temperature rise. These
    %       equations are rough estimates from the Stern review.
    % 4. These costs are added and put as the z-axis on a plot
    %
    % costs of ships and maintenance are taken from various literature sources
    % that you must find.
    %--------------------------------------------------------------------------
    """
    lut_flag = True
    new_figs = False
    f=1.0;
    D_aerosol=45.;               # size of the aerosol particles generated (nm)
    
    if(lut_flag):
        (lut1,lut2,lut3,lut4)=lutMCB.doAnalysis()
        # calculate the average size, number weighted
        numer=np.sum(batchRunsMCB.N_aer1* \
        	np.exp(np.log(batchRunsMCB.Dm)+0.5*batchRunsMCB.logSig**2))
        D_aerosol=numer/np.sum(batchRunsMCB.N_aer1)*1e9
        #D_aerosol=lutMCB.Dm*1e9
    
    frac_tot1=f*17.5/100;     # fraction of earth geoengineered
    num_ship=f*1500;           # 1500 ships needed for 17.5% seeding
    time1=1;                   # implement for 1 year, 
    cost_energy=3.9e-8;         # cost of energy - dollars per joule 
    ship_lifetime=30 # 30 year lifetime
    cost_fac=4e9*(1./ship_lifetime+0.025) # 4 billion, but 2.5% of that per year maintenance
	
    RAYLEIGH_JET=1;
    TAYLOR_CONE=2; # not added this PSD yet
    SUPERCRITICAL=3;
    EFFERVESCENT=4;
    
    if batchRunsMCB.spray_method==batchRunsMCB.COOPER:
    	method=EFFERVESCENT;
    if batchRunsMCB.spray_method==batchRunsMCB.HARRISON_EFFERVESCENCE:
    	method=EFFERVESCENT;
    if batchRunsMCB.spray_method==batchRunsMCB.HARRISON_DELAVAL:
    	method=SUPERCRITICAL;
    if batchRunsMCB.spray_method==batchRunsMCB.EDMUND_SCF:
    	method=SUPERCRITICAL;
    if batchRunsMCB.spray_method==batchRunsMCB.RAYLEIGH:
    	method=RAYLEIGH_JET;
    if batchRunsMCB.spray_method==batchRunsMCB.TAYLOR_CONE:
    	method=TAYLOR_CONE;


    """
    "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % 1. Dependence of cloud albedo on mass and size of NaCl
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    # fit parameters
    beta=[0.356539140609705 ,  1.435442271779103 ,  \
        0.000831243139779 ,  7.202626231940690 , -0.002074725642290]; 

    mr=np.logspace(-14,-4,100)
    mr1=mr
    if new_figs:
    	fig1=plt.figure() #(name='albedo relationship');
    else:
    	plt.figure(fig1.number)
		
    D=100.e-9
    A=cloud_albedo(mr,D,beta);
    if (lut_flag):
        A=lut2
        mr1=lutMCB.NaClMR
        
    plt.plot(mr1,A);
    plt.xscale('log')
#     plt.ylabel('NaCl mode diameter (nm)');
    plt.xlabel('NaCl mixing ratio')
    plt.ylabel('albedo of geoengineered clouds')
    plt.savefig('/tmp/' + username + '/albedo_vs_NaCl_and_diameter.png')
    #--------------------------------------------------------------------------



    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % 2.a climate calculation 
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    frac_tot,co2ppm=np.mgrid[0:1+0.001:0.001,400:1510:10]
    frac_tot=np.logspace(-5,0,100)
    co2ppm=np.linspace(400,1500,111)
    frac_tot,co2ppm=np.meshgrid(frac_tot,co2ppm)
    frac_tot=frac_tot.transpose()
    co2ppm=co2ppm.transpose()
    # make the albedo 0.04 higher than baseline
    (mr_peak,A_peak)=maxAmaxMR(mr1,A)
#     albedo_geog=A[0]+0.02;
    albedo_geog=A_peak

    num_ship=frac_tot*1500;
    # in order to achieve albedo_geog we need to have the following mr
    if lut_flag:
    	#A(mr) and mr(A)
        scint=interp1d(A,mr1,fill_value='extrapolate')
        scint1=interp1d(mr1,A,fill_value='extrapolate')
        # find the mixing ratio that achieved albedo_geog
        mr=np.min(10**multi_root.multi_root( \
        	solve1,[np.log10(mr1[0]),np.log10(mr1[-1])]))
        mr=mr_peak
    else:
        scint=interp1d(A,mr,fill_value='extrapolate')
        scint1=interp1d(mr,A,fill_value='extrapolate')
        mr=np.min(10**multi_root.multi_root( \
        	solve1,[np.log10(mr1[0]),np.log10(mr1[-1])]))
    


    # now calculate the flow rate
    h=1000;                         # height of the boundary layer (meters)
    tau=3.*86400;                   # time scale aerosol fall out of atmosphere
                                    # 3 days
    Q=4.*np.pi*6.4e6**2.*frac_tot*17.5/100.*h*mr*1.2/(tau);  # equation B2
                                    # kg of NaCl per second

    Qsea=(Q)/(35);                 # m^3 of sea water per second


    #tsurface = greenhouse_func_geo.surface_temp(albedo_geog,frac_tot,co2ppm )
    
    (r,c)=np.shape(co2ppm)
    tsurface = np.zeros((r,c))  
    for i in range(r):
    	for j in range(c):  
    		levenson_2024.PCO2 = co2ppm[i][j]*levenson_2024.Ps*1e-6
    		levenson_2024.ngeou = frac_tot[i][j]
    		levenson_2024.Asg=A[0]
    		levenson_2024.Asgn=albedo_geog
    		roots=sco.fsolve(levenson_2024.do_calc, 289.)
    		tsurface[i][j]=roots[0]

    #tsurface=albedo_geo(albedo_geog,frac_tot,co2ppm);   # call climate model 
                                                        # for all co2s
    #--------------------------------------------------------------------------


    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % 2.ai look at the error in Stern review calculations
    % lower, actual, and upper bounds
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    if new_figs:
    	fig2,axes=plt.subplots()
    else:
    	plt.figure(fig2.number)
    	
    # cost of ships - 5 billion, plus 10% of cost per year - taken from
    # literature
    cost=cost_fac*frac_tot * time1
    for j in range(3):
        # loss to society, which is dependent on temperature rise
        if (j==0):  # estimated from Stern review - lower bound
                loss_in_GDP=(0.5*(tsurface-tsurface[0,0])+2)* \
                    1e12/100;
                col='r';lin=1;ls1='--'
        elif(j==1):     # estimated from Stern review - upper bound
                loss_in_GDP=(5*(tsurface-tsurface[0,0])-15.)* \
                    1e12/100;
                col='r';lin=1;ls1='--'
        elif(j==2):     # estimated from Stern review - center
                loss_in_GDP=(2*(tsurface-tsurface[0,0])-3.5)* \
                    1e12/100;
                loss_in_GDP=scenario(tsurface-tsurface[0,0], kind='central')* \
                    1e12/100;
                col='r';lin=3;ls1='-'
                
        # this is loss every year
        loss_in_GDP=np.maximum(loss_in_GDP,0)*time1;
    
    if method==RAYLEIGH_JET:
    	# Rayleigh jet instability - equation A8
    	power_of_sprayers=(0.45/(D_aerosol*1e-9)+3.2e6)*Qsea;
    elif method==TAYLOR_CONE:
    	# Taylor cone jets - equation A12
    	Qs=5.6e-10;
    	L=1e-3;
    	sig=72e-3;
    	power_of_sprayers= \
    	(128.*8.9e-4*L*Qs/(np.pi*16e-6**4) + 4.*sig/16e-6+1.2e9)*Qsea;
    elif method==SUPERCRITICAL:
    	# Supercritical flow - equation A15
    	power_of_sprayers=4e9*Qsea;
    elif method==EFFERVESCENT:
    	# Effervescent spraying - equation A17
    	power_of_sprayers=3.4e8*Qsea;
    	
    # total cost of the scheme:
    cost=cost+power_of_sprayers*cost_energy*time1*86400*365.25;
    
    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % this calculates the point where the curve starts to increase quickly
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
#         r,c,=np.shape(loss_in_GDP);
#         costval=np.zeros(r);
#         co2val=np.zeros(r);
#         for i in range(r):
#            ind,=np.where(loss_in_GDP[i,:]+cost[i,:]== \
#                np.min(loss_in_GDP[i,:]+cost[i,:]));
#        
#            co2val[i]=co2ppm[i,ind[-1]];
#            costval[i]=cost[i,ind[-1]];
#         #----------------------------------------------------------------------
# 
#         plt.plot(costval,co2val,color=col,linewidth=lin,linestyle=ls1);
#     
#     plt.xlabel('Cost of implementation (dollars)')
#     plt.ylabel('CO_2 mixing ratio (ppm)')
#     plt.title('Point at which the gross cost drastically increases with errors')
#     plt.savefig('/tmp/' + username + '/co2_vs_cost_sweet_spot.png')
    #--------------------------------------------------------------------------
    Z=loss_in_GDP+cost
    ind,=np.where(co2ppm[0,:]<=800)
    ind=ind[-1]
    if new_figs:
	    axC=axes
    axC.plot(cost[:,ind],Z[:,ind])
    if new_figs:
	    axT = axC.twinx()
    axT.plot(cost[:,ind],tsurface[:,ind]-tsurface[0,0],'--')

    axC.set_xlabel('Cost of implementation (dollars)')
    axC.set_ylabel('Cost to society (including negative impact of $\\Delta T$ dollars)')
    axT.set_ylabel('$\\Delta T$ (global temp, dashed)')
    axC.set_yscale('log')
    #axC.legend(['Copper Eff','Harrison Eff','Harrison de Laval','Edmund','Rayleigh','Taylor Cone'])
    plt.savefig('/tmp/' + username + '/co2_doubling.png')
    plt.title('Doubling CO$_2$ to 800 ppm')
	
    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % color plot of cost
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    fig3=plt.figure() #('name','cost vs co2 vs delta T');
    plt.pcolor(cost,co2ppm,Z, \
        norm=colors.LogNorm());
    h=plt.colorbar();
#     h.set_ticks(np.mgrid[8:14:1])
#     h.set_ticklabels(10**np.mgrid[8:14:1])
#     set(h,'ytick',8:1:13,'yticklabel',10.^(8:1:13));
    plt.xlabel('Cost of implementation (dollars)')
    plt.ylabel('CO_2 mixing ratio (ppm)')
    if method==RAYLEIGH_JET:
        plt.title('Rayleigh Jet Spraying')
    elif method==TAYLOR_CONE:
        plt.title('Taylor Cone Spraying')
    elif method==SUPERCRITICAL:
        plt.title('Supercritical Spraying')
    elif method==EFFERVESCENT:
        plt.title('Effervescent Spraying')
        
    h.set_label('Cost of implementing+cost to society (dollars)')
    plt.savefig('/tmp/' + username + '/co2_vs_cost_vs_gross_cost.png')
    #--------------------------------------------------------------------------


