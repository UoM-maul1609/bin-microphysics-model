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

username=getpass.getuser()

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
    f=1.0;
    D_aerosol=45.;               # size of the aerosol particles generated (nm)
    
    if(lut_flag):
        (lut1,lut2)=lutMCB.doAnalysis()
        D_aerosol=lutMCB.Dm*1e9
    
    frac_tot1=f*17.5/100;     # fraction of earth geoengineered
    num_ship=f*1500;           # 1500 ships needed for 17.5% seeding
    time1=50;                   # implement for 50 years
    cost_energy=3.9e-8;         # cost of energy - dollars per joule 

    RAYLEIGH_JET=1;
    TAYLOR_CONE=2;
    SUPERCRITICAL=3;
    EFFERVESCENT=4;

    method=RAYLEIGH_JET;


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
    plt.figure() #(name='albedo relationship');
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
    albedo_geog=0.44;

    num_ship=frac_tot*1500;
    # in order to achieve albedo_geog we need to have the following mr
    if lut_flag:
        scint=interp1d(A,mr1,fill_value='extrapolate')
    else:
        scint=interp1d(A,mr,fill_value='extrapolate')
    
    mr=scint(albedo_geog)


    # now calculate the flow rate
    h=1000;                         # height of the boundary layer (meters)
    tau=3.*86400;                   # time scale aerosol fall out of atmosphere
                                    # 3 days
    Q=4.*np.pi*6.4e6**2.*frac_tot*17.5/100.*h*mr*1.2/(tau);  # equation B2
                                    # kg of NaCl per second

    Qsea=(Q)/(35);                 # m^3 of sea water per second


    tsurface = greenhouse_func_geo.surface_temp(albedo_geog,frac_tot,co2ppm )
    #tsurface=albedo_geo(albedo_geog,frac_tot,co2ppm);   # call climate model 
                                                        # for all co2s
    #--------------------------------------------------------------------------


    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % 2.ai look at the error in Stern review calculations
    % lower, actual, and upper bounds
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    plt.figure() #('name','cost vs co2 vs delta T');

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
                col='r';lin=3;ls1='-'
                
        loss_in_GDP=np.maximum(loss_in_GDP,0);
    
    
    
        # cost of ships - 5 billion, plus 10% of cost per year - taken from
        # literature
        cost=4e9*(frac_tot+0.025*time1*frac_tot);
#         cost=0.
    
    
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
        cost=cost+power_of_sprayers*cost_energy*time1*86400*365; 

    
    
        """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % this calculates the point where the curve starts to increase quickly
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        """
        r,c,=np.shape(loss_in_GDP);
        costval=np.zeros(r);
        co2val=np.zeros(r);
        for i in range(r):
           ind,=np.where(loss_in_GDP[i,:]+cost[i,:]== \
               np.min(loss_in_GDP[i,:]+cost[i,:]));
       
           co2val[i]=co2ppm[i,ind[-1]];
           costval[i]=cost[i,ind[-1]];
        #----------------------------------------------------------------------

        plt.plot(costval,co2val,color=col,linewidth=lin,linestyle=ls1);
    
    plt.xlabel('Cost of implementation (dollars)')
    plt.ylabel('CO_2 mixing ratio (ppm)')
    plt.title('Point at which the gross cost drastically increases with errors')
    plt.savefig('/tmp/' + username + '/co2_vs_cost_sweet_spot.png')
    #--------------------------------------------------------------------------



    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % color plot of cost
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    plt.figure() #('name','cost vs co2 vs delta T');
    Z=loss_in_GDP+cost
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


