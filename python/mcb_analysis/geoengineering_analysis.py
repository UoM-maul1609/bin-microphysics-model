import matplotlib.pyplot as plt
import numpy as np


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
    f=1.0;
    D_aerosol=45;               # size of the aerosol particles generated (nm)
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
    beta=[0.356539140609705   1.435442271779103   ...
        0.000831243139779   7.202626231940690  -0.002074725642290]; 

    mr=np.logspace(-14,-4,100)
    plt.figure(name='albedo relationship');
    A=cloud_albedo(mr,D,beta);
    pcolor(mr,D,A);shading interp;
    h=colorbar;
    set(gca,'yscale','log','xscale','log');
    ylabel('NaCl mode diameter (nm)');
    xlabel('NaCl mixing ratio')
    ylabel(h,'albedo of geoengineered clouds')
    print -dpng images/albedo_vs_NaCl_and_diameter.png
    #--------------------------------------------------------------------------



    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % 2.a climate calculation 
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    [frac_tot,co2ppm]=meshgrid(0:0.001:1,400:10:1500);
    albedo_geog=0.44;

    num_ship=frac_tot.*1500;
    % in order to achieve albedo_geog=0.5 we need to have the following mr
    warning off;
    mr=griddata(A,D,mr,albedo_geog,D_aerosol);
    warning on;


    % now calculate the flow rate
    h=1000;                         % height of the boundary layer (meters)
    tau=3.*86400;                   % time scale aerosol fall out of atmosphere
                                    % 3 days
    Q=4.*pi.*6.4e6.^2.*frac_tot.*17.5./100.*h.*mr.*1.2./(tau);  % equation B2
                                    % kg of NaCl per second

    Qsea=(Q)./(35);                 % m^3 of sea water per second


    tsurface=albedo_geo(albedo_geog,frac_tot,co2ppm);   % call climate model 
                                                        % for all co2s
    #--------------------------------------------------------------------------


    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % 2.ai look at the error in Stern review calculations
    % lower, actual, and upper bounds
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    figure('name','cost vs co2 vs delta T');

    for j=1:3
        # loss to society, which is dependent on temperature rise
        switch (j)
            case (1)    # estimated from Stern review - lower bound
                loss_in_GDP=(217.32.*(tsurface-tsurface(1,1))-543.3).* ...
                    1e12./100;
                col='r--';lin=1;
            case(2)     # estimated from Stern review - upper bound
                loss_in_GDP=(72.44.*(tsurface-tsurface(1,1))-36.22).* ...
                    1e12./100;
                col='r--';lin=1;
            case(3)     # estimated from Stern review - center
                loss_in_GDP=(144.88.*(tsurface-tsurface(1,1))-253.54).* ...
                    1e12./100;
                col='r-';lin=3;
        end
        loss_in_GDP=max(loss_in_GDP,0);
    
    
    
        # cost of ships - 5 billion, plus 10% of cost per year - taken from
        # literature
        cost=5e9.*frac_tot+0.1.*5e9.*time1.*frac_tot;

    
    
        switch (method)
            case RAYLEIGH_JET
                # Rayleigh jet instability - equation A8
                power_of_sprayers=(0.45/(D_aerosol.*1e-9)+3.2e6).*Qsea;

            case TAYLOR_CONE
                # Taylor cone jets - equation A12
                Qs=5.6e-10;
                L=1e-3;
                sig=72e-3;
                power_of_sprayers=...
                    (128.*8.9e-4.*L.*Qs/(pi.*16e-6^4) + 4.*sig/16e-6+1.2e9).*Qsea;

            case SUPERCRITICAL
                # Supercritical flow - equation A15
                power_of_sprayers=4e9.*Qsea;

            case EFFERVESCENT
                # Effervescent spraying - equation A17
                power_of_sprayers=3.4e8.*Qsea;
    
        end    
    
    
    

    
        # total cost of the scheme:
        cost=cost+power_of_sprayers.*cost_energy.*time1.*86400.*365; 

    
    
        """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % this calculates the point where the curve starts to increase quickly
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        """
        [r,c]=size(loss_in_GDP);
        costval=zeros([r,1]);
        co2val=zeros([r,1]);
        for i=1:r
           ind=find(loss_in_GDP(i,:)+cost(i,:)== ...
               min(loss_in_GDP(i,:)+cost(i,:)));
       
           co2val(i)=co2ppm(i,ind);
           costval(i)=cost(i,ind);
        end
        %----------------------------------------------------------------------

        plot(costval,co2val,col,'linewidth',lin); hold on;
    end
    xlabel('Cost of implementation (dollars)')
    ylabel('CO_2 mixing ratio (ppm)')
    title('Point at which the gross cost drastically increases with errors')
    print -dpng images/co2_vs_cost_sweet_spot.png
    #--------------------------------------------------------------------------



    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % color plot of cost
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    figure('name','cost vs co2 vs delta T');
    pcolor(cost,co2ppm,real(log10(loss_in_GDP+cost)));shading interp;
    h=colorbar;
    set(h,'ytick',8:1:13,'yticklabel',10.^(8:1:13));
    xlabel('Cost of implementation (dollars)')
    ylabel('CO_2 mixing ratio (ppm)')
    ylabel(h,'Cost of implementing+cost to society (dollars)')
    print -dpng images/co2_vs_cost_vs_gross_cost.png
    #--------------------------------------------------------------------------



