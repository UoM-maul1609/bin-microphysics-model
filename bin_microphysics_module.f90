	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>code to allocate arrays, and call activation 
	module bmm
    use numerics_type
    use numerics, only : find_pos, poly_int
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the bin microphysics model

    implicit none
        ! constants for the bin microphysics model

        real(wp), parameter :: r_gas=8.314_wp, molw_a=29.e-3_wp,molw_water=18.e-3_wp, &
                                cp=1005.0_wp, cpv=1870._wp, cpw=4.27e3_wp, cpi=2104.6_wp, &
                                grav=9.81_wp, &
        						lv=2.5e6_wp, ls=2.837e6_wp, lf=ls-lv, ttr=273.15_wp, &
        						joules_in_an_erg=1.0e-7_wp,joules_in_a_cal=4.187e0_wp, &
        						rhow=1000._wp, ra=r_gas/molw_a,rv=r_gas/molw_water , &
        						eps1=ra/rv, rhoice=910._wp, &
        						mass_fragment1=pi/6._wp*rhoice*10.e-6_wp**3._wp, &
        						mass_fragment2=mass_fragment1, &
        						mass_fragment3=mass_fragment1
        						

        type parcel
            ! variables for bin model
            integer(i4b) :: n_bins1,n_modes,n_comps, n_bin_mode, n_bin_modew, n_bin_mode1, &
                            n_sound, ice_flag, sce_flag
            real(wp) :: dt
            real(wp), dimension(:,:), allocatable :: q_sound
            real(wp), dimension(:), allocatable :: t_sound, z_sound, rh_sound, &
                                                    p_sound, theta_q_sound
            real(wp) :: z,p,t,w,rh, qinit, t_cbase, q_cbase, p_cbase, z_cbase, &
                        t_ctop, q_ctop, p_ctop, z_ctop, theta_q_cbase, theta_q_ctop, &
                        x_ent, theta_q
                        
                        
            ! liquid water
            real(wp), dimension(:), allocatable :: d, maer, npart, rho_core, &
                            rh_eq, rhoat, dw, da_dt, ndrop, npartall
            real(wp), dimension(:,:), allocatable :: mbin, mbinall, rhobin, &
                                        nubin,molwbin,kappabin ! all bins x all comps                                
            ! variables for ODE:                    
            integer(i4b) :: neq, itol, ipr, ite, iz, iw, irh, &
                            itask, istate, iopt, mf, lrw, liw
            integer(i4b), dimension(:), allocatable :: iwork
            integer(i4b), dimension(1) :: ipar
            real(wp) :: tt, tout
            real(wp), dimension(:), allocatable :: y, yold, atol, rwork
            real(wp), dimension(1) :: rpar
            real(wp), dimension(1) :: rtol
            
            ! ice water
            real(wp), dimension(:), allocatable :: dice, maerice, npartice, rho_coreice, &
                            rh_eqice, rhoatice, dwice, da_dtice, nice, &
                            phi, rhoi, nump, rime
            real(wp), dimension(:,:), allocatable :: mbinice, rhobinice, &
                                        nubinice,molwbinice,kappabinice ! all bins x all comps  
                                        
                                        
            ! general
            integer(i4b) :: imoms
            real(wp), allocatable, dimension(:,:) :: moments, mbinedges,ecoll,ecoal
            real(wp), allocatable, dimension(:) :: momtemp, vel
            integer(i4b), allocatable, dimension(:) :: momenttype
            integer(i4b), dimension(:,:), allocatable :: indexc                            
                                          
            ! variables for ODE:                    
            integer(i4b) :: neqice, itolice, ipri, itei, izi, iwi, irhi, &
                            itaskice, istateice, ioptice, mfice, lrwice, liwice
            integer(i4b), dimension(:), allocatable :: iworkice
            integer(i4b), dimension(1) :: iparice
            real(wp) :: ttice, toutice, totaddto
            real(wp), dimension(1) :: rtolice
            real(wp), dimension(:), allocatable :: yice, yoldice, atolice, rworkice
            real(wp), dimension(1) :: rparice
            
            
            logical :: break_flag=.false.
            
        end type parcel

        type sounding
            ! variables for grid
            integer(i4b) :: n_levels
            real(wp), dimension(:,:), allocatable :: q
            real(wp), dimension(:), allocatable :: theta, p, z, rh
        end type sounding


        type io
            ! variables for io
            integer(i4b) :: ncid, varid, x_dimid, bin_dimid, bin2_dimid, bin3_dimid, &
                            mode_dimid, comp_dimid, y_dimid, z_dimid, &
                            dimids(2), a_dimid, xx_dimid, yy_dimid, &
                            zz_dimid, i_dimid, j_dimid, k_dimid, nq_dimid, nprec_dimid
            integer(i4b) :: icur=1
            logical :: new_file=.true.
        end type io


        ! declare a parcel type
        type(parcel) :: parcel1
        ! declare a sounding type
        type(sounding) :: sounding1
        ! declare an io type
        type(io) :: io1
        

        ! some namelist variables
        logical :: micro_init=.true., adiabatic_prof=.false., vert_ent=.false.
        real(wp) :: ent_rate, dmina,dmaxa
        real(wp) :: zinit,tpert,winit,winit2, amplitude2, tau2, &
                    tinit,pinit,rhinit,z_ctop, alpha_therm, alpha_cond, &
                    alpha_therm_ice, alpha_dep
        integer(i4b) :: microphysics_flag=0, kappa_flag,updraft_type, vent_flag, &
                        sce_flag=0,ice_flag=0, bin_scheme_flag=1
        logical :: use_prof_for_tprh, hm_flag, mode1_flag, mode2_flag
        integer(i4b) :: break_flag
        real(wp) :: dz,dt, runtime, t_thresh
        ! sounding spec
        real(wp) :: psurf, tsurf
        integer(i4b), parameter :: nlevels_r=1000
        integer(i4b), parameter :: nq=3
        integer(i4b) :: n_levels_s, idum, n_sel
        real(wp) :: mult, rh_act
        real(wp), allocatable, dimension(:,:) :: q_read !nq x nlevels_r
        real(wp), allocatable, dimension(:) :: theta_read,rh_read,  z_read
        ! aerosol setup
        integer(i4b) :: n_intern, n_mode,n_sv,sv_flag,n_bins,n_comps
        ! aerosol_spec
        real(wp), allocatable, dimension(:,:) :: n_aer1,d_aer1,sig_aer1, mass_frac_aer1
        real(wp), allocatable, dimension(:) ::  molw_core1,density_core1,nu_core1, &
                                        kappa_core1
        real(wp), allocatable, dimension(:) :: org_content1, molw_org1, kappa_org1, &
                                    density_org1, delta_h_vap1,nu_org1,log_c_star1
        
        

        ! variables for model
        real(wp) :: theta_surf,theta_init, &
            theta_q_sat,t1old, p111, w_cb, n_dummy, d_dummy, x2old=1.0_wp
        logical :: set_theta_q_cb_flag=.true.

        ! Chen and Lamb (1994) Gamma variable fit (scaled and centred logarithm)
        integer(i4b), parameter :: n_cl=18
        real(wp), dimension(n_cl), parameter :: gam_cl=[-0.072328469664620_wp, &
            -0.324623262465577_wp, 0.363138099937540_wp, 3.323089908344732_wp, &
            0.874844989423720_wp, &
            -13.554426432462339_wp, -9.810322482346461_wp, 27.846739088352344_wp, &
            26.480447842355410_wp,&
             -29.890199206698309_wp, -32.327548996894521_wp, 15.827423311652167_wp, &
             18.466605783503052_wp, -4.158566361058538_wp, -5.039533848938808_wp, &
             1.477272813054374_wp, 1.038600921563425_wp, -0.457007828432810_wp]
        real(sp), dimension(2), parameter :: gam_mu_cl=[260.163817050062335_wp, &
                                                    8.274747821396463_wp]


        character (len=200) :: outputfile='output', scefile='input'

	!private 
	!public :: read_in_bmm_namelist, initialise_bmm_arrays, bmm_driver, io1
    public
    private :: outputfile
	contains	
	
		
		
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>allocate arrays for activation code
	!>@param[in] n_intern: number of aerosol modes of same kind
	!>@param[in] n_mode: number of aerosol modes
	!>@param[in] n_sv: number of organic / volatility modes
	!>@param[in] n_bins: number of size bins in a mode
	!>@param[in] n_comps: number of different compositions in a mode
	!>@param[in] nq: number of q-variables in sounding
	!>@param[in] n_levels_s: number of levels in sounding
	!>@param[inout] q_read, theta_read, rh_read, z_read: sounding	
	!>@param[inout] n_aer1: number conc. in modes
	!>@param[inout] d_aer1: diameter in modes
	!>@param[inout] sig_aer1: geo std in modes
	!>@param[inout] mass_frac_aer1:mass_fraction of each component
	!>@param[inout] molw_core1:molw in core
	!>@param[inout] density_core1: solute density
	!>@param[inout] nu_core1: van hoff factor
	!>@param[inout] kappa_core1: kappa parameter
	!>@param[inout] org_content1: organic content in vol bins
	!>@param[inout] molw_org1: molw in volatility bins
	!>@param[inout] kappa_org1: kappa in volatility bins
	!>@param[inout] density_org1: density in volatility bins
	!>@param[inout] delta_h_vap1: enthalpy in volatility bins
	!>@param[inout] nu_org1: van hoff factor in volatility bins
	!>@param[inout] log_c_star1: log_c_star in volatility bins
	subroutine allocate_arrays(n_intern,n_mode,n_sv,n_bins,n_comps,nq,n_levels_s, &
		                    q_read,theta_read,rh_read,z_read, &
		                    n_aer1,d_aer1,sig_aer1,mass_frac_aer1, molw_core1, &
		                    density_core1, nu_core1, kappa_core1, &
		                    org_content1,molw_org1,kappa_org1,density_org1, &
		                    delta_h_vap1,nu_org1,log_c_star1)
		use numerics_type
		implicit none
		integer(i4b), intent(in) :: n_intern, n_mode, n_sv, n_bins,n_comps, nq, &
		                            n_levels_s
		real(wp), dimension(:), allocatable, intent(inout) :: theta_read,rh_read,z_read, &
		                        org_content1,molw_org1,kappa_org1, &
		                        density_org1,delta_h_vap1,nu_org1,log_c_star1
		real(wp), dimension(:,:), allocatable, intent(inout) :: q_read, &
		                        n_aer1,d_aer1,sig_aer1,mass_frac_aer1
		real(wp), dimension(:), allocatable, intent(inout) :: molw_core1,density_core1, &
		                        nu_core1,kappa_core1
		
		integer(i4b) :: AllocateStatus
		allocate( q_read(1:nq,1:n_levels_s), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( theta_read(1:n_levels_s), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( rh_read(1:n_levels_s), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( z_read(1:n_levels_s), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		
		allocate( n_aer1(1:n_intern,1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( d_aer1(1:n_intern,1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( sig_aer1(1:n_intern,1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
	
		allocate( mass_frac_aer1(1:n_mode,1:n_comps), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( molw_core1(1:n_comps), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_core1(1:n_comps), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_core1(1:n_comps), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( kappa_core1(1:n_comps), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

		allocate( org_content1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( molw_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( kappa_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( density_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( delta_h_vap1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( nu_org1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( log_c_star1(1:n_sv), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

	end subroutine allocate_arrays
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! read in the namelist                                                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>read in the data from the namelists for the bin microphysics model
	!>@param[in] nmlfile
	subroutine read_in_bmm_namelist(nmlfile)
		implicit none
        character (len=200), intent(in) :: nmlfile
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /sounding_spec/ psurf, tsurf,  &
                    q_read, theta_read, rh_read, z_read
        ! define namelists for environment
        namelist /run_vars/ outputfile, scefile,runtime, dt, &
                    zinit,tpert,use_prof_for_tprh,winit,winit2,amplitude2, &
                    tinit,pinit,rhinit, &
                    microphysics_flag, ice_flag, bin_scheme_flag, sce_flag, &
                    hm_flag, break_flag, mode1_flag, mode2_flag, vent_flag, &
                    kappa_flag, updraft_type,t_thresh, adiabatic_prof, vert_ent, &
                    z_ctop, ent_rate,n_levels_s, &
                    alpha_therm,alpha_cond,alpha_therm_ice,alpha_dep
        namelist /aerosol_setup/ n_intern,n_mode,n_sv,sv_flag, n_bins,n_comps
        namelist /aerosol_spec/ n_aer1,d_aer1,sig_aer1, dmina,dmaxa, &
                                mass_frac_aer1, molw_core1, &
                                density_core1,nu_core1,kappa_core1, &
                                org_content1, molw_org1,kappa_org1, &
                                density_org1, delta_h_vap1,nu_org1, log_c_star1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists	and allocate arrays								   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        read(8,nml=aerosol_setup)
        ! allocate memory / init
		call allocate_arrays(n_intern,n_mode,n_sv,n_bins,n_comps,nq,n_levels_s, &
		                    q_read,theta_read,rh_read,z_read, &
		                    n_aer1,d_aer1,sig_aer1,mass_frac_aer1, molw_core1, &
		                    density_core1, nu_core1, kappa_core1, &
		                    org_content1,molw_org1,kappa_org1,density_org1, &
		                    delta_h_vap1,nu_org1,log_c_star1)
        
        read(8,nml=sounding_spec)
        read(8,nml=aerosol_spec)
        close(8)
        tau2 = 2._wp*pi/winit2*amplitude2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine read_in_bmm_namelist







	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! initialise arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>interpolates the sounding to the grid
    subroutine initialise_bmm_arrays(psurf, tsurf, q_read, theta_read, rh_read, z_read, &
                    runtime, dt, zinit, tpert, use_prof_for_tprh, winit, tinit, pinit, &
                    rhinit, microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                    kappa_flag, updraft_type, adiabatic_prof, vert_ent, z_ctop, &
                    ent_rate, n_levels_s, alpha_therm, alpha_cond, alpha_therm_ice, &
                    alpha_dep, n_intern, n_mode, n_sv, sv_flag, n_bins, n_comps, &
                    n_aer1,d_aer1,sig_aer1,dmina,dmaxa,mass_frac_aer1,molw_core1, &
                    density_core1, nu_core1, kappa_core1, org_content1, molw_org1, &
                    kappa_org1, density_org1, delta_h_vap1,nu_org1, log_c_star1, &
                    sce_flag)
    use numerics_type
    use numerics, only : find_pos, poly_int, zeroin, fmin,vode_integrate

    implicit none
    logical, intent(in) :: use_prof_for_tprh, adiabatic_prof, vert_ent
    integer(i4b), intent(in) :: microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                    kappa_flag, updraft_type, n_levels_s,n_intern, n_mode, n_sv, &
                    sv_flag, n_bins, n_comps, sce_flag
    real(wp), intent(in) :: psurf, tsurf, runtime, dt, zinit, tpert, winit, tinit, &
                    pinit, rhinit, alpha_therm, alpha_cond, alpha_therm_ice, &
                    alpha_dep, dmina,dmaxa, z_ctop, ent_rate
    real(wp), dimension(1:n_levels_s), intent(in) :: theta_read, rh_read, z_read
    real(wp), dimension(1:nq,1:n_levels_s), intent(in) :: q_read
    real(wp), dimension(1:n_intern,1:n_mode), intent(in) :: n_aer1,d_aer1,sig_aer1
    real(wp), dimension(1:n_mode,1:n_comps), intent(in) :: mass_frac_aer1
    real(wp), dimension(1:n_comps), intent(in) :: molw_core1, density_core1, &
                                                nu_core1, kappa_core1
    real(wp), dimension(1:n_sv), intent(in) :: org_content1, molw_org1, kappa_org1, &
                                density_org1, delta_h_vap1, nu_org1, log_c_star1
    
    
    
    
    real(wp) :: num, ntot, number_per_bin, test, var1, &
                eps2, z1, z2, htry, hmin, var, dummy
    real(wp), dimension(1) :: p1, z11
    real(wp) :: p11, p22, rm, cpm
    integer(i4b) :: i,j,k, AllocateStatus, iloc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set variables and allocate arrays in parcel                                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parcel1%n_sound=n_levels_s   
    parcel1%n_bins1=n_bins    
    parcel1%n_modes=n_mode
    parcel1%n_comps=n_comps
    parcel1%n_bin_modew=n_bins*n_mode
    parcel1%n_bin_mode1=(n_bins+1)*n_mode
    parcel1%z=zinit
    parcel1%p=pinit
    parcel1%t=tinit
    parcel1%w=winit
    parcel1%rh=rhinit
    parcel1%dt=dt

    parcel1%ice_flag=ice_flag
    parcel1%n_bin_mode=&
        parcel1%n_bins1*n_mode*(1+parcel1%ice_flag)     ! for all the liquid and ice    
    parcel1%imoms=ice_flag*5                            ! phi, nmon, vol, rim, unf
    
    allocate( parcel1%d(1:parcel1%n_bin_mode1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    ! same bin edges used for ice
    allocate( parcel1%mbinedges(1:parcel1%n_bins1+1,1:parcel1%n_modes), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%maer(1:parcel1%n_bin_modew), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%npart(1:parcel1%n_bin_modew), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%npartall(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%mbin(1:parcel1%n_bin_modew,1:n_comps+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%mbinall(1:parcel1%n_bin_mode,1:n_comps+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%rho_core(1:parcel1%n_modes), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    
    allocate( parcel1%momtemp(1:parcel1%n_bin_mode), &
        STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%moments(1:parcel1%n_bin_mode,1:n_comps+parcel1%imoms), &
        STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%momenttype(1:n_comps+parcel1%imoms), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    

    allocate( parcel1%rhobin(1:parcel1%n_bin_modew,1:n_comps), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%nubin(1:parcel1%n_bin_modew,1:n_comps), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%molwbin(1:parcel1%n_bin_modew,1:n_comps), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%kappabin(1:parcel1%n_bin_modew,1:n_comps), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

    allocate( parcel1%rh_eq(1:parcel1%n_bin_modew), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%rhoat(1:parcel1%n_bin_modew), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%dw(1:parcel1%n_bin_modew), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%da_dt(1:parcel1%n_bin_modew), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%ndrop(1:parcel1%n_bin_modew), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

    allocate( parcel1%ecoal(1:parcel1%n_bin_mode,1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%ecoll(1:parcel1%n_bin_mode,1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%indexc(1:parcel1%n_bin_mode,1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%vel(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

    allocate( parcel1%q_sound(1:parcel1%n_sound,1:nq), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%z_sound(1:parcel1%n_sound), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%t_sound(1:parcel1%n_sound), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%p_sound(1:parcel1%n_sound), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%rh_sound(1:parcel1%n_sound), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%theta_q_sound(1:parcel1%n_sound), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the density of aerosol particles within a mode                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,n_mode
        var1=sum(mass_frac_aer1(i,:)/ density_core1)
        parcel1%rho_core(i) = 1._wp/var1
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set-up size distribution                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=1,n_mode
        idum=k ! this is sent through to zbrent to select the correct mode
        ! find total number in mode between dmina and dmaxa:
        call lognormal_n_between_limits(n_aer1(:,k),d_aer1(:,k),sig_aer1(:,k), &
                                        n_intern,dmina,dmaxa, num)
        !print *,num
        ! set up variables for parcel model
        ntot=num
        number_per_bin=ntot/real(n_bins,wp)
        parcel1%npart(1+(k-1)*n_bins:(k)*n_bins)=number_per_bin
        parcel1%d(1+(k-1)*(n_bins+1))=dmina
        do i=1,n_bins
            d_dummy=parcel1%d(i+(k-1)*(n_bins+1))
            n_dummy=number_per_bin*(1._wp-1.e-5_wp)
            parcel1%d(i+1+(k-1)*(n_bins+1))= zeroin(&
                        d_dummy*0.9_wp,dmaxa*2._wp,find_upper_diameter, 1.e-30_wp)
        enddo
        parcel1%d((k)*(n_bins+1))=dmaxa ! nail it to end point - round off
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! aerosol mass - total                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=1,parcel1%n_modes
        do j=1,parcel1%n_bins1
            i=j+(k-1)*(n_bins+1)
            parcel1%maer(j+(k-1)*(n_bins))= &
                pi/6._wp*(0.5_wp*(parcel1%d(i+1)+parcel1%d(i)))**3 * &
                parcel1%rho_core(k)
        enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the mass of each component in a bin, including water (Koehler eq)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,parcel1%n_bins1
        do j=1,parcel1%n_modes
            do k=1,parcel1%n_comps
                parcel1%mbin(i+(j-1)*n_bins,k)= &  
                        parcel1%maer(i+(j-1)*n_bins)*mass_frac_aer1(j,k)
                ! density in each bin:
                parcel1%rhobin(i+(j-1)*n_bins,k)=density_core1(k)
                ! nu in each bin:
                parcel1%nubin(i+(j-1)*n_bins,k)=nu_core1(k)
                ! molw in each bin:
                parcel1%molwbin(i+(j-1)*n_bins,k)=molw_core1(k)
                ! kappa in each bin:
                parcel1%kappabin(i+(j-1)*n_bins,k)=kappa_core1(k)
            enddo
        enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! get initial conditions (based on sounding)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! integrate dp/dz=-p/ra/t
    ! dp/dz=-p/(ra*theta*(p/100000.)**0.287)
    parcel1%q_sound=q_read    
    parcel1%z_sound=z_read  
    parcel1%t_sound(1)=theta_read(1)*(psurf/1.e5_wp)**(ra/cp) 
    parcel1%p_sound(1)=psurf
    
    eps2=1.e-5_wp
    htry=10._wp
    hmin=1.e-2_wp
    p1=psurf

    ! integrate hydrostatic equation
    do i=2,n_levels_s
        call vode_integrate(p1,z_read(i-1),z_read(i),eps2,htry,hmin,hydrostatic1b)
        parcel1%t_sound(i)=theta_read(i)*(p1(1)/1.e5_wp)**(ra/cp)  
        parcel1%p_sound(i)=p1(1)      
    enddo
    
    ! calculate RH
    do i=1,n_levels_s
        parcel1%rh_sound(i)=parcel1%q_sound(1,i)/ &
            (eps1*svp_liq(parcel1%t_sound(i)) / &
            (parcel1%p_sound(i)- svp_liq(parcel1%t_sound(i))) )
    enddo
    
    ! calculate theta_q
    do i=1,n_levels_s
        parcel1%theta_q_sound(i)= &
            calc_theta_q3(parcel1%t_sound(i),parcel1%p_sound(i),parcel1%q_sound(1,i))
        
    enddo
    ! interpolate to find parcel conditions
    if (use_prof_for_tprh) then
        ! interpolate to find theta
        iloc=find_pos(parcel1%z_sound(1:n_levels_s),parcel1%z)
        iloc=min(n_levels_s-1,iloc)
        iloc=max(1,iloc)
        ! linear interp t
        call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%t_sound(iloc:iloc+1), &
                    min(parcel1%z,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%t=var +tpert
        ! linear interp rh
        call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%p_sound(iloc:iloc+1), &
                    min(parcel1%z,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%p=var 
        ! linear interp rh
        call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%rh_sound(iloc:iloc+1), &
                    min(parcel1%z,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%rh=var 
        print *,'t,p,rh from sounding: ', parcel1%t, parcel1%p, parcel1%rh
    endif
    parcel1%qinit=parcel1%rh*eps1*svp_liq(parcel1%t)/(parcel1%p-svp_liq(parcel1%t))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! get cloud-base conditions (for entrainment process)                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (vert_ent) then
	    ! find cloud-base pressure:
	    ! cloud base qv
	    parcel1%q_cbase=parcel1%qinit
		! calculate the dry adiabat:
		theta_init=parcel1%t*(1.e5_wp/parcel1%p)**(ra/cp)
		! calculate p required so that qs(t,p) = q init
		parcel1%p_cbase=zeroin(parcel1%p, 100._wp,cloud_base, 1.e-30_wp)
		rm=ra+parcel1%qinit*rv
		cpm=cp+parcel1%qinit*cpv
		parcel1%t_cbase=theta_init*(parcel1%p_cbase/1.e5_wp)**(rm/cpm)
		parcel1%theta_q_cbase= &
		    calc_theta_q3(parcel1%t_cbase,parcel1%p_cbase,parcel1%q_cbase)
		print *,'Cloud-base t, p: ',parcel1%t_cbase,parcel1%p_cbase
		
		
		
		! now, find the height of cloud-base along dry adiabat:
		theta_surf=tsurf*(1.e5_wp/psurf)**(ra/cp)
		p11=parcel1%p
		z11(1)=parcel1%z
		p22=parcel1%p_cbase
		htry=p22-p11
		eps2=1.e-5_wp
		call vode_integrate(z11,p11,p22,eps2,htry,hmin,hydrostatic1)
		parcel1%z_cbase=z11(1)
        ! cloud-top properties from sounding:
        parcel1%z_ctop=z_ctop
        ! interpolate to find t
        iloc=find_pos(parcel1%z_sound(1:n_levels_s),parcel1%z_ctop)
        iloc=min(n_levels_s-1,iloc)
        iloc=max(1,iloc)
        ! linear interp t
        call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%t_sound(iloc:iloc+1), &
                    min(parcel1%z_ctop,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%t_ctop=var
        ! linear interp p
        call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%p_sound(iloc:iloc+1), &
                    min(parcel1%z_ctop,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%p_ctop=var
        ! linear interp q
        call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%rh_sound(iloc:iloc+1), &
                    min(parcel1%z_ctop,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%q_ctop=var*eps1*svp_liq(parcel1%t_ctop) / &
            (parcel1%p_ctop-svp_liq(parcel1%t_ctop))
		parcel1%theta_q_ctop= &
		    calc_theta_q3(parcel1%t_ctop,parcel1%p_ctop,parcel1%q_ctop)
		print *,'Cloud-top t, p: ',parcel1%t_ctop,parcel1%p_ctop
	endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! put water on bin, using koehler equation                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    select case(kappa_flag)
        case(0)
            do i=1,parcel1%n_bin_modew
                n_sel=i
                rh_act=0._wp !min(parcel1%rh,0.999_wp)
                mult=-1._wp
                ! has to be less than the peak moles of water at activation
                test=fmin(1.e-50_wp,1.e1_wp, koehler02,1.e-30_wp)
                rh_act=min(parcel1%rh,0.999_wp)
                mult=1._wp
                d_dummy=zeroin(1.e-30_wp, test, koehler02,1.e-30_wp)*molw_water 
                parcel1%mbin(i,n_comps+1)= d_dummy
            enddo
!             call koehler01(parcel1%t,parcel1%mbin(:,n_comps+1),&
!                 parcel1%mbin,parcel1%rhobin,&
!                 parcel1%nubin,parcel1%molwbin,parcel1%n_bin_modew,&
!                 parcel1%rh_eq,parcel1%rhoat,parcel1%dw) 
!             print *,parcel1%rh_eq
        case(1)
            do i=1,parcel1%n_bin_modew
                n_sel=i
                rh_act=0._wp !min(parcel1%rh,0.999_wp)
                mult=-1._wp
                ! has to be less than the peak moles of water at activation
                test=fmin(1.e-50_wp,1.e1_wp, kkoehler02,1.e-30_wp)
                rh_act=min(parcel1%rh,0.999_wp)
                mult=1._wp
                d_dummy=zeroin(1.e-30_wp, test, kkoehler02,1.e-30_wp)*molw_water 
                parcel1%mbin(i,n_comps+1)= d_dummy
            enddo
!             call kkoehler01(parcel1%t,parcel1%mbin(:,n_comps+1),&
!                 parcel1%mbin,parcel1%rhobin,&
!                 parcel1%kappabin,parcel1%molwbin,parcel1%n_bin_modew,&
!                 parcel1%rh_eq,parcel1%rhoat,parcel1%dw) 
!             print *,parcel1%rh_eq
        case default
            print *,'error kappa flag'
            stop
    end select
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set-up ODE variables                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parcel1%neq=parcel1%n_bin_modew+5 ! p,t,rh,z,w
    parcel1%tt=0._wp
    parcel1%tout=parcel1%tt+parcel1%dt
    parcel1%itol=2
    parcel1%rtol=1.e-4_wp
    allocate( parcel1%y(parcel1%neq), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"
    allocate( parcel1%yold(parcel1%neq), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"
    allocate( parcel1%atol(parcel1%neq), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"
    
    parcel1%atol(1:parcel1%n_bin_modew)=1.e-25_wp
    
    parcel1%ipr=parcel1%n_bin_modew+1 ! pressure
    parcel1%ite=parcel1%n_bin_modew+2 ! temperarture
    parcel1%irh=parcel1%n_bin_modew+3 ! rh
    parcel1%iz =parcel1%n_bin_modew+4 ! altitude
    parcel1%iw =parcel1%n_bin_modew+5 ! vertical wind
    
    parcel1%atol(parcel1%ipr)=10._wp
    parcel1%atol(parcel1%ite)=1.e-4_wp
    parcel1%atol(parcel1%irh)=1.e-8_wp
    parcel1%atol(parcel1%iz) =2.e-2_wp
    parcel1%atol(parcel1%iw) =2.e-2_wp
    
    if(parcel1%iw .ne. parcel1%neq) stop "*** problem with array lengths ***"
    parcel1%itask=1
    parcel1%istate=1
    parcel1%iopt=1
    parcel1%mf=22
    
    parcel1%lrw=22+9*parcel1%neq+2*parcel1%neq**2
    allocate( parcel1%rwork(parcel1%lrw), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"

    parcel1%liw=30+parcel1%neq
    allocate( parcel1%iwork(parcel1%liw), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"
    
    ! extra input variables:
    parcel1%iwork=0
    parcel1%rwork=0._wp
    parcel1%iwork(6) = 100 ! max steps
    parcel1%iwork(7) = 10 ! max message printed per problem
    parcel1%iwork(5) = 5 ! order
    parcel1%rwork(5) = 0._wp !1.e-3_wp ! initial time-step
    parcel1%rwork(6) = dt ! max time-step
    parcel1%rwork(7) = 0._wp !1.e-9_wp ! min time-step allowed
    parcel1%rwork(14) = 2._wp ! tolerance scale factor
    
    ! put water in solution vector and set p, t, rh, z, w
    parcel1%y(1:parcel1%n_bin_modew)=parcel1%mbin(:,n_comps+1)
    parcel1%y(parcel1%ipr)=parcel1%p
    parcel1%y(parcel1%ite)=parcel1%t
    parcel1%y(parcel1%irh)=parcel1%rh
    parcel1%y(parcel1%iz) =parcel1%z
    parcel1%y(parcel1%iw) =parcel1%w
    ! do not print messages
    call xsetf(0)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(ice_flag .eq. 1) then
        ! allocation:
        allocate( parcel1%dice(1:parcel1%n_bin_mode1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%maerice(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%npartice(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%mbinice(1:parcel1%n_bin_modew,1:n_comps+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%rho_coreice(1:parcel1%n_modes), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

        allocate( parcel1%rhobinice(1:parcel1%n_bin_modew,1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%nubinice(1:parcel1%n_bin_modew,1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%molwbinice(1:parcel1%n_bin_modew,1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%kappabinice(1:parcel1%n_bin_modew,1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

        allocate( parcel1%rh_eqice(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%rhoatice(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%dwice(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%da_dtice(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%nice(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        
        allocate( parcel1%phi(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%rhoi(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%nump(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"     
        allocate( parcel1%rime(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        
        parcel1%phi=1._wp
        parcel1%rhoi=rhoice
        parcel1%nump=1._wp
        parcel1%rime=0._wp
        
                
        parcel1%rho_coreice(:) = parcel1%rho_core(:)
        
        parcel1%npartice=0._wp
        parcel1%dice=parcel1%d
        parcel1%maerice=parcel1%maer
        parcel1%mbinice=parcel1%mbin
        parcel1%rhobinice=parcel1%rhobin
        parcel1%nubinice=parcel1%nubin
        parcel1%molwbinice=parcel1%molwbin
        parcel1%kappabinice=parcel1%kappabin
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! set-up ODE variables                                                         !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        parcel1%neqice=parcel1%n_bin_modew+4 ! p,t,rh,w
        parcel1%ttice=0._wp
        parcel1%toutice=parcel1%tout
        parcel1%itolice=2
        parcel1%rtolice=1.e-3_wp
        allocate( parcel1%yice(parcel1%neqice), stat = allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory ***"
        allocate( parcel1%yoldice(parcel1%neqice), stat = allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory ***"
        allocate( parcel1%atolice(parcel1%neqice), stat = allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory ***"
    
        parcel1%atolice(1:parcel1%n_bin_modew)=1.e-25_wp
    
        parcel1%ipri=parcel1%n_bin_modew+1 ! pressure
        parcel1%itei=parcel1%n_bin_modew+2 ! temperature
        parcel1%irhi=parcel1%n_bin_modew+3 ! rh
        parcel1%iwi =parcel1%n_bin_modew+4 ! vertical wind
    
        parcel1%atolice(parcel1%ipri)=10._wp
        parcel1%atolice(parcel1%itei)=1.e-4_wp
        parcel1%atolice(parcel1%irhi)=1.e-8_wp
        parcel1%atolice(parcel1%iwi) =2.e-2_wp
    
        if(parcel1%iwi .ne. parcel1%neqice) stop "*** problem with array lengths ***"
        parcel1%itaskice=1
        parcel1%istateice=1
        parcel1%ioptice=1
        parcel1%mfice=22
        parcel1%lrwice=22+9*parcel1%neqice+2*parcel1%neqice**2
        allocate( parcel1%rworkice(parcel1%lrwice), stat = allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory ***"

        parcel1%liwice=30+parcel1%neqice
        allocate( parcel1%iworkice(parcel1%liwice), stat = allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory ***"
    
        ! extra input variables:
        parcel1%iworkice=0
        parcel1%rworkice=0._wp
        parcel1%iworkice(6) = 100 ! max steps
        parcel1%iworkice(7) = 10 ! max message printed per problem
        parcel1%iworkice(5) = 5 ! order
        parcel1%rworkice(5) = 0._wp !1.e-3_wp ! initial time-step
        parcel1%rworkice(6) = dt ! max time-step
        parcel1%rworkice(7) = 0._wp !1.e-9_wp ! min time-step allowed
        parcel1%rworkice(14) = 2._wp ! tolerance scale factor
    
        ! put water in solution vector and set p, t, rh, z, w
        parcel1%yice(1:parcel1%n_bin_modew)=parcel1%mbinice(:,n_comps+1)
        parcel1%yice(parcel1%ipri)=parcel1%p
        parcel1%yice(parcel1%itei)=parcel1%t
        parcel1%yice(parcel1%irhi)=parcel1%rh
        parcel1%yice(parcel1%iwi) =parcel1%w
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise conserved moments                                                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parcel1%moments=0._wp
    do j=1,parcel1%n_comps
        ! above the aerosol
        do i=1,parcel1%n_bin_modew
            ! aerosol moments
            parcel1%moments(i,j)=parcel1%npart(i)*parcel1%mbin(i,j)
        enddo
    enddo
    parcel1%momenttype(1:parcel1%n_comps)=1 ! 1 is mass, 2 is number

    ! ** ice **
    ! same for the ice
    ! additional moments - 2 general ones and 3 just for ice maybe just do 5
    if (parcel1%ice_flag.eq.1) then
    
        do j=1,parcel1%n_comps
            ! above the aerosol
            do i=1,parcel1%n_bin_modew
                ! aerosol moments
                parcel1%moments(i+parcel1%n_bin_modew,j)= &
                    parcel1%npartice(i)*parcel1%mbinice(i,j)
            enddo
        enddo
    
        ! extra ice moments
        do i=1,parcel1%n_bin_modew
            ! ice moments: phi, nmon, vol, rim, unf
            ! phi: 1*n
            parcel1%moments(i+parcel1%n_bin_modew,parcel1%n_comps+1)=parcel1%npartice(i)
            ! nmon: 1*n
            parcel1%moments(i+parcel1%n_bin_modew,parcel1%n_comps+2)=parcel1%npartice(i)
            ! vol: mass/rho
            parcel1%moments(i+parcel1%n_bin_modew,parcel1%n_comps+3)=parcel1%npartice(i)* &
                parcel1%mbinice(i,parcel1%n_comps+1)/rhoice
        enddo  
        
              
        do i=1,parcel1%n_bin_modew
            ! ice moments: phi, nmon, vol, rim, unf
            ! rim: mass
            parcel1%moments(i,parcel1%n_comps+4)=parcel1%npart(i)* &
                parcel1%mbin(i,parcel1%n_comps+1)
            ! unf: mass
            parcel1%moments(i,parcel1%n_comps+5)=parcel1%npart(i)* &
                parcel1%mbin(i,parcel1%n_comps+1)
        enddo        
        parcel1%momenttype(parcel1%n_comps+1:parcel1%n_comps+parcel1%imoms)=[2,2,1,1,1]
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end subroutine initialise_bmm_arrays
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! MAP SCE                                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>maps the sce variables onto BMM
    subroutine write_sce_to_bmm(n_bin_mode,n_bin_modew,n_binst,n_mode, n_comps, n_moments, &
                    ice_flag, &
                    npart, moments, mbin, vel, indexc,ecoll,mbinedges)
    implicit none
    integer(i4b), intent(in) :: n_bin_mode, n_bin_modew, &
        n_binst, n_mode, n_comps, n_moments, ice_flag
    real(wp), dimension(n_bin_mode), intent(in) :: npart,vel
    real(wp), dimension(n_bin_mode,n_moments), intent(in) :: moments
    real(wp), dimension(n_bin_mode,n_comps+1), intent(in) :: mbin
    real(wp), dimension(n_binst+1,n_mode), intent(in) :: mbinedges
    integer(i4b), dimension(n_bin_mode,n_bin_mode), intent(in) :: indexc
    real(wp), dimension(n_bin_mode,n_bin_mode), intent(in) :: ecoll


    
    parcel1%moments=moments
    parcel1%npart=npart(1:n_bin_modew)
    parcel1%npartall=npart
    parcel1%y(1:n_bin_modew)=mbin(1:n_bin_modew,n_comps+1)
    parcel1%maer=sum(mbin(1:n_bin_modew,1:n_comps),2)
    parcel1%mbin=mbin(1:n_bin_modew,1:n_comps+1)
    parcel1%mbinedges=mbinedges
    parcel1%mbinall(:,:)=mbin
    
    parcel1%indexc=indexc
    parcel1%ecoll=ecoll
    parcel1%vel=vel

    if(ice_flag.eq.1) then
        parcel1%npartice=npart(n_bin_modew+1:2*n_bin_modew)
        parcel1%yice(1:n_bin_modew)=mbin(n_bin_modew+1:2*n_bin_modew,n_comps+1)
        parcel1%maerice=sum(mbin(1+n_bin_modew:2*n_bin_modew,1:n_comps),2)
        parcel1%mbinice=mbin(1+n_bin_modew:2*n_bin_modew,1:n_comps+1)
        
        ! starting, so no need to set ice moments here
        
        ! other aerosol properties were set in set-up
        ! if drops create ice crystals we need to add phi and nmon
        parcel1%moments(1:n_bin_modew,n_comps+1)=npart(1:n_bin_modew)
        parcel1%moments(1:n_bin_modew,n_comps+2)=npart(1:n_bin_modew)
    endif
        
    end subroutine write_sce_to_bmm
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! integrate lognormal                                                          !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>finds the number of aerosol particles between two limits
	!>@param[in] n_aer1,d_aer1, sig_aer1: lognormal parameters
	!>@param[in] n_intern: number of internal modes
	!>@param[in] dmin,dmax: max / min diameters
	!>@param[inout] num: number of aerosol particles
    subroutine lognormal_n_between_limits(n_aer1,d_aer1,sig_aer1,n_intern,dmin,dmax, &
                                        num)
    implicit none
    real(wp), intent(in) :: dmin,dmax
    real(wp), intent(in), dimension(n_intern) :: n_aer1,d_aer1,sig_aer1
    integer(i4b), intent(in) :: n_intern
    real(wp), intent(inout) :: num
    
    integer(i4b) :: i
       
    num=0._wp                                 
    do i=1,n_intern
        num=num+n_aer1(i)*(0.5_wp*erfc(-log(dmax/d_aer1(i))/sqrt(2._wp)/sig_aer1(i) ) - &
            0.5_wp*erfc(-log(dmin/d_aer1(i))/sqrt(2._wp)/sig_aer1(i) ))
    enddo
    
    end subroutine lognormal_n_between_limits
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate the upper diameter of bin edge                                     !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>called using zbrent to find the upper bin edge in size distribution
	!>@param[in] x: dmax guess
	!>@return find_upper_bin_edge: zero when root found
    function find_upper_diameter(x)
        use numerics_type
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: find_upper_diameter, num
        
        call lognormal_n_between_limits(n_aer1(:,idum),d_aer1(:,idum),sig_aer1(:,idum), &
                                    n_intern,d_dummy,x, num)
        find_upper_diameter=num-n_dummy
        
    end function find_upper_diameter
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate the height of cloud base                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>called using zbrent to find the pressure of cloud base
	!>@param[in] p: pressure
	!>@return cloud_base: zero when root found
	function cloud_base(p)
	use numerics_type
	implicit none
	real(wp), intent(in) :: p
	real(wp) :: t,qs, cloud_base, rm, cpm
	
	rm=ra+rv*parcel1%qinit
	cpm=cp+cpv**parcel1%qinit
	
	t=theta_init*(p/1.e5_wp)**(rm/cpm)
	qs=eps1*svp_liq(t)/(p-svp_liq(t))
	
	cloud_base=qs-parcel1%q_cbase
	
	end function cloud_base
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	subroutine hydrostatic1(p,z,dzdp)
	use numerics_type
	implicit none
	real(wp), intent(in) :: p
	real(wp), dimension(:), intent(in) :: z
	real(wp), dimension(:), intent(out) :: dzdp
	real(wp) :: t
	
	t=theta_surf*(p/1.e5_wp)**(ra/cp)
	dzdp(1)=-(ra*t) / (grav*p)
	
	end subroutine hydrostatic1

	subroutine hydrostatic1a(z,p,dpdz)
	use numerics_type
	implicit none
	real(wp), intent(in) :: z
	real(wp), dimension(:), intent(in) :: p
	real(wp), dimension(:), intent(out) :: dpdz
	real(wp) :: t
	
	t=theta_surf*(p(1)/1.e5_wp)**(ra/cp)
	dpdz(1)=-(grav*p(1)) / (ra*t) 
	
	end subroutine hydrostatic1a

	subroutine hydrostatic1b(z,p,dpdz)
	use numerics_type
	implicit none
	real(wp), intent(in) :: z
	real(wp), dimension(:), intent(in) :: p
	real(wp), dimension(:), intent(out) :: dpdz
	real(wp) :: t, var, dummy, theta
	integer(i4b) :: iloc
	
	! interpolate to find theta
    iloc=find_pos(z_read(1:n_levels_s),z)
    iloc=min(n_levels_s-1,iloc)
    iloc=max(1,iloc)
    ! linear interp theta
    call poly_int(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
                min(z,z_read(n_levels_s)), var,dummy)
    theta=var     
                
	t=theta*(p(1)/1.e5_wp)**(ra/cp)
	dpdz(1)=-(grav*p(1)) / (ra*t) 
	
	end subroutine hydrostatic1b

	subroutine hydrostatic2(p,z,dzdp)
	use numerics_type
	use numerics, only : zeroin
	implicit none
	real(wp), intent(in) :: p
	real(wp), dimension(:), intent(in) :: z
	real(wp), dimension(:), intent(out) :: dzdp
	real(wp) :: t
	
	p111=p
	t=theta_surf*(p111/1.e5_wp)**(ra/cp)
	t=zeroin(t,t1old*1.01_wp,calc_theta_q,1.e-5_wp)
	
!	print *,'hi',t,calc_theta_q(t)
	! find the temperature by iteration
	dzdp(1)=-(ra*t) / (grav*p)
	
	end subroutine hydrostatic2

	subroutine hydrostatic2a(z,p,dpdz)
	use numerics_type
	use numerics, only : zeroin
	implicit none
	real(wp), intent(in) :: z
	real(wp), dimension(:), intent(in) :: p
	real(wp), dimension(:), intent(out) :: dpdz
	real(wp) :: t
	
	p111=p(1)
	t=theta_surf*(p111/1.e5_wp)**(ra/cp)
	t=zeroin(t,t1old*1.01_wp,calc_theta_q,1.e-5_wp)
!	print *,'hi',t,calc_theta_q(t)
	! find the temperature by iteration
	dpdz(1)=-(grav*p(1))/(ra*t)
	
	end subroutine hydrostatic2a
	

	function calc_theta_q(t111)
	use numerics_type
	implicit none
	real(wp), intent(in) :: t111
	real(wp) :: calc_theta_q
	real(wp) :: qs,rm,cpm
	qs=eps1*svp_liq(t111)/(p111-svp_liq(t111))
	rm=ra+rv*qs
	cpm=cp+cpv*qs
	calc_theta_q=t111*(1.e5_wp/p111)**(rm/cpm)*exp(lv*qs/cpm/t111)-theta_q_sat

	end function calc_theta_q     

	function calc_theta_q2(p)
	use numerics_type
	implicit none
	real(wp), intent(in) :: p
	real(wp) :: calc_theta_q2
	real(wp) :: ws
	ws=eps1*svp_liq(t1old)/(p-svp_liq(t1old))
	calc_theta_q2=t1old*(1e5_wp/p)**(ra/cp)*exp(lv*ws/cp/t1old)-theta_q_sat

	end function calc_theta_q2    

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! theta q                                                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the moist potential temperature
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[in] q: total water
	!>@return calc_theta_q3: moist potential temperature
	function calc_theta_q3(t,p,q)
	use numerics_type
	implicit none
	real(wp), intent(in) :: t,p,q
	real(wp) :: calc_theta_q3
	real(wp) :: qs, rm, cpm, rh
	qs=eps1*svp_liq(t)/(p-svp_liq(t))
	rm=ra+rv*q
	cpm=cp+cpv*qs !+cpw*max(q-qs,0._wp)
	rh=q/qs
	calc_theta_q3=t*(1.e5_wp/p)**(rm/cpm)*exp(lv*q/(t*(cpm))) * &
	    rh**(-q*rv/(cp+q*cpv))

	end function calc_theta_q3    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over liquid                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over liquid water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_liq: saturation vapour pressure over liquid water
	function svp_liq(t)
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp_liq
		svp_liq = 100._wp*6.1121_wp* &
			  exp((18.678_wp - (t-ttr)/ 234.5_wp)* &
			  (t-ttr)/(257.14_wp + (t-ttr)))
	end function svp_liq
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! saturation vapour pressure over ice                                          !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the saturation vapour pressure over ice water according to buck fit
	!>@param[in] t: temperature
	!>@return svp_ice: saturation vapour pressure over ice water
	function svp_ice(t)
		use numerics_type
		implicit none
		real(wp), intent(in) :: t
		real(wp) :: svp_ice
		svp_ice = 100._wp*6.1115_wp* &
            exp((23.036_wp - (t-ttr)/ 333.7_wp)* &
            (t-ttr)/(279.82_wp + (t-ttr)))

	end function svp_ice


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! diffusivity of water vapour in air										   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the diffusivity of water vapour in air
	!>@param[in] t: temperature, p: pressure
	!>@return dd: diffusivity of water vapour in air
    function dd(t,p)
      use numerics_type
      implicit none
      real(wp), intent(in) :: t, p
      real(wp) :: dd, t1
      t1=max(t,200._wp)
      dd=2.11e-5_wp*(t1/ttr)**1.94_wp*(101325._wp/p)
    end function dd

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! conductivity of water vapour												   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the thermal conductivity of air
	!>@param[in] t: temperature
	!>@return ka: thermal conductivity of air
    function ka(t)
      use numerics_type
      implicit none
      real(wp), intent(in) :: t
      real(wp) :: ka, t1
      t1=max(t,200._wp)
      ka=(5.69_wp+0.017_wp*(t1-ttr))*1.e-3_wp*joules_in_a_cal
    end function ka

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! viscosity of air - page 417 pruppacher and klett							   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the viscosity of air
	!>@param[in] t: temperature
	!>@return viscosity_air: viscosity of air
    function viscosity_air(t)
        use numerics_type
        implicit none
        real(wp), intent(in) :: t
        real(wp) :: viscosity_air
        real(wp) :: tc

        tc = t-ttr
        tc = max(tc,-200._wp)

        if( tc.ge.0._wp) then
            viscosity_air = (1.718_wp+0.0049_wp*tc) * 1e-5_wp ! the 1d-5 converts from poise to si units
        else
            viscosity_air = (1.718_wp+0.0049_wp*tc-1.2e-5_wp*tc**2) * 1e-5_wp
        end if
    end function viscosity_air
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! surface tension of water - pruppacher and klett							   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the surface tension of water
	!>@param[in] t: temperature
	!>@return surface_tension: surface_tension of water
    function surface_tension(t)	
        use numerics_type
        real(wp), intent(in) :: t
        real(wp) :: tc, surface_tension

        tc=t-ttr
        tc = max(tc,-40._wp)

        ! pruppacher and klett pg 130 
        surface_tension = 75.93_wp + 0.115_wp * tc + 6.818e-2_wp * tc**2 + &
                          6.511e-3_wp * tc**3 + 2.933e-4_wp * tc**4 + &
                          6.283e-6_wp * tc**5 + 5.285e-8_wp * tc**6
        if(tc.ge.0._wp) then
            surface_tension = 76.1_wp - 0.155_wp*tc
        end if
    
        surface_tension = surface_tension*joules_in_an_erg ! convert to j/cm2 
        surface_tension = surface_tension*1.e4_wp ! convert to j/m2 

    !    surface_tension=72d-3
        !sigma = 75.93_wp * joules_in_an_erg*1d4
    end function surface_tension
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the wet diameter                                				   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the wet diameter of a particle
	!>@param[in] t: temperature
	!>@param[in] mwat: mass of water
	!>@param[in] mbin: mass of aerosol components in each bin
	!>@param[in] rhobin: density of each component
	!>@param[in] sz: length of array
	!>@param[inout] dw: wet diameter
    subroutine wetdiam(mwat,mbin,rhobin,sz,dw) 
      use numerics_type
      implicit none
      real(wp), dimension(:), intent(in) :: mwat
      real(wp), dimension(:,:), intent(in) :: mbin,rhobin
      integer(i4b), intent(in) :: sz
      real(wp), dimension(:),intent(inout) :: dw
      
      real(wp), dimension(sz) :: rhoat

      ! calculate the diameter and radius
      rhoat(:)=mwat(:)/rhow+sum(mbin(:,1:n_comps)/rhobin(:,:),2)
      rhoat(:)=(mwat(:)+sum(mbin(:,1:n_comps),2))/rhoat(:);
  
      ! wet diameter:
      dw(:)=((mwat(:)+sum(mbin(:,1:n_comps),2))*6._wp/(pi*rhoat(:)))**(1._wp/3._wp)
      
    end subroutine wetdiam
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the equilibrium humidity over a particle        				   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the equilibrium humidity over a particle
	!>@param[in] t: temperature
	!>@param[in] mwat: mass of water
	!>@param[in] mbin: mass of aerosol components in each bin
	!>@param[in] rhobin: density of each component
	!>@param[in] nubin: van hoff factor in each bin
	!>@param[in] molwbin: molecular weight in each bin
	!>@param[in] sz: length of array
	!>@param[inout] rh_eq: equilibrium humidity
	!>@param[inout] rhoat: density of particle
	!>@param[inout] dw: wet diameter
    subroutine koehler01(t,mwat,mbin,rhobin,nubin,molwbin,sz,rh_eq,rhoat,dw) 
      use numerics_type
      implicit none
      real(wp), dimension(:), intent(in) :: mwat
      real(wp), dimension(:,:), intent(in) :: mbin,rhobin,nubin,molwbin
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz) :: nw
      real(wp), dimension(:),intent(inout) :: rh_eq,rhoat, dw
      real(wp), intent(in) :: t
      real(wp) :: sigma

      ! calculate the diameter and radius
      nw(:)=mwat(:)/molw_water
      rhoat(:)=mwat(:)/rhow+sum(mbin(:,1:n_comps)/rhobin(:,:),2)
      rhoat(:)=(mwat(:)+sum(mbin(:,1:n_comps),2))/rhoat(:);
  
      ! wet diameter:
      dw(:)=((mwat(:)+sum(mbin(:,1:n_comps),2))*6._wp/(pi*rhoat(:)))**(1._wp/3._wp)
  
      ! calculate surface tension
      sigma=surface_tension(t)

      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      rh_eq(:)=exp(4._wp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))* &
           (nw(:))/(nw(:)+sum(mbin(:,1:n_comps)/molwbin(:,:)*nubin(:,:),2) ) 

    end subroutine koehler01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the equilibrium humidity over a particle, using k-koehler          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the equilibrium humidity over a particle using K-koehler theory
	!>@param[in] t: temperature
	!>@param[in] mwat: mass of water
	!>@param[in] mbin: mass of aerosol components in each bin
	!>@param[in] rhobin: density of each component
	!>@param[in] kappabin: kappa in each bin
	!>@param[in] molwbin: molecular weight in each bin
	!>@param[in] sz: length of array
	!>@param[inout] rh_eq: equilibrium humidity
	!>@param[inout] rhoat: density of particle
	!>@param[inout] dw: wet diameter
    subroutine kkoehler01(t,mwat,mbin,rhobin,kappabin,molwbin,sz,rh_eq,rhoat,dw) 
      use numerics_type
      implicit none
      real(wp), dimension(:), intent(in) :: mwat
      real(wp), dimension(:,:), intent(in) :: mbin,rhobin,kappabin,molwbin
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz) :: nw,dd,kappa
      real(wp), dimension(:),intent(inout) :: rh_eq,rhoat, dw
      real(wp), intent(in) :: t
      real(wp) :: sigma
      ! calculate the diameter and radius
      nw(:)=mwat(:)/molw_water
      rhoat(:)=mwat(:)/rhow+sum(mbin(:,1:n_comps)/rhobin(:,:),2)
      rhoat(:)=(mwat(:)+sum(mbin(:,1:n_comps),2))/rhoat(:);
  
      ! wet diameter:
      dw(:)=((mwat(:)+sum(mbin(:,1:n_comps),2))* 6._wp/(pi*rhoat(:)))**(1._wp/3._wp)
  
      dd(:)=((sum(mbin(:,1:n_comps)/rhobin(:,:),2))*6._wp/(pi))**(1._wp/3._wp) ! dry diameter
                                  ! needed for eqn 6, petters and kreidenweis (2007)
  
      kappa(:)=sum(mbin(:,1:n_comps)/rhobin(:,:)*kappabin(:,:),2) &
               / sum(mbin(:,1:n_comps)/rhobin(:,:),2)
               ! equation 7, petters and kreidenweis (2007)

      ! calculate surface tension
      sigma=surface_tension(t)

      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      rh_eq(:)=exp(4._wp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))* &
           (dw(:)**3-dd(:)**3)/(dw(:)**3-dd(:)**3*(1._wp-kappa))
           ! eq 6 petters and kreidenweis (acp, 2007)
    end subroutine kkoehler01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! koehler equations                                                     	   !
    ! this is coded so it can be called with a root-finder, to find the inverse	   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the equilibrium humidity over a particle using koehler theory
	!>@param[in] nw: number of moles of water
	!>@return koehler02: equilibrium rh - but called via root-finder, so nw is returned
    function koehler02(nw) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: nw
      real(wp) :: massw
      real(wp) :: rhoat, dw,koehler02
      real(wp) :: sigma
      
      

      ! wet diameter:
      massw=nw*molw_water
      rhoat=massw/rhow+sum(parcel1%mbin(n_sel,1:n_comps) / &
            parcel1%rhobin(n_sel,1:n_comps))
      rhoat=(massw+parcel1%maer(n_sel))/rhoat;
      dw=((massw+parcel1%maer(n_sel))* 6._wp/(pi*rhoat))**(1._wp/3._wp)
  
      ! calculate surface tension
      sigma=surface_tension(parcel1%t)
  
      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      koehler02=mult*(exp(4._wp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           (nw)/(nw+sum(parcel1%mbin(n_sel,1:n_comps)/ &
           parcel1%molwbin(n_sel,1:n_comps) * &
           parcel1%nubin(n_sel,1:n_comps)) ))-rh_act

    end function koehler02
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! kappa-koehler equations                                                      !
    ! this is coded so it can be called with a root-finder, to find the inverse	   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the equilibrium humidity over a particle using K-koehler theory
	!>@param[in] nw: number of moles of water
	!>@return kkoehler02: equilibrium rh - but called via root-finder, so nw is returned
    function kkoehler02(nw) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: nw
      real(wp) :: massw
      real(wp) :: rhoat, dw,dd,kappa,kkoehler02
      real(wp) :: sigma

      ! wet diameter:
      massw=nw*molw_water
      rhoat=massw/rhow+sum(parcel1%mbin(n_sel,1:n_comps) / &
            parcel1%rhobin(n_sel,1:n_comps))
      rhoat=(massw+parcel1%maer(n_sel))/rhoat;
      dw=((massw+parcel1%maer(n_sel))* 6._wp/(pi*rhoat))**(1._wp/3._wp)
  
      ! calculate surface tension
      sigma=surface_tension(parcel1%t)
  
      dd=(sum(parcel1%mbin(n_sel,1:n_comps) / parcel1%rhobin(n_sel,:))* &
          6._wp/(pi))**(1._wp/3._wp) ! dry diameter
                                  ! needed for eqn 6, petters and kreidenweis (2007)
  
      kappa=sum(parcel1%mbin(n_sel,1:n_comps) / parcel1%rhobin(n_sel,1:n_comps)* &
               parcel1%kappabin(n_sel,:)) &
               / sum(parcel1%mbin(n_sel,1:n_comps) / parcel1%rhobin(n_sel,1:n_comps))
               ! equation 7, petters and kreidenweis (2007)

      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      kkoehler02=mult*(exp(4._wp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           (dw**3-dd**3)/(dw**3-dd**3*(1._wp-kappa)))-rh_act
           ! eq 6 petters and kreidenweis (acp, 2007)
    end function kkoehler02
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! koehler equations                                                     	   !
    ! this is coded so it can be called with a root-finder, to find the inverse	   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the equilibrium humidity over a particle using koehler theory
	!> n_sel must be set to the mode, d_dummy the water in the bin to return the 
	!> aerosol dry mass
	!>@param[in] nw: number of moles of water
	!>@return koehler03: equilibrium rh - but called via root-finder, so mbin is returned
    function koehler03(mbin) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: mbin ! mass of aerosol particle
      real(wp) :: massw
      real(wp) :: rhoat, dw,nw,koehler03
      real(wp) :: sigma
      
      ! calculate the diameter and radius
      nw=d_dummy/molw_water ! moles of water
      rhoat=d_dummy/rhow+mbin* sum(mass_frac_aer1(n_sel,1:n_comps)/ &
                           density_core1(1:n_comps))
      rhoat=(d_dummy+mbin)/rhoat;
      dw=((d_dummy+mbin)* 6._wp/(pi*rhoat))**(1._wp/3._wp)
  
      ! calculate surface tension
      sigma=surface_tension(parcel1%t)
  
      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      koehler03=mult*(exp(4._wp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           (nw)/(nw+mbin* &
           sum(mass_frac_aer1(n_sel,1:n_comps)/ &
           molw_core1(1:n_comps)* &
           nu_core1(1:n_comps)) ))-rh_act


    end function koehler03
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! K-koehler equations                                                     	   !
    ! this is coded so it can be called with a root-finder, to find the inverse	   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the equilibrium humidity over a particle using koehler theory
	!> n_sel must be set to the mode, d_dummy the water in the bin to return the 
	!> aerosol dry mass
	!>@param[in] nw: number of moles of water
	!>@return kkoehler03: equilibrium rh - but called via root-finder, so mbin is returned
    function kkoehler03(mbin) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: mbin ! mass of aerosol particle
      real(wp) :: massw
      real(wp) :: rhoat, dw,nw,dd,kappa, kkoehler03
      real(wp) :: sigma
      ! calculate the diameter and radius
      nw=d_dummy/molw_water ! moles of water
      rhoat=d_dummy/rhow+mbin* sum(mass_frac_aer1(n_sel,1:n_comps)/ &
                           density_core1(1:n_comps))
      rhoat=(d_dummy+mbin)/rhoat;
      dw=((d_dummy+mbin)* 6._wp/(pi*rhoat))**(1._wp/3._wp)
  
      dd=((sum(mbin*mass_frac_aer1(n_sel,1:n_comps)/ density_core1(1:n_comps),1))* &
          6._wp/(pi))**(1._wp/3._wp) ! dry diameter
                                  ! needed for eqn 6, petters and kreidenweis (2007)
  
      kappa=sum(mbin*mass_frac_aer1(n_sel,1:n_comps) &
               / density_core1(1:n_comps)* &
               kappa_core1(1:n_comps),1) &
               / sum(mbin*mass_frac_aer1(n_sel,1:n_comps) &
               /density_core1(1:n_comps),1)
               ! equation 7, petters and kreidenweis (2007)
           
      ! calculate surface tension
      sigma=surface_tension(parcel1%t)
  
      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      kkoehler03=mult*(exp(4._wp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           (dw**3-dd**3)/(dw**3-dd**3*(1._wp-kappa)))-rh_act
           ! eq 6 petters and kreidenweis (acp, 2007)
    end function kkoehler03
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the terminal velocity of cloud drops                               !
    ! see pruppacher and klett                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the terminal velocity of water drops
	!>@param[inout] vel, nre, cd: terminal velocity, reynolds number and drag coefficient
	!>@param[in] diam, rhoat, t, p: diameter, density, temperature and pressure
	!>@param[in] sz: size of the array to calculate terminal velocities
    subroutine terminal01(vel,diam,rhoat, t,p,nre,cd,sz)
      use numerics_type
      implicit none
      real(wp), intent(in) :: t, p
      real(wp), dimension(:), intent(in) :: diam, rhoat
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz), intent(inout) :: nre,cd,vel 
      real(wp) :: eta, sigma, physnum, phys6, mfpath, tc, rhoa
      real(wp), dimension(sz) :: diam2, mass, bondnum, bestnm, x,y


        tc=t-ttr
        vel = 0._wp ! zero array
        rhoa = p / (ra * t) ! density of air
        diam2=diam ! temporary array that can be changed

        eta = viscosity_air(t)

        nre = 0._wp ! zero array
        where(diam2.gt.7000.e-6_wp)
            diam2=7000.e-6_wp
        end where
        mass=pi/6._wp*diam2**3*rhoat
    
        sigma = surface_tension(t)
        
        ! regime 3:  eqns 5-12, 10-146 & 10-148 from p & k 
        physnum = (sigma**3._wp) * (rhoa**2_wp) / ((eta**4._wp) * grav * (rhow - rhoa))		
        phys6 = physnum**(1._wp / 6._wp)
        where(diam2.gt.1070.e-6_wp) 
            bondnum = (4._wp/3._wp)*grav * (rhow - rhoa) * (diam2**2) / sigma

            x = log(bondnum*phys6)
            y = -5.00015_wp + 5.23778_wp * x - 2.04914_wp * x * x + 0.475294_wp * (x**3) &
                - 0.542819e-1_wp * (x**4._wp) + 0.238449e-2_wp * (x**5)

            nre = phys6 * exp(y)

            vel = eta * (nre)/ (rhoa * diam2)

            cd = 8._wp * mass * grav * rhoa/(pi * ((diam2 / 2._wp)* eta)**2)
            cd = cd	/ (nre**2) 
        end where

        ! regime 2:  eqns 10-142, 10-145 & 10-146 from p & k 
        where(diam2.le.1070.e-6_wp.and.diam2.gt.20.e-6_wp)
            bestnm = 32._wp * ((diam2 / 2._wp)**3) * (rhow - rhoa) * rhoa * &
                      grav / (3._wp * eta**2)
            x = log(bestnm)
            y = -3.18657_wp + 0.992696_wp * x - 0.153193e-2_wp * x * x &
                -0.987059e-3_wp * (x**3) - 0.578878e-3_wp * (x**4) &
                + 0.855176e-4_wp * (x**5) - 0.327815e-5_wp * (x**6)
            nre =  exp(y)
            vel = eta * nre / (2._wp * rhoa * (diam2 / 2._wp))
            cd = bestnm/(nre**2)
        end where

        ! regime 1:  eqns 10-138, 10-139 & 10-140 from p & k 
        mfpath = 6.6e-8_wp * (101325_wp / p) * (t / 293.15_wp)
        where(diam2.le.20.e-6_wp) 
            vel = 2._wp * ((diam2 / 2._wp)**2) * grav * (rhow - rhoa) / (9._wp * eta)
            vel = vel * (1._wp + 1.26_wp * mfpath / (diam2 / 2._wp))
            nre = vel * rhoa * diam2 / eta

            cd = 8._wp * mass * grav * rhoa/(pi * ((diam2 / 2._wp)* eta)**2)
            cd = cd	/ (nre**2) 
         end where

        where(isnan(vel))
          vel=0._wp
        end where
    end subroutine terminal01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the ventilation coefficient for cloud drops                        !
    ! see  pruppacher and klett - page 538-541                                     !
    ! original reference: pruppacher and rasmussen 1979                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the ventilation factor for water drops
	!>@param[inout] fv, fh: ventilation factors for vapour and heat
	!>@param[in] diam, rhoat, t, p: diameter, density, temperature and pressure
	!>@param[in] sz: size of the array to calculate ventilation factors
    subroutine ventilation01(diam, rhoat,t, p, fv, fh,sz)
      use numerics_type
      implicit none
      real(wp), intent(in) :: t, p
      real(wp), dimension(:), intent(in) :: diam, rhoat
      real(wp), dimension(sz), intent(inout) :: fv,fh
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz) :: nre,cd,vel,calc
      real(wp) :: d1,k1,rhoa, eta, nu, nsc1,nsc2
      ! density of air
      rhoa = p/ra/t
      ! diffusivity of water vapour in air
      d1 = dd(t,p)
      ! conductivity of air
      k1 = ka(t)
      ! viscosity of air
      eta=viscosity_air(t)
      ! kinematic viscosity
      nu = eta / rhoa
      ! schmitt numbers:
      nsc1 = nu / d1
      nsc2 = nu / k1

      ! terminal velocity of water drops
      call terminal01(vel,diam,rhoat, t,p,nre,cd,sz)
  
      ! mass ventilation - use dv+++++++++
      calc = (nsc1**(1._wp/3._wp)) * sqrt(nre)
      where(calc.gt.51.4_wp)
        calc=51.4_wp
      end where

      where(calc.lt.1.4_wp)
        fv=1.00_wp+0.108_wp*calc**2
      elsewhere
        fv=0.78_wp+0.308_wp*calc
      end where
      !-----------------------------------
    
      ! heat ventilation - use ka---------
      calc = (nsc2**(1._wp/3._wp)) * sqrt(nre)
      where(calc.gt.51.4_wp)
        calc=51.4_wp
      end where

      where(calc.lt.1.4_wp)
        fh=1.00_wp+0.108_wp*calc**2
      elsewhere
        fh=0.78_wp+0.308_wp*calc
      end where
      !-----------------------------------
    end subroutine ventilation01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the terminal velocity of ice particles                             !
    ! see heymsfield and westbrook (2010, jas)                                     !
    ! who corrected a bias in the fall speeds derived by mitchell and others       !
    ! this method also works reasonably well for sub-100 micron crystals           !
    ! see westbrook qj paper for latter point                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the terminal velocity of ice crystals
	!>@param[inout] vel, nre: terminal velocity and reynolds number
	!>@param[in] mwat, t, p, phi, rhoi, nump, rime
	!>@param[in] sz: size of the array to calculate terminal velocities
    subroutine terminal02(vel,mwat, t,p,phi,rhoi,nump,rime,nre,sz)
      use numerics_type
      implicit none
      real(wp), intent(in) :: t, p
      real(wp), dimension(:), intent(in) :: mwat, phi, rhoi, nump, rime
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz), intent(inout) :: nre,vel 
      real(wp) :: eta, rhoa
      real(wp), dimension(sz) :: dmax,drime,area,ar, x
      integer(i4b) :: i	
      ! air properties
      eta = viscosity_air(t)
      rhoa = p/ra/t
    
      ! calculate the maximum dimension of the particle
!       call maxdimension01(mwat-rime,rhoi,phi,nump,rime,dmax,drime,sz)
!   
!       ! calculate the area of the particle
!       call areaaggregates01(area,mwat-rime,rhoi,phi,nump,dmax,drime,sz)
!       ! area ratio
!       ar=area/(pi/4._wp* (dmax**2._wp))
!       ar=min(max(ar,0.1_wp),1._wp)
      dmax=(mwat/rhoice*6._wp/pi)**(1/3)
      ar=1._wp
  
      ! heymsfield and westbrook
      x=rhoa*8._wp*mwat*grav/( (eta**2._wp)*pi*(ar**0.5_wp))
      nre=(8.0_wp**2._wp)/4._wp* &
          ( (sqrt(1._wp+(4._wp*sqrt(x))/( (8._wp**2._wp)*sqrt(0.35_wp)))-1._wp)**2._wp)
      vel=eta*(nre)/(rhoa*dmax)
    
      ! viscous regime
      where(nre.lt.1._wp) 
        vel = grav*mwat / (6._wp*pi*eta*0.465_wp*dmax*(ar**0.5_wp))
      end where

      where(isnan(vel)) 
            vel=0._wp
      end where
    end subroutine terminal02
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate ventilation coefficients for ice crystals                          !
    ! see page 553 p+k                                                             !
    ! original reference: ji 1991 and wang and ji 1992                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the ventilation coefficients for ice crystals
	!>@param[inout] fv, fh: ventilation coefficients for mass and heat
	!>@param[in] mwat, t, p, phi, rhoi, nump, rime
	!>@param[in] sz: size of the array to calculate terminal velocities
    subroutine ventilation02(mwat, t, p,phi, rhoi, nump,rime,fv, fh,sz)
      use numerics_type

      implicit none
      real(wp), intent(in) :: t, p
      real(wp), dimension(:), intent(in) :: mwat, phi, rhoi, nump, rime
      real(wp), dimension(sz), intent(inout) :: fv,fh
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz) :: nre,vel,nre2,x
      real(wp) :: d1,k1,rhoa, eta, nu, nsc1,nsc2, calc1, calc2
      ! density of air
      rhoa = p/ra/t
      ! diffusivity of water vapour in air
      d1 = dd(t,p)
      ! conductivity of air
      k1 = ka(t)
      ! viscosity of air
      eta=viscosity_air(t)
      ! kinematic viscosity
      nu = eta / rhoa
      ! schmitt numbers:
      nsc1 = nu / d1
      nsc2 = nu / k1

      ! terminal velocity of ice crystals
      call terminal02(vel,mwat, t,p,phi,rhoi,nump,rime,nre,sz)

      ! mass ventilation - use dv; heat ventilation - use ka +++++++
      calc1 = nsc1**(1._wp/3._wp)
      calc2 = nsc2**(1._wp/3._wp)
  
      ! columns
      nre2=min(nre,20._wp)
      where(phi.gt.1.0_wp)  
        x = calc1*sqrt(nre2)	
        fv = 1.0_wp - 0.000668_wp*x/4._wp + 2.39402_wp*((x/4._wp)**2._wp) + &
             0.73409_wp*((x/4._wp)**3._wp)-0.73911_wp*((x/4._wp)**4._wp)
        x = calc2*sqrt(nre2);	
        fh = 1.0_wp - 0.000668_wp*x/4._wp + 2.39402_wp*((x/4._wp)**2._wp) + &
             0.73409_wp*((x/4._wp)**3._wp)-0.73911_wp*((x/4._wp)**4._wp)
      end where
      !--------
  
      ! plates
      nre2=min(nre,120._wp)
      where(phi.le.1._wp) 
        x = calc1*sqrt(nre2)	
        fv = 1.0_wp - 0.06042_wp*x/10._wp + 2.79820_wp*((x/10._wp)**2._wp) - &
             0.31933_wp*((x/10._wp)**3._wp)-0.06247_wp*((x/10._wp)**4._wp)
        x = calc2*sqrt(nre2)	
        fh = 1.0_wp - 0.06042_wp*x/10._wp + 2.79820_wp*((x/10._wp)**2._wp) - &
             0.31933_wp*((x/10._wp)**3._wp)-0.06247_wp*((x/10._wp)**4._wp)
      end where
      !-------
  
      ! broad-branched crystals
      !nre2=min(nre,120d0) ! already done above
      where(phi.lt.0.2_wp.and.rhoi.le.500._wp) 
        x = calc1*sqrt(nre2)	
        fv = 1.0_wp + 0.35463_wp*x/10._wp + 3.55333_wp*((x/10._wp)**2._wp)
        x = calc2*sqrt(nre2)
        fh = 1.0_wp + 0.35463_wp*x/10._wp + 3.55333_wp*((x/10._wp)**2._wp)
      end where
      ! -----------------------
    end subroutine ventilation02
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate growth rate of a cloud droplet					  				   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates diffusional growth rate for water drops
	!>@param[in] t,p, rh: temperature, pressure, and rh
	!>@param[in] rh_eq, rhoat, diam: equilibrium rh, density, diameter
	!>@param[in] sz: size of array
	!>@return dropgrowthrate01: growth rate of drops
    function dropgrowthrate01(t,p,rh,rh_eq,rhoat,diam,sz) 
      use numerics_type
      implicit none
      real(wp), intent(in) :: t, p, rh
      real(wp), dimension(:), intent(in) :: rh_eq,rhoat, diam
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz) :: dropgrowthrate01
      real(wp), dimension(sz) :: rad, dstar,kstar,fv,fh
      real(wp) :: d1,k1,rhoa
  
      rad=diam/2._wp
      ! density of air
      rhoa=p/ra/t
      ! diffusivity of water vapour in air
      d1=dd(t,p)
      ! thermal conductivity of air
      k1=ka(t)
      ! ventilation coefficient
      fv=1._wp
      fh=1._wp
      if(vent_flag.eq.1) then
        call ventilation01(diam, rhoat,t, p, fv, fh,sz)
      end if
  
      ! modify diffusivity and conductivity
      dstar=d1*fv/(rad/(rad+0.7_wp*8.e-8_wp)+d1*fv/rad/alpha_cond*sqrt(2._wp*pi/rv/t))
      kstar=k1*fh/(rad/(rad+2.16e-7_wp)+k1*fh/rad/alpha_therm/cp/rhoa*sqrt(2._wp*pi/ra/t))
  
      ! 455 jacobson and 511 pruppacher and klett
      dropgrowthrate01=dstar*lv*rh_eq*svp_liq(t)* &
                       rhoat/kstar/t*(lv*molw_water/t/r_gas-1._wp) 
      dropgrowthrate01=dropgrowthrate01+rhoat*r_gas*t/molw_water  
      dropgrowthrate01=dstar*(rh-rh_eq)*svp_liq(t)/rad/dropgrowthrate01
                 
    end function dropgrowthrate01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate growth rate of an ice crystal					  				   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function icegrowthrate01(t,p,rh_ice,rh_eq,mwat,mbin,rhobin,phi, rhoi, nump,rime,sz) 
      use numerics_type

      implicit none
      real(wp), intent(in) :: t, p, rh_ice
      real(wp), dimension(:), intent(in) :: rh_eq,mwat, phi, rhoi, nump, rime
      real(wp), dimension(:,:), intent(in) :: mbin,rhobin
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz) :: icegrowthrate01
      real(wp), dimension(sz) :: rad, dstar,kstar,rhoat,diam,fv,fh
      real(wp) :: d1,k1,rhoa
      integer(i4b) :: i

      ! calculate the capacitance - get rid of yiceold as messy
      rad=capacitance01(mwat(:),phi,rhoi,nump,rime,sz)
      ! density of air
      rhoa=p/ra/t
      ! diffusivity of water vapour in air
      d1=dd(t,p)
      ! thermal conductivity of air
      k1=ka(t)
      ! ventilation coefficient
      fv=1._wp
      fh=1._wp
      if(vent_flag.eq.1) then
        call ventilation02(mwat, t, p,phi, &
             rhoi,nump, rime,fv, fh,sz)
      end if
      ! modify diffusivity and conductivity
      dstar=d1*fv/(rad/(rad+0.7_wp*8e-8_wp)+d1*fv/rad/alpha_dep*sqrt(2._wp*pi/rv/t))
      kstar=k1*fh/(rad/(rad+2.16e-7_wp)+k1*fh/rad/alpha_therm_ice/cp/rhoa*sqrt(2._wp*pi/ra/t))
  
      ! 473 jacobson 
      icegrowthrate01=dstar*ls*rh_eq*svp_ice(t)/ &
                       kstar/t*(ls*molw_water/t/r_gas-1._wp) 
      icegrowthrate01=icegrowthrate01+r_gas*t/molw_water  
      icegrowthrate01=4._wp*pi*rad*dstar*(rh_ice-rh_eq)*svp_ice(t)/icegrowthrate01
                 
    end function icegrowthrate01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the capacitance of an ice crystal				  				   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculates the capacitance of oblate (longer a) and prolate (longer c) 
    ! spheroids. a is the length of the prism axis (half of) and c the basal 
    ! (half of). see page 1214 of chen and lamb, jas, 1994.
    function capacitance01(mwat,phi,rhoice,numi,rimemass,sz)
      use numerics_type
      implicit none
      real(wp), dimension(:), intent(in) :: mwat, phi, rhoice, numi, rimemass
      real(wp), dimension(sz) :: capacitance01,vol,a,c,ecc,dmax,drime
  
      integer(i4b), intent(in) :: sz
  
      vol=mwat/rhoice
  
      a=( 3._wp*vol/(4._wp*pi*phi) )**(1._wp/3._wp)
      c=a*phi
  
      where(phi.lt.1._wp)
        ecc=sqrt(1._wp-phi**2._wp)
        capacitance01=a*ecc/asin(ecc)
      elsewhere
        ecc=sqrt(1._wp-phi**(-2._wp))
        capacitance01=c*ecc/log((1._wp+ecc)*phi)
      end where
    
      where(abs(phi-1._wp).lt.1.e-4_wp)
        capacitance01=a
      end where
      ! westbrook et al. (2008, jas): capacitance of aggregates is 0.25 times the 
      ! maximum dimension of the ice particle
!       call maxdimension01(mwat-rimemass,rhoice,phi,numi,rimemass,dmax,drime,sz)
!       where(numi.ge.2._wp)
!          capacitance01=0.25_wp*dmax
!       end where
  
    end function capacitance01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! koop et al 2000 nucleation rate                                              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates nucleation rate of ice in supercooled water according to koop etal (2000)
	!>@param[in] aw, t, p: activity of water, temperature, and pressure
	!>@param[in] sz: size of array
	!>@return koopnucrate: nucleation rate
    function koopnucrate(aw,t,p,sz)
          use numerics_type
          implicit none
          integer(i4b), intent(in) :: sz
          real(wp), parameter :: k_water_amb=1.6_wp, &
                    dk_water_dp=-8.8_wp,k_ice_amb=0.22_wp,dk_ice_dp=-0.17_wp
          real(wp) :: pg,t,p, & 
                      integral3, mudiff0
          real(wp), intent(in), dimension(:) :: aw
          real(wp), dimension(sz) :: koopnucrate,deltaaw,logj
      
          pg=p/1.e9_wp

      
          integral3=(-230.76_wp - 0.1478_wp * t + 4099.2_wp * t**(-1) + &
                     48.8341_wp * log(t) ) * &
            (pg - 0.5_wp * (k_water_amb + dk_water_dp * pg)* pg**2 &
             - (1._wp/6._wp) * dk_water_dp * pg**3 ) &
             - (19.43_wp - 2.2e-3_wp * t + 1.08e-5_wp * t**2 ) * &
            (pg - 0.5_wp * (k_ice_amb + dk_ice_dp * pg) * pg**2 - &
            (1._wp/6._wp) * dk_ice_dp * pg**3 )

          mudiff0 = 210368._wp + 131.438_wp * t - 3.32373e6_wp * t**(-1)  &
                   - 41729.1_wp * log(t)
    
          ! delta activity
          deltaaw = aw * exp(integral3 / r_gas / t) - exp(mudiff0 / r_gas / t)
        

          ! nucleation rate
          logj = -906.7_wp + 8502._wp * deltaaw - 26924._wp * deltaaw**2 + 29180._wp * &
                deltaaw**3

          koopnucrate = (10._wp**logj) * 1.e6_wp;	! nucleation rate in m^-3 s^-1 
    end function koopnucrate
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! derivatives for a warm parcel model                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates rates of change for a warm parcel model
	!>@param[in] neq: length of solution vector
	!>@param[in] tt: time
	!>@param[in] y: solution vector
	!>@param[inout] ydot: derivates calculated
	!>@param[in] rpar: real data coming in
	!>@param[in] ipar: integer data coming in
    subroutine fparcelwarm(neq, tt, y, ydot, rpar, ipar)
        use numerics_type
        use numerics, only : dfsid1,find_pos

        implicit none
        real(wp), intent(inout) :: tt
        real(wp), intent(inout), dimension(neq) :: y, ydot
        integer(i4b), intent(inout) :: neq
        real(wp), intent(inout) :: rpar
        integer(i4b), intent(inout) :: ipar

        ! local variables
        real(wp) :: wv=0._wp, wl=0._wp, wi=0._wp, rm, cpm, &
                  drv=0._wp, dri=0._wp,dri2=0._wp, &
                  rh,t,p,err,sl, w, &
                  te, qve, pe, var, dummy, rhoe, rhop, b

        integer(i4b) :: i, j,iloc, ipart, ipr, ite, irh, iz,iw

        
        ipart=parcel1%n_bin_modew
        ipr=parcel1%ipr
        ite=parcel1%ite
        irh=parcel1%irh
        iz =parcel1%iz
        iw =parcel1%iw

        if ((updraft_type==3).and.(tt>t_thresh)) then
            y(iw)=winit2*cos(2._wp*pi*(tt-t_thresh)/tau2)
        endif

        rh=y(irh)
        t=y(ite)
        p=y(ipr)
        w=y(iw)
    

        ! check there are no negative values
        where(y(1:ipart).le.0.e1_wp)
            y(1:ipart)=1.e-22_wp
        end where


        ! calculate mixing ratios from rh, etc
        sl=svp_liq(t)*rh/(p-svp_liq(t)) ! saturation ratio
        sl=(sl*p/(1._wp+sl))/svp_liq(t)
        wv=eps1*rh*svp_liq(t) / (p-svp_liq(t)) ! vapour mixing ratio
        wl=sum(parcel1%npart*y(1:ipart))             ! liquid mixing ratio

        ! calculate the moist gas constants and specific heats
        rm=ra+wv*rv
        cpm=cp+wv*cpv+wl*cpw+wi*cpi

        ! now calculate derivatives
        ! adiabatic parcel model
        ydot(iz )=w                         ! vertical wind
        ydot(ipr)=-p/rm/t*grav*ydot(iz)      ! hydrostatic equation

        ! calculate equilibrium rhs
        select case (kappa_flag)
            case (0)
              call koehler01(t,y(1:ipart),parcel1%mbin(:,1:n_comps), &
                   parcel1%rhobin(:,1:n_comps), parcel1%nubin(:,1:n_comps), &
                   parcel1%molwbin(:,1:n_comps),ipart, &
                   parcel1%rh_eq,parcel1%rhoat, parcel1%dw) 
            case (1)
              call kkoehler01(t,y(1:ipart),parcel1%mbin(:,1:n_comps), &
                   parcel1%rhobin(:,1:n_comps), parcel1%kappabin(:,1:n_comps), &
                   parcel1%molwbin(:,1:n_comps),ipart, &
                   parcel1%rh_eq,parcel1%rhoat, parcel1%dw)
        case default
            print *,'error kappa_flag'
        end select
        
        
        
        

        ! particle growth rate - radius growth rate
        parcel1%da_dt=dropgrowthrate01(t,p,sl,parcel1%rh_eq, &
            parcel1%rhoat,parcel1%dw,ipart)
        ! do not bother if number concentration too small
        do i=1,ipart
            if(isnan(parcel1%da_dt(i)).or.parcel1%npart(i).le. 1.e-9_wp) then
              parcel1%da_dt(i)=0._wp
            endif
        enddo


        
        ! mass growth rate
        ydot(1:ipart)=pi*parcel1%rhoat*parcel1%dw**2 * parcel1%da_dt
        ! change in vapour content
        drv = -sum(ydot(1:ipart)*parcel1%npart)
        if((.not. adiabatic_prof) .and. (.not. vert_ent)) then ! entraining?
            !calculate the environmental p, qv, te, density
            ! parcel p, density
            ! buoyancy...
            ! locate position
            iloc=find_pos(parcel1%z_sound(1:n_levels_s),y(iz))
            iloc=min(n_levels_s-1,iloc)
            iloc=max(1,iloc)
            ! linear interp p
            call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%p_sound(iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            pe=var
            ! linear interp qv
            call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%q_sound(1,iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            qve=var
            ! linear interp te
            call poly_int(parcel1%z_sound(iloc:iloc+1), parcel1%t_sound(iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            te=var
            ! env density:
            rhoe=pe/(rm*te)
            ! parcel density:
            rhop=p/(rm*t)
            !buoyancy
            if((parcel1%z_sound(n_levels_s) .lt. y(iz)) .or. &
                (parcel1%z_sound(1) .gt. y(iz))) then
                b=0._wp
            else
                b=grav*(rhoe-rhop)/rhoe
            endif
            ! forcing
            drv=drv+w*ent_rate*(qve-wv)

        endif
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in temperature of parcel                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ydot(ite)=rm/p*ydot(ipr)*t/cpm  ! temperature change: expansion
        ydot(ite)=ydot(ite)-lv/cpm*drv ! temp change: condensation
        if((.not. adiabatic_prof) .and. (.not. vert_ent)) then ! entraining?
            ydot(ite)=ydot(ite)+w*ent_rate*(te-y(ite) + lv/cpm*(qve-wv))
            !ydot(iw) = b -w*ent_rate*y(iw)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         if(y(iz) .gt. parcel1%z_cbase) then
!             parcel1%theta_q=calc_theta_q3(t,p,wv+wl)
!             parcel1%x_ent=(parcel1%theta_q-parcel1%theta_q_cbase) / &
!                 (parcel1%theta_q_ctop-parcel1%theta_q_cbase)
!             parcel1%x_ent=max(0._wp,parcel1%x_ent)
!             parcel1%x_ent=min(1._wp,parcel1%x_ent)
!             
!             ydot(ite)=ydot(ite)+&
!                 parcel1%x_ent/dt*(parcel1%t_ctop-y(ite) + &
!                 lv/cpm*(parcel1%q_ctop-wv)) +&
!                 (1._wp-parcel1%x_ent)/dt*(parcel1%t_cbase-y(ite) + &
!                 lv/cpm*(parcel1%q_cbase-wv))
!         endif 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in rh of parcel                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ydot(irh)=(p-svp_liq(t))*svp_liq(t)*drv
        ydot(irh)=ydot(irh)+svp_liq(t)*wv*ydot(ipr)
        ydot(irh)=ydot(irh)-wv*p*dfsid1(svp_liq,t,1.e0_wp,1.e-8_wp,err)*ydot(ite)
        ydot(irh)=ydot(irh) / (eps1*svp_liq(t)**2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
    end subroutine fparcelwarm
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! derivatives for a cold parcel model                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates rates of change for a the ice part of the parcel model
	!>@param[in] neq: length of solution vector
	!>@param[in] tt: time
	!>@param[in] y: solution vector
	!>@param[inout] ydot: derivates calculated
	!>@param[in] rpar: real data coming in
	!>@param[in] ipar: integer data coming in
    subroutine fparcelcold(neq, tt, y, ydot, rpar, ipar)
        use numerics_type
        use numerics, only : dfsid1,find_pos

        implicit none
        real(wp), intent(inout) :: tt
        real(wp), intent(inout), dimension(neq) :: y, ydot
        integer(i4b), intent(inout) :: neq
        real(wp), intent(inout) :: rpar
        integer(i4b), intent(inout) :: ipar

        ! local variables
        real(wp) :: wv=0._wp, wl=0._wp, wi=0._wp, rm, cpm, &
                  drv=0._wp, dri=0._wp,dri2=0._wp, &
                  rh,t,p,err,sl, w, &
                  te, qve, pe, var, dummy, rhoe, rhop, b, rh_ice

        integer(i4b) :: i, j,iloc, ipartice, ipr, ite, irh, iz,iw

        ipartice=parcel1%n_bin_modew
        ipr=parcel1%ipri
        ite=parcel1%itei
        irh=parcel1%irhi
        iw =parcel1%iwi
        
        ydot(iw)=0._wp

        rh=y(irh)
        t=y(ite)
        p=y(ipr)
        w=y(iw)
    

        ! check there are no negative values
        where(y(1:ipartice).le.0.e1_wp)
            y(1:ipartice)=1.e-22_wp
        end where


        ! calculate mixing ratios from rh, etc
        sl=svp_liq(t)*rh/(p-svp_liq(t)) ! saturation ratio
        sl=(sl*p/(1._wp+sl))/svp_liq(t)
        wv=eps1*rh*svp_liq(t) / (p-svp_liq(t)) ! vapour mixing ratio
        wl=sum(parcel1%npart*parcel1%y(1:ipartice))          ! liquid mixing ratio
        wi=sum(parcel1%npartice*y(1:ipartice))             ! liquid mixing ratio
        rh_ice = wv / ( eps1*svp_ice(t) / (p-svp_ice(t) ) ) ! rh over ice

        ! calculate the moist gas constants and specific heats
        rm=ra+wv*rv
        cpm=cp+wv*cpv+wl*cpw+wi*cpi

        ! now calculate derivatives
        ! adiabatic parcel model
        ydot(ipr)=0._wp      ! hydrostatic equation

        
        ! particle growth rate - mass growth rate
        parcel1%rh_eq=1._wp
        ydot(1:ipartice)=icegrowthrate01(t,p,rh_ice,parcel1%rh_eq,y(1:ipartice), &
            parcel1%mbinice(:,1:n_comps),parcel1%rhobinice,&
            parcel1%phi,parcel1%rhoi,parcel1%nump,parcel1%rime,ipartice) 
        
        ! do not bother if number concentration too small
        do i=1,ipartice
            if(isnan(ydot(i)).or.parcel1%npartice(i).le. 1.e-9_wp) then
              ydot(i)=0._wp
            endif
        enddo

        ! change in vapour content
        drv = -sum(ydot(1:ipartice)*parcel1%npartice)
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in temperature of parcel                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ydot(ite)=rm/p*ydot(ipr)*t/cpm  ! temperature change: expansion
        ydot(ite)=ydot(ite)-ls/cpm*drv ! temp change: sublimation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in rh of parcel                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ydot(irh)=(p-svp_liq(t))*svp_liq(t)*drv
        ydot(irh)=ydot(irh)+svp_liq(t)*wv*ydot(ipr)
        ydot(irh)=ydot(irh)-wv*p*dfsid1(svp_liq,t,1.e0_wp,1.e-8_wp,err)*ydot(ite)
        ydot(irh)=ydot(irh) / (eps1*svp_liq(t)**2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
    end subroutine fparcelcold
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! jacobian for a warm parcel model : dummy subroutine                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine jparcelwarm(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
          use numerics_type
          implicit none
          real(wp) :: t
          real(wp), dimension(neq) :: y
          real(wp), dimension(nrpd, neq) :: pd
          integer(i4b) :: neq, ml, mu, nrpd
          real(wp) :: rpar
          integer(i4b) :: ipar
      
    end subroutine jparcelwarm
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! INP source function                                                          !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the number concentration of INPs using the DeMott et al. (2010)
	!> parameterisation (Predicting global atmospheric ice...
	!>                    https://doi.org/10.1073/pnas.0910818107)
	!>@param[in] t, naer05: temperature, number concentration of aerosols > 0.5 um
	!>@return demott_2010: number concentration of INPs
	function demott_2010(t,naer05)
		implicit none
		real(wp), intent(in) :: t,naer05
		real(wp) :: demott_2010
		real(wp) :: tc
		tc=ttr-t
		! equation 1 from
		! https://www.pnas.org/content/107/25/11217
		! number per std m^3
		demott_2010=min(0.0594_wp*(tc)**3.33_wp * &
		    (naer05/1.e6_wp)**(0.0264_wp*tc+0.0033_wp),naer05)
		
	end function demott_2010
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ice nucleation                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine noncollisional_iceformation(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin, &
                         moments,medges, &
                         t,p,nbins1,ncomps,nbinw,nmoms,nmodes,yice,rh,dt,&
                         sce_flag,mode1_flag) 
    use numerics_type
    use sce, only : calculate_mode1
    implicit none
    real(wp), intent(inout) :: t
    real(wp), intent(in) :: p,rh,dt
    real(wp), dimension(nbinw), intent(inout) :: npart,npartice
    real(wp), dimension(nbinw), intent(in) :: mwat
    real(wp), dimension(nbinw,ncomps), intent(in) :: &
                                          rhobin,nubin,kappabin,molwbin
    real(wp), dimension(2*nbinw,nmoms), intent(inout) :: moments
    integer(i4b), intent(in) :: ncomps,nbinw, nmoms, nmodes, nbins1
    real(wp), dimension(nbinw,ncomps), intent(in) :: mbin2
    real(wp), dimension(nbinw,ncomps+1), intent(inout) :: mbin2_ice
    real(wp), dimension(nbins1+1,nmodes), intent(in) :: medges
    integer(i4b), intent(in) :: sce_flag
    logical, intent(in) :: mode1_flag
    

    real(wp), dimension(nbinw) :: nw,aw,jw,dn01,m01,ns,dw,dd,kappa,rhoat
    real(wp), dimension(nbinw,ncomps) :: dmaer01
    real(wp), dimension(ncomps) :: dmaer01a
    real(wp), dimension(nmoms) :: momtemp

    real(wp), intent(inout), dimension(nbinw) :: yice

    integer(i4b) :: i,j,k, inew, it, ib
    real(wp) :: fracinliq, fracinice, naer05, nprimary
    real(wp) :: n, nt, nb, mt, mb, mnew, nleft, mttot, mbtot, mleft, mall,  &
            mtot_orig, mtot_mt, mtot_mb

    ! (1) calculate the homogeneous nucleation of ice in SC water
    ! (2) calculate the primary nucleation
    ! (3) calculate the mode-1 freezing fragmentation
    ! (4) transfer the moments to the ice phase

    if(t.gt.ttr) return

    m01 = yice*npartice
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! (1) first calculate the ice formation over dt using koop et al. 2000   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! number of moles of water
    nw = mwat / molw_water
    ! activity of water
    select case(kappa_flag)
        case(0)
            aw(:)=(nw(:))/(nw(:)+sum(mbin2(:,:)/molwbin(:,:)*nubin(:,:),2) )
        case(1)
            rhoat(:)=mwat(:)/rhow+sum(mbin2(:,:)/rhobin(:,:),2)
            rhoat(:)=(mwat(:)+sum(mbin2(:,:),2))/rhoat(:);

            dw(:)=((mwat(:)+sum(mbin2(:,:),2))*6._wp/(pi*rhoat(:)))**(1._wp/3._wp)

            dd(:)=((sum(mbin2(:,:)/rhobin(:,:),2))* &
                    6._wp/(pi))**(1._wp/3._wp) ! dry diameter
                          ! needed for eqn 6, petters and kreidenweis (2007)
            kappa(:)=sum((mbin2(:,:)+1.e-60_wp)/rhobin(:,:)*kappabin(:,:),2) &
                    / sum((mbin2(:,:)+1.e-60_wp)/rhobin(:,:),2)
            ! equation 7, petters and kreidenweis (2007)
            aw=(dw**3-dd**3)/(dw**3-dd**3*(1._wp-kappa)) ! from eq 6,p+k(acp,2007)
        case default
            print *,'error kappa_flag'
            stop
    end select
    ! koop et al. (2000) nucleation rate - due to homogeneous nucleation.
    jw(:)=koopnucrate(aw,t,p,nbinw)
    dn01=0._wp
    ! the number of ice crystals nucleated by homogeneous nucleation:
    dn01(:)=dn01(:)+abs( npart(:)*(1._wp-exp(-jw(:)*mwat(:)/rhow*dt)) )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! (2) second calculate the ice formation over dt using DeMott et al. 2010!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the dry size of the aerosol particle                         
    dd(:)=((sum(mbin2(:,:)/rhobin(:,:),2))*6._wp/(pi))**(1._wp/3._wp) ! dry diameter
    naer05=0._wp
    ! add up the total number > 500 nm
    do i=1,nbinw
        if ((dd(i).gt.0.5e-6_wp).and.(rh.ge.1._wp)) naer05=naer05+npart(i)-dn01(i)
    enddo  
    ! calculate the number of ice crystals according to DeMott et al. 2010
    nprimary=demott_2010(t,naer05)
    ! ensure the number nucleated is less than the number that can nucleate
    nprimary=min(naer05,nprimary) 
    ! this is ice moments for total number of monomers (ncomps+2)
    nprimary=max(nprimary-sum(dn01+moments(nbinw+1:2*nbinw,ncomps+2)),0._wp)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! deplete the aerosols larger than 0.5 microns
    ! (maybe in future will want to deplete the largest first, but at the moment
    !       it is doing it for all modes - in reverse order)
    do i=nbinw,1,-1
        if (nprimary.le.0._wp) exit
        if ((dd(i).gt.0.5e-6_wp).and.(rh.ge.1._wp)) then
            ! it can nucleate: reduce number of primary ice
            dn01(i)=dn01(i)+min(nprimary,npart(i)-dn01(i))
            nprimary=nprimary-min(npart(i)-dn01(i),nprimary)
        endif
    enddo    
    ! limit to number of drops / unfrozen aerosol
    dn01(:)=min(dn01(:),npart(:))
    if(t.gt.ttr) dn01=0._wp
    ! we now have the total number of ice crystals to 'nucleate'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! (3) third calculate fragmentation of these freezing drops              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,nmodes
        do j=1,nbins1
            ! reference to nbinw
            k=j+(i-1)*(nbins1)
            
            ! if these are freezing then calculate mode-1
            ! total number, number of tiny, number of big
            !   and corresponding masses for this freezing event
            mall = mwat(k)*dn01(k)
            mnew = mwat(k)
            nleft = dn01(k)
            n=0._wp
            nt=0._wp
            nb=0._wp
            mt=0._wp
            mb=0._wp
            if (dn01(k) > 0._wp) then
                if(mode1_flag) &
                    call calculate_mode1(mwat(k),0._wp,t,n,nt,nb,mb,mt)
                ! adjust the mass and numbers of the freezing drops in 
                ! the three modes so mass is conserved
                mnew = mnew-nt*mt-nb*mb
                n=n*dn01(k) ! the number concentration of these fragments
                nt=nt*dn01(k)
                nb=nb*dn01(k)
                !nleft=nleft-n ! ice mass left over to go in new bin
            endif     
            ! total mass going into ice bins
            mttot=mt*nt
            mbtot=mb*nb 
            mleft=nleft*mnew
                        
            ! number conc. of liquid bins:
            npart(k) = npart(k) - dn01(k)
            ! fraction of liquid left
            fracinliq = npart(k) / max(npart(k)+dn01(k),1.e-30_wp)
            ! aerosol moments to go into ice
            momtemp(1:ncomps) = (1._wp-fracinliq)*moments(k,1:ncomps)
            ! scale aerosol moments according to this fraction
            moments(k,1:ncomps)=moments(k,1:ncomps)*fracinliq
            
            
            
            
            ! find the bins where mnew, mt, and mb need to go.
            inew    = find_medge(medges,mnew,nbins1,nmodes,i)
            it      = find_medge(medges,mt,nbins1,nmodes,i)
            ib      = find_medge(medges,mb,nbins1,nmodes,i)
            
            ! for this ice bin we need to create three new ice bins in
            ! inew, it, and ib
            if(mall>0._wp) then
                ! aerosol moments going into these bins (by mass)
                moments(inew+nbinw,1:ncomps)=moments(inew+nbinw,1:ncomps)+ &
                                            momtemp(:)*mleft/mall
                moments(it+nbinw,1:ncomps)=moments(it+nbinw,1:ncomps)+ &
                                            momtemp(:)*mttot/mall
                moments(ib+nbinw,1:ncomps)=moments(ib+nbinw,1:ncomps)+ &
                                            momtemp(:)*mbtot/mall
            endif            
            ! the ice mass in these bins
            m01(inew)=m01(inew) + mleft
            m01(it)=m01(it) + mttot
            m01(ib)=m01(ib) + mbtot

            ! number in these bins
            npartice(inew) = npartice(inew) + nleft
            npartice(it) = npartice(it) + nt
            npartice(ib) = npartice(ib) + nb
            
            
            ! now the ice moments: phi, nmon, vol, rim, unfro
            ! phi
            moments(inew+nbinw,ncomps+1)   =moments(inew+nbinw,ncomps+1)+nleft
            moments(it+nbinw,ncomps+1)     =moments(it+nbinw,ncomps+1)+nt
            moments(ib+nbinw,ncomps+1)     =moments(ib+nbinw,ncomps+1)+nb
            ! nmon
            moments(inew+nbinw,ncomps+2)   =moments(inew+nbinw,ncomps+2)+nleft
            moments(it+nbinw,ncomps+2)     =moments(it+nbinw,ncomps+2)+nt
            moments(ib+nbinw,ncomps+2)     =moments(ib+nbinw,ncomps+2)+nb
            ! vol
            moments(inew+nbinw,ncomps+3) =moments(inew+nbinw,ncomps+3)+mleft/rhoice
            moments(it+nbinw,ncomps+3)   =moments(it+nbinw,ncomps+3)+mttot/rhoice
            moments(ib+nbinw,ncomps+3)   =moments(ib+nbinw,ncomps+3)+mbtot/rhoice
            
            
            ! mass of aerosol going to ice - needs to be shared over the bins 
            ! where the ice goes
            !dmaer01a(:)=mbin2(k,:)*dn01(k)
            
            ! ice moments
            ! mass of ice 
            
        enddo
    enddo
    
    where((m01.gt.0._wp).and.(npartice.gt.0._wp))
        yice=m01/npartice
    end where
    ! aerosol mass in ice bins
!     mbin2_ice(:,1:ncomps)=dmaer01(:,:)/(1.e-50_wp+spread(npartice,2,ncomps))
    mbin2_ice(:,1:ncomps)=moments(nbinw+1:2*nbinw,1:ncomps)/(1.e-50_wp+spread(npartice,2,ncomps))
    !moments(1+nbinw:2*nbinw,1:ncomps)=dmaer01(:,:)
    mbin2_ice(:,ncomps+1)=yice
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! latent heat of fusion (these freeze so this is how much heat is released):
    t=t+lf/cp*sum(mwat(:)*dn01(:))
    end subroutine noncollisional_iceformation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! helper function                                                              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function find_medge(medges,m,nbins1,nmodes,imode)
        use numerics_type
        implicit none
        real(wp), dimension(nbins1+1,nmodes), intent(in) :: medges
        real(wp), intent(in) :: m
        integer(i4b), intent(in) :: nbins1,nmodes,imode
        
        integer(i4b) :: find_medge
        integer(i4b) :: i,j
        do i = 1,nbins1
            if(medges(i,imode)>=m) exit
        enddo
        i=i-1
        find_medge=i+(imode-1)*(nbins1)
        
    
    end function find_medge    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ice nucleation                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine icenucleation(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin, &
                         moments,t,p,sz,sz2,sz3,yice,rh,dt) 
      use numerics_type
      use sce, only : calculate_mode1
      implicit none
      real(wp), intent(inout) :: t
      real(wp), intent(in) :: p,rh,dt
      real(wp), dimension(sz2), intent(inout) :: npart,npartice
      real(wp), dimension(sz2), intent(in) :: mwat
      real(wp), dimension(sz2,sz), intent(in) :: &
                                              rhobin,nubin,kappabin,molwbin
      real(wp), dimension(2*sz2,sz3), intent(inout) :: moments
      integer(i4b), intent(in) :: sz,sz2, sz3
      real(wp), dimension(sz2,sz), intent(in) :: mbin2
      real(wp), dimension(sz2,sz+1), intent(inout) :: mbin2_ice

      real(wp), dimension(sz2) :: nw,aw,jw,dn01,m01,ns,dw,dd,kappa,rhoat
      real(wp), dimension(sz2,sz) :: dmaer01
      
      real(wp), intent(inout), dimension(sz2) :: yice
      
      integer(i4b) :: i
      real(wp) :: fracinliq, fracinice, naer05, nprimary
      
      
      ! todo: 
      ! (1) break this into a do loop to make it easier
      ! (2) for each drop that freezes calculate the mode 1 break-up 
      ! (3) put the fragments and moments in the new categories
      !     they go into the same mode as the "mode" of the drop that freezes
      !     maybe loop backwards from the current size
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! first calculate the ice formation over dt using koop et al. 2000       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! number of moles of water
      nw(:)=mwat(:)/molw_water
      ! activity of water
      select case(kappa_flag)
        case(0)
          aw(:)=(nw(:))/(nw(:)+sum(mbin2(:,:)/molwbin(:,:)*nubin(:,:),2) )
        case(1)
          rhoat(:)=mwat(:)/rhow+sum(mbin2(:,:)/rhobin(:,:),2)
          rhoat(:)=(mwat(:)+sum(mbin2(:,:),2))/rhoat(:);
  
          dw(:)=((mwat(:)+sum(mbin2(:,:),2))*6._wp/(pi*rhoat(:)))**(1._wp/3._wp)
  
          dd(:)=((sum(mbin2(:,:)/rhobin(:,:),2))* &
                 6._wp/(pi))**(1._wp/3._wp) ! dry diameter
                              ! needed for eqn 6, petters and kreidenweis (2007)
          kappa(:)=sum((mbin2(:,:)+1.e-60_wp)/rhobin(:,:)*kappabin(:,:),2) &
                 / sum((mbin2(:,:)+1.e-60_wp)/rhobin(:,:),2)
                 ! equation 7, petters and kreidenweis (2007)
          aw=(dw**3-dd**3)/(dw**3-dd**3*(1._wp-kappa)) ! from eq 6,p+k(acp,2007)
        case default
          print *,'error kappa_flag'
          stop
      end select
      ! koop et al. (2000) nucleation rate - due to homogeneous nucleation.
      jw(:)=koopnucrate(aw,t,p,sz2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dn01=0._wp
      
      ! the number of ice crystals nucleated by homogeneous nucleation:
      dn01(:)=dn01(:)+abs( npart(:)*(1._wp-exp(-jw(:)*mwat(:)/rhow*dt)) )


      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculate the dry size of the aerosol particle                         !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dd(:)=((sum(mbin2(:,:)/rhobin(:,:),2))* &
             6._wp/(pi))**(1._wp/3._wp) ! dry diameter
      naer05=0._wp
      do i=1,sz2
        if ((dd(i).gt.0.5e-6_wp).and.(rh.ge.1._wp)) naer05=naer05+npart(i)-dn01(i)
      enddo  
      ! calculate the number of ice crystals
      nprimary=demott_2010(t,naer05)
      ! ensure the number nucleated is less than the number that can nucleate
      nprimary=min(naer05,nprimary) 
      nprimary=max(nprimary-sum(dn01+moments(sz2+1:2*sz2,sz+2)),0._wp)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! deplete the aerosols larger than 0.5 microns
      ! maybe in future will want to deplete the largest first
      do i=sz2,1,-1
        if (nprimary.le.0._wp) exit
        if ((dd(i).gt.0.5e-6_wp).and.(rh.ge.1._wp)) then
            ! it can nucleate
            ! reduce number of primary ice
            dn01(i)=dn01(i)+min(nprimary,npart(i)-dn01(i))
            nprimary=nprimary-min(npart(i)-dn01(i),nprimary)
        endif
      enddo
      
      
      
      dn01(:)=min(dn01(:),npart(:))
      
      if(t.gt.ttr) then
          dn01=0._wp
      endif
      !!!!
      ! total aerosol mass in each bin added together:
      dmaer01(:,:)=(mbin2_ice(:,1:sz)*(spread(npartice(:),2,sz)+1.e-50_wp)+ &
                      mbin2(:,:)*spread(dn01(:),2,sz) ) 
      ! total water mass that will be in the ice bins:
      m01=(yice*npartice+mwat(:)*dn01(:)) 

      ! number conc. of liquid bins:
      npart(:)=npart(:)-dn01(:)
      ! number conc. of ice bins:
      npartice(:)=npartice(:)+dn01(:)
      ! new ice mass in bin:
      m01=m01/(npartice) 
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! transfer moments to ice                                                !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,sz2
          ! aerosol
          fracinliq=npart(i)/max(npart(i)+dn01(i),1.e-30_wp)
          fracinice=npartice(i)/max(npartice(i)-dn01(i),1.e-30_wp)
          moments(i,1:sz) = moments(i,1:sz)*fracinliq
          moments(i+sz2,1:sz) = moments(i+sz2,1:sz)*fracinice
          ! phi, nmon, vol, rim, unf
          ! phi is 1 times number concentration
          moments(i+sz2,sz+1) = moments(i+sz2,sz+1)+dn01(i) 
          ! nmon is 1 times number concentration
          moments(i+sz2,sz+2) = moments(i+sz2,sz+2)+dn01(i) 
          ! volume is number concentration times volume of these crystals
          if (m01(i) .gt. 0._wp) then
              moments(i+sz2,sz+3) = moments(i+sz2,sz+3)+dn01(i)*mwat(i)/rhoice
          else
              moments(i+sz2,sz+3)=0._wp
          endif
      enddo
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      where((m01.gt.0._wp).and.(npartice.gt.0._wp))
        yice=m01
!       elsewhere
!         yice=yice
      end where
      
      ! aerosol mass in ice bins
      mbin2_ice(:,1:sz)=dmaer01(:,:)/(1.e-50_wp+spread(npartice,2,sz))
      moments(1+sz2:2*sz2,1:sz)=dmaer01(:,:)
      mbin2_ice(:,sz+1)=yice
      
      ! latent heat of fusion:
      t=t+lf/cp*sum(mwat(:)*dn01(:))
    end subroutine icenucleation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! calculate new volume and phi                                                       !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief calculates ice growth model of Chen and Lamb (1994) 
	!>@param[in] t,qv, qvsat,rhoa,dm,gamma_t,dep_density
	!>@param[inout] v,phi
    subroutine chen_and_lamb_prop(dm,gamma_t,v,phi, dep_density)
        implicit none
        real(wp), intent(in) :: dm, gamma_t,dep_density
        real(wp), intent(inout) :: v, phi

        real(wp) :: deltaV,v_old,rgamma_tp2,ln_vn_vo
        integer(i4b) :: i
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! increment to volume of crystals - see equation 41
        ! note that this will be per kg of air, rather than crystal but, since we are 
        ! taking the ratio to determine c and a-axes, it should not matter
        deltaV=dm/dep_density
        v_old=v
        v=v+deltaV
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! solving equations 43 and 43 over dV
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rgamma_tp2=1._sp/(gamma_t+2._sp)
        ln_vn_vo=log(v/v_old)
        phi=phi*exp((gamma_t-1._sp)*rgamma_tp2*ln_vn_vo)       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    end subroutine chen_and_lamb_prop
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! moving centre binning                                                              !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief rebins according to the moving centre scheme - Jacobson's book
    subroutine moving_centre(n_bin_mode,n_bin_modew,n_binst,n_mode, &
                    n_comps, n_moments, &
                    npart, masses, moments, mbin,mbinedges) 
        implicit none
        integer(i4b), intent(in) :: n_bin_mode, n_bin_modew, n_binst, n_mode, &
                                n_comps, n_moments
        real(wp), dimension(n_binst+1,n_mode), intent(in) :: mbinedges
        real(wp), dimension(n_bin_modew), intent(inout) :: npart, masses
        real(wp), dimension(n_bin_modew,n_moments), intent(inout) :: moments
        real(wp), dimension(n_bin_modew,n_comps+1), intent(inout) :: mbin
        
        logical, dimension(n_bin_modew) :: moment_exists
        integer(i4b) :: i,j,thismode, thisbin,newplace
        real(wp), dimension(n_bin_modew,n_moments) :: momtemp
        real(wp), dimension(n_bin_modew) :: nparttemp, totmass
        
        momtemp=0._wp ! zero
        nparttemp=0._wp ! zero
        totmass=0._wp ! zero

        ! loop over the "warm" moments
        do i=1,n_bin_modew
            thismode=(i-1)/n_binst+1                ! this is the mode of i
            thisbin=modulo(i-1,n_binst)+1           ! this is the bin
            if(npart(i).gt.0._wp) then
                moment_exists(i)=.true.
                ! if the mass is in the right bin, just add it
                if ((masses(i).gt.mbinedges(thisbin,thismode)).and. &
                    (masses(i).le.mbinedges(thisbin+1,thismode))) then
                    momtemp(i,:)=momtemp(i,:)+moments(i,:)
                    nparttemp(i)=nparttemp(i)+npart(i)
                    totmass(i)=totmass(i)+npart(i)*masses(i)
                else
                    ! if the current mass is not in the correct bin
                    ! find the bin it should be in
                    do j=1,n_binst
                        if ((masses(i).gt.mbinedges(j,thismode)).and. &
                            (masses(i).le.mbinedges(j+1,thismode))) then
                        
                            newplace=(thismode-1)*n_binst+j
!                             print *,newplace, j, i,thisbin,thismode,masses(i), mbinedges(j,thismode), &
!                                 mbinedges(j+1,thismode)
                            ! add the moment from i to newplace
                            momtemp(newplace,:)=momtemp(newplace,:)+moments(i,:)
                            nparttemp(newplace)=nparttemp(newplace)+npart(i)
                            totmass(newplace)=totmass(newplace)+masses(i)*npart(i)
                            
                        endif
                    enddo
                    
                endif                
                
            endif
        
        enddo
        
        
        ! now, set the masses, moments in second loop
        do i=1,n_bin_modew
            thismode=(i-1)/n_binst+1                ! this is the mode of i
            thisbin=modulo(i-1,n_binst)+1           ! this is the bin

            if(nparttemp(i).gt.0._wp) then
                moments(i,:)=momtemp(i,:)
                npart(i)=nparttemp(i)
                masses(i)=totmass(i) / npart(i)
                mbin(i,1:n_comps)=moments(i,1:n_comps)/npart(i)
                if (.not.((masses(i).gt.mbinedges(thisbin,thismode)).and. &
                    (masses(i).le.mbinedges(thisbin+1,thismode)))) then
                    
                    masses(i)=0.5_wp*(mbinedges(thisbin,thismode)+ &
                                    mbinedges(thisbin+1,thismode))
                    npart(i)=totmass(i) / masses(i)
                endif                
                mbin(i,n_comps+1)=masses(i)
            else
                moments(i,:)=0._wp
                npart(i)=0._wp
                masses(i)=0.5_wp*(mbinedges(thisbin,thismode)+ &
                                    mbinedges(thisbin+1,thismode))
                mbin(i,1:n_comps)=0._wp
                mbin(i,n_comps+1)=masses(i)
            endif
        enddo
        
        
        
        
    end subroutine moving_centre
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Chen and Lamb (1994) ancillary variables                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief calculates ancillary variables for ice growth model of Chen and Lamb (1994) 
	!>@param[in] t,qv, qvsat,rhoa
	!>@param[inout] v,phi,gamma_t,dep_density
    subroutine chen_and_lamb_anc(t,qv,qvsat,rhoa,gamma_t, dep_density)
        implicit none
        real(wp), intent(in) :: t,qv,qvsat,rhoa
        real(wp), intent(inout) :: gamma_t,dep_density

        real(wp) :: delta_rho,t1
        integer(i4b) :: i
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! calculate the inherent growth ratio - this is from a 17th order polynomial
        gamma_t=0._sp
        t1=min(max(t,243.15),273.15) ! range of fit
        do i=1,n_cl
            gamma_t=gamma_t+((t1-gam_mu_cl(1))/gam_mu_cl(2))**(n_cl-i)*gam_cl(i)
        enddo
        gamma_t=10._sp**gamma_t
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! equation 42 from Chen and Lamb (1994, JAS: The Theoretical Basis for 
        !   Parameterisation of Ice Crystal Habits)
        delta_rho=(qv-qvsat)*rhoa*1000._sp ! g/m^3
        dep_density=rhoice*exp(-3._sp*max(delta_rho-0.05_sp,0._sp)/gamma_t)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine chen_and_lamb_anc
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Chen and Lamb (1994) capacitance factors                                           !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief calculates the ratio of capacitance to that of an equivalent spehre
	!  Chen and Lamb (1994) 
	!>@param[in] phi
	!>@return cap_fac
    function chen_and_lamb_cap_fac(phi)
        implicit none
        real(wp), intent(in) :: phi
        real(wp) :: chen_and_lamb_cap_fac
        real(wp) :: fac1,fac2,ecc

        
        ! factor to convert between R and a - derived from equating volume of sphere to 
        ! volume of spheroid and taking the ratio of a / r
        fac1=(1._sp/(phi))**(1/3)
        
        ! factor to convert between a and capacitance
        if(phi<0.99_sp) then
            ! see equation 39 of Chen and Lamb (1994)
            ecc=sqrt(1._sp-phi**2)
            fac2=ecc/asin(ecc)
        elseif(phi>1.01_sp) then
            ! see equation 40 of Chen and Lamb (1994)
            ecc=sqrt(1._sp-(1._sp/phi)**2)
            fac2=1._sp/phi/log((1._sp+ecc)*phi)
        else
            fac2=1._sp
        endif
        
        ! total factor
        chen_and_lamb_cap_fac=fac1*fac2
        
    end function chen_and_lamb_cap_fac
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! perform mass balance                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculate current mass
    subroutine mass_balance(neq,neqice,y,yice,npart,npartice, &
                        mass1,n_bin_modew,irh,ite,ipr,ice_flag)
    implicit none
    integer(i4b), intent(in) :: neq,neqice,n_bin_modew,irh,ite,ipr,ice_flag
    real(wp), dimension(neq), intent(in) :: y
    real(wp), dimension(n_bin_modew), intent(in) :: npart, npartice
    real(wp), dimension(neqice), intent(in) :: yice
    real(wp), intent(inout) :: mass1
    
        ! total water before:
        mass1=sum(npart*y(1:n_bin_modew))+ &
            y(irh)*eps1* &
            svp_liq(y(ite)) / &
            (y(ipr)-svp_liq(y(ite)))
        if(ice_flag .eq. 1) then
            mass1=mass1+sum(npartice*yice(1:n_bin_modew))
        endif
        
    end subroutine mass_balance 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       
       
       
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! adjust the relative humidity for mass balance                                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>subtract the difference from the water vapour and compute new RH
    subroutine adjust_relative_humidity(mass1,mass2,vapour_mass,t,p, rh)
        implicit none
        real(wp), intent(in) :: mass1, mass2, t, p
        real(wp), intent(inout) :: vapour_mass, rh
        
        real(wp) :: deltam

        deltam=mass2-mass1
        vapour_mass=rh*eps1*svp_liq(t) / (p-svp_liq(t))
        ! adjust to conserve:
        vapour_mass=vapour_mass-deltam
        rh=vapour_mass / ( eps1*svp_liq(t) / (p-svp_liq(t)) )
            
    end subroutine adjust_relative_humidity
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the ice properties and vapour growth conditions                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>using moments, calculate the ice growth properties
    subroutine ice_vapour_growth_properties(rhoa,qvsat,qv,&
        neqice,yice,n_bin_modew,npartice,phi,nump,rhoi, &
        n_bin_mode,n_moms,n_comps,moments, &
        gamma_t,dep_density,t,p,rhi)
    
        implicit none
        real(wp), intent(in) :: t,p,rhi
        real(wp), intent(inout) :: rhoa, qvsat,qv,gamma_t,dep_density
        integer(i4b), intent(in) :: neqice,n_bin_modew,n_bin_mode,n_moms,n_comps
        real(wp), dimension(neqice), intent(in) :: yice
        real(wp), dimension(n_bin_modew), intent(inout) :: npartice, phi, nump, rhoi
        real(wp), dimension(n_bin_mode,n_moms), intent(in) :: moments
        
        
        rhoa = p / (ra*t)                             ! density of air
        qv = rhi*eps1*svp_liq(t) / (p-svp_liq(t))        ! qv_sat
        qvsat = (eps1*svp_ice(t) / (p-svp_ice(t)))

        call chen_and_lamb_anc(t, &
                qv,qvsat,rhoa,gamma_t, dep_density) ! set the deposition density
        
        phi=moments(n_bin_modew+1:n_bin_mode,n_comps+1) &
            / moments(n_bin_modew+1:n_bin_mode,n_comps+2)
        
        nump=moments(n_bin_modew+1:n_bin_mode,n_comps+2) &
            / npartice

        rhoi=(npartice*yice(1:n_bin_modew)-&
            moments(n_bin_modew+1:n_bin_mode,n_comps+4)) &
            / moments(n_bin_modew+1:n_bin_mode,n_comps+3)

    end subroutine ice_vapour_growth_properties         
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Cloud-top entrainment according to sanchez et al 2017                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>using moments, calculate the ice growth properties
    subroutine cloud_top_entrainment(n_levels_s, n_sound, neq, n_bin_modew, &
            rh, t, p, z, z_cbase, theta_q_ctop, q_ctop, theta_q_cbase, q_cbase, theta_q, &
            x_ent, y, npart, z_sound, theta_q_sound, set_theta_q_cb_flag)   
        use numerics, only : zeroin, dvode
        implicit none
        integer(i4b), intent(in) :: n_levels_s, n_sound, neq, n_bin_modew
        real(wp), intent(in) :: p,z, z_cbase, theta_q_ctop, q_ctop
        real(wp), intent(inout) :: theta_q_cbase, q_cbase, theta_q, x_ent, rh, t
        real(wp), dimension(neq), intent(in) :: y
        real(wp), dimension(n_bin_modew), intent(in) :: npart
        real(wp), dimension(n_sound), intent(in) :: z_sound, theta_q_sound
        logical, intent(inout) :: set_theta_q_cb_flag
        
        real(wp) :: vapour_mass, liquid_mass, mass2, var, dummy, cpm, x1, x2, x2old
        integer(i4b) :: iloc
        
        ! vapour mass in the parcel:
        vapour_mass=rh*eps1* svp_liq(t) / (p-svp_liq(t))

        ! liquid mass in the parcel:
        liquid_mass=sum(npart*y(1:n_bin_modew))
        ! total water in the parcel:
        mass2=liquid_mass+ vapour_mass

        ! this is the theta q at cloud base, given the total water content
        if (set_theta_q_cb_flag) then
            theta_q_cbase=calc_theta_q3(t, p, mass2)
                                    
            ! locate position
            iloc=find_pos(z_sound(1:n_levels_s),z)
            iloc=min(n_levels_s-1,iloc)
            iloc=max(1,iloc)
            ! linear interp theta_q
            call poly_int(z_sound(iloc:iloc+1), theta_q_sound(iloc:iloc+1), &
                    min(z,z_sound(n_levels_s)), var,dummy)        
            theta_q_cbase=var

                                    
            q_cbase=mass2 ! total water is all vapour at cloud base
            if(z .gt. z_cbase) then
                set_theta_q_cb_flag=.false.
            endif
        endif
        
        ! do the entrainment only above cloud-base
        if(.not. set_theta_q_cb_flag) then
            cpm=cp+vapour_mass*cpv+liquid_mass*cpw
        
            ! the new value of theta_q
            theta_q=calc_theta_q3(t,  p, mass2)
            ! locate position
            iloc=find_pos(z_sound(1:n_levels_s),z)
            iloc=min(n_levels_s-1,iloc)
            iloc=max(1,iloc)
            ! linear interp theta_q
            call poly_int(z_sound(iloc:iloc+1), theta_q_sound(iloc:iloc+1), &
                    min(z,z_sound(n_levels_s)), var,dummy)        
            theta_q=var
                                        
            ! equation 5 (Sanchez et al, 2017, acp)
            x_ent=(theta_q-theta_q_cbase) / &
                    (theta_q_ctop-theta_q_cbase)
            x1=min(max(x_ent,0._wp),1._wp)
            x2=1._wp-x1
            
            !print *,parcel1%theta_q,parcel1%theta_q_cbase,parcel1%theta_q_ctop
            ! entrainment of vapour:
            vapour_mass= &!vapour_mass+&
                x1*(q_ctop)+&
                x2*(q_cbase)-liquid_mass
            ! entrainment of liquid:
            !parcel1%npart=parcel1%npart*x2/x2old

            ! entrainment of theta_q (equation 4, Sanchez et al, 2017, ACP):
            theta_q=x1*theta_q_ctop + &
                x2*theta_q_cbase
            
            p111=p
            theta_q_sat=theta_q
            t=zeroin(150._wp, theta_q_sat, calc_theta_q,1.e-30_wp)

            ! rh
            rh=vapour_mass / &
                 ( eps1*svp_liq(t) / &
                 (p-svp_liq(t)) )
            x2old=max(x2,1.e-20_wp)
            print *,x1,x2
        endif
        
        
    end subroutine cloud_top_entrainment        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   
           
                
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update volume and shape in moments                                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>update volume and shape of each crystal
    subroutine update_volume_and_shape(n_bin_modew,n_bin_mode,n_moments,n_comps, &
        momtemp,moments,neqice,yice,yoldice,gamma_t,dep_density,npartice)
    
        implicit none
        integer(i4b), intent(in) :: n_bin_modew,n_bin_mode,n_moments,n_comps, neqice
        real(wp), intent(in) :: gamma_t, dep_density
        real(wp), dimension(neqice), intent(in) :: yice, yoldice
        real(wp), dimension(n_bin_mode,n_moments), intent(inout) :: moments
        real(wp), dimension(n_bin_mode), intent(inout) :: momtemp
        real(wp), dimension(n_bin_modew), intent(in) :: npartice
        
        integer(i4b) :: i
                
        do i=1,n_bin_modew
            ! new total volume
            momtemp(i)= &
                moments(n_bin_modew+i,n_comps+3)+ &
                    (yice(i)-yoldice(i))/dep_density*npartice(i)
            ! new phi
            if (momtemp(i)> 0._wp) then
                moments(n_bin_modew+i,n_comps+1)= moments(n_bin_modew+i,n_comps+1)*&
                    exp((gamma_t-1._wp)/(gamma_t+2._wp)*&
                        log(momtemp(i) / moments(n_bin_modew+i,n_comps+3)))
            else
                moments(n_bin_modew+i,n_comps+1)= npartice(i)
            endif
        enddo
        ! save new volume
        moments(n_bin_modew+1:n_bin_mode,n_comps+3)= momtemp(1:n_bin_modew)                
    end subroutine update_volume_and_shape         
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! reduce the rime mass in proportion during evaporation                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>rime mass reduces in proportion to total mass
    subroutine reduce_rime(ipart, massnew, massold, rimemass)
        implicit none
        integer(i4b), intent(in) :: ipart
        real(wp), intent(in), dimension(ipart) :: massold, massnew
        real(wp), intent(inout), dimension(ipart) :: rimemass
        
        rimemass = rimemass * min(massnew, massold)/massold            
            
    end subroutine reduce_rime
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
                
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! one time-step of the bin-microphysics                                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates one time-step of bin-microphysics
    subroutine bin_microphysics(func1,func2,func3,func4)
    use numerics_type
    use numerics, only : zeroin, dvode
    implicit none
    real(wp) :: mass1, mass2, deltam, vapour_mass, liquid_mass, x1,x2 , cpm, &
        var, dummy, gamma_t, dep_density, rhoa, qv, qvsat
    integer(i4b) :: iloc, i
    
    interface
        subroutine func1(neq, tt, y, ydot, rpar, ipar)
            use numerics_type
            use numerics, only : dfsid1,find_pos

            implicit none
            real(wp), intent(inout) :: tt
            real(wp), intent(inout), dimension(neq) :: y, ydot
            integer(i4b), intent(inout) :: neq
            real(wp), intent(inout) :: rpar
            integer(i4b), intent(inout) :: ipar
        end subroutine func1
    end interface
    interface
        subroutine func2(neq, tt, y, ydot, rpar, ipar)
            use numerics_type
            use numerics, only : dfsid1,find_pos

            implicit none
            real(wp), intent(inout) :: tt
            real(wp), intent(inout), dimension(neq) :: y, ydot
            integer(i4b), intent(inout) :: neq
            real(wp), intent(inout) :: rpar
            integer(i4b), intent(inout) :: ipar
        end subroutine func2
    end interface
    interface
        subroutine func3(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin,moments, &
                         t,p,sz,sz2,sz3,yice,rh,dt) 
            use numerics_type
            implicit none
            real(wp), intent(inout) :: t
            real(wp), intent(in) :: p,rh,dt
            real(wp), dimension(sz2), intent(inout) :: npart,npartice
            real(wp), dimension(sz2), intent(in) :: mwat
            real(wp), dimension(sz2,sz), intent(in) :: mbin2, &
                                                  rhobin,nubin,kappabin,molwbin
            real(wp), dimension(2*sz2,sz3), intent(inout) :: moments
            integer(i4b), intent(in) :: sz,sz2, sz3
            real(wp), dimension(sz2,sz+1), intent(inout) :: mbin2_ice
            real(wp), intent(inout), dimension(sz2) :: yice
        end subroutine func3
    end interface
    interface
        subroutine func4(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin,moments,medges, &
                         t,p,nbins1,ncomps,nbinw,nmoms,nmodes,yice,rh,dt,sce_flag, &
                         mode1_flag) 
            use numerics_type
            implicit none
            real(wp), intent(inout) :: t
            real(wp), intent(in) :: p,rh,dt
            real(wp), dimension(nbinw), intent(inout) :: npart,npartice
            real(wp), dimension(nbinw), intent(in) :: mwat
            real(wp), dimension(nbinw,ncomps), intent(in) :: mbin2, &
                                                  rhobin,nubin,kappabin,molwbin
            real(wp), dimension(2*nbinw,nmoms), intent(inout) :: moments
            integer(i4b), intent(in) :: ncomps,nbinw, nmoms, nmodes, nbins1
            real(wp), dimension(nbinw,ncomps+1), intent(inout) :: mbin2_ice
            real(wp), intent(inout), dimension(nbinw) :: yice
            real(wp), intent(in), dimension(nbins1+1,nmodes) :: medges
            integer(i4b), intent(in) :: sce_flag
            logical, intent(in) :: mode1_flag
        end subroutine func4
    end interface
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mass balance                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(adiabatic_prof) then
        parcel1%yold=parcel1%y ! store old
        call mass_balance(parcel1%neq,parcel1%neqice,parcel1%y,&
                    parcel1%yice,parcel1%npart,parcel1%npartice, &
                        mass1,parcel1%n_bin_modew,&
                        parcel1%irh,parcel1%ite,parcel1%ipr,ice_flag)
    endif    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ODE solver                                                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parcel1%tout=parcel1%tt+parcel1%dt
    do while (parcel1%tt .lt. parcel1%tout)
        parcel1%istate=1
        
        call dvode(func1,parcel1%neq,parcel1%y,parcel1%tt,parcel1%tout,&
                   parcel1%itol,parcel1%rtol,parcel1%atol,&
                   parcel1%itask,parcel1%istate,parcel1%iopt,&
                   parcel1%rwork,parcel1%lrw,&
                   parcel1%iwork,parcel1%liw,jparcelwarm, &
                   parcel1%mf,parcel1%rpar,parcel1%ipar)
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (sce_flag.gt.0) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Moving Centre binning                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call moving_centre(parcel1%n_bin_mode,parcel1%n_bin_modew,&
                parcel1%n_bins1,parcel1%n_modes, parcel1%n_comps, &
                parcel1%imoms+parcel1%n_comps, parcel1%npart, &
                parcel1%y(1:parcel1%n_bin_modew), &
                parcel1%moments(1:parcel1%n_bin_modew,:), &
                parcel1%mbin,parcel1%mbinedges)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif

    if(ice_flag .eq. 1) then
        ! ice part of the parcel model
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Ice nucleation                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call func4(parcel1%npart(1:parcel1%n_bin_modew), &
                parcel1%npartice(1:parcel1%n_bin_modew), &
                parcel1%y(1:parcel1%n_bin_modew), &
                parcel1%mbin(:,1:n_comps), &
                parcel1%mbinice(:,1:n_comps+1), &
                parcel1%rhobin(:,1:n_comps), &
                parcel1%nubin(:,1:n_comps), &
                parcel1%kappabin(:,1:n_comps), &
                parcel1%molwbin(:,1:n_comps), &
                parcel1%moments(1:2*parcel1%n_bin_modew,1:parcel1%imoms), &
                parcel1%mbinedges(:,:), &
                parcel1%y(parcel1%ite), &
                parcel1%y(parcel1%ipr),&
                parcel1%n_bins1, &
                n_comps,parcel1%n_bin_modew,parcel1%imoms+n_comps,parcel1%n_modes, &
                parcel1%yice(1:parcel1%n_bin_modew), &
                parcel1%y(parcel1%irh), parcel1%dt,sce_flag, &
                mode1_flag) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ODE solver for ice                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        parcel1%toutice=parcel1%ttice+parcel1%dt
        parcel1%yice(parcel1%ipri)=parcel1%y(parcel1%ipr)
        parcel1%yice(parcel1%itei)=parcel1%y(parcel1%ite)
        parcel1%yice(parcel1%irhi)=parcel1%y(parcel1%irh)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! needed for ice crystal shape                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        parcel1%yoldice=parcel1%yice                            ! set the old value of ice
        
        call ice_vapour_growth_properties(rhoa,qvsat,qv,&
            parcel1%neqice,parcel1%yice,parcel1%n_bin_modew,parcel1%npartice,&
            parcel1%phi,parcel1%nump,parcel1%rhoi, &
            parcel1%n_bin_mode,parcel1%n_comps+parcel1%imoms,&
            parcel1%n_comps,parcel1%moments, &
            gamma_t,dep_density,parcel1%y(parcel1%ite),parcel1%yice(parcel1%ipri),&
            parcel1%yice(parcel1%irhi))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        do while (parcel1%ttice .lt. parcel1%toutice)
            parcel1%istateice=1
            call dvode(func2,parcel1%neqice,parcel1%yice,parcel1%ttice,parcel1%toutice,&
                           parcel1%itolice,parcel1%rtolice,parcel1%atolice,&
                           parcel1%itaskice,parcel1%istateice,parcel1%ioptice,&
                           parcel1%rworkice,parcel1%lrwice,&
                           parcel1%iworkice,parcel1%liwice,jparcelwarm, &
                           parcel1%mfice,parcel1%rparice,parcel1%iparice)
        enddo
        ! rime also evaporates off
        call reduce_rime(parcel1%n_bin_modew,parcel1%yice(1:parcel1%n_bin_modew), &
            parcel1%yoldice(1:parcel1%n_bin_modew), &
            parcel1%moments(parcel1%n_bin_modew+1:parcel1%n_bin_mode,parcel1%n_comps+4))

        parcel1%y(parcel1%ipr)=parcel1%yice(parcel1%ipri)
        parcel1%y(parcel1%ite)=parcel1%yice(parcel1%itei)
        parcel1%y(parcel1%irh)=parcel1%yice(parcel1%irhi)
        
        
        ! update volume, and shape
        call update_volume_and_shape(parcel1%n_bin_modew,parcel1%n_bin_mode,&
            parcel1%n_comps+parcel1%imoms,parcel1%n_comps, &
            parcel1%momtemp,parcel1%moments,parcel1%neqice, &
            parcel1%yice,parcel1%yoldice,gamma_t,dep_density,parcel1%npartice)
        ! set number of monomers to come from liquid too?            
        !parcel1%moments(1:parcel1%n_bin_modew,n_comps+2)=parcel1%npart
        ! set mass to come from liquid
        parcel1%moments(1:parcel1%n_bin_modew,n_comps+4)=parcel1%npart* &
            parcel1%y(1:parcel1%n_bin_modew)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Moving Centre binning                                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (sce_flag.gt.0) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Moving Centre binning                                            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call moving_centre(parcel1%n_bin_mode,parcel1%n_bin_modew,&
                    parcel1%n_bins1,parcel1%n_modes, parcel1%n_comps, &
                    parcel1%imoms+parcel1%n_comps, parcel1%npartice, &
                    parcel1%yice(1:parcel1%n_bin_modew), &
                    parcel1%moments(1+parcel1%n_bin_modew:parcel1%n_bin_mode,:), &
                    parcel1%mbinice,parcel1%mbinedges)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    endif



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! stop the simulation if parcel is above cloud-top                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if((parcel1%y(parcel1%iz) .gt. parcel1%z_ctop)  .and. &
        vert_ent) parcel1%break_flag=.true.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! vertical entrainment outside of solver (see Sanchez et al, 2017, ACP)!
    ! (lateral entrainment is inside solver)                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(vert_ent .and. (.not. adiabatic_prof) .and. &
        (parcel1%y(parcel1%iz) .gt. parcel1%z_cbase) ) then
        call cloud_top_entrainment(n_levels_s, parcel1%n_sound, parcel1%neq, &
            parcel1%n_bin_modew, &
            parcel1%y(parcel1%irh), parcel1%y(parcel1%ite), &
            parcel1%y(parcel1%ipr), parcel1%y(parcel1%iz), parcel1%z_cbase, &
            parcel1%theta_q_ctop, parcel1%q_ctop, parcel1%theta_q_cbase, &
            parcel1%q_cbase, parcel1%theta_q, &
            parcel1%x_ent, parcel1%y, parcel1%npart, &
            parcel1%z_sound, parcel1%theta_q_sound, set_theta_q_cb_flag)   
        

        
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mass balance                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(adiabatic_prof) then
        ! total water after:
        call mass_balance(parcel1%neq,parcel1%neqice,parcel1%y,&
                    parcel1%yice,parcel1%npart,parcel1%npartice, &
                        mass2,parcel1%n_bin_modew,&
                        parcel1%irh,parcel1%ite,parcel1%ipr,ice_flag)
        
        ! adjust the humidity - could also adjust temperature - not currently done
        call adjust_relative_humidity(mass1,mass2,vapour_mass, &
            parcel1%y(parcel1%ite), parcel1%y(parcel1%ipr), &
            parcel1%y(parcel1%irh))
        
    endif    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    end subroutine bin_microphysics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    
    
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! HELPER ROUTINE                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check(status)
    use netcdf
    use numerics_type
    integer(I4B), intent ( in) :: status

    if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped"
    end if
    end subroutine check
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! output to netcdf                                                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>output 1 time-step of model
	!>@param[inout] new_file
    subroutine output(new_file,outputfile)

    use numerics_type
    use netcdf

    implicit none
    logical, intent(inout) :: new_file
    character (len=*),intent(in) :: outputfile
    real(sp) :: phi
    ! output to netcdf file
    if(new_file) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! open / create the netcdf file                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call check( nf90_create(outputfile, NF90_CLOBBER, io1%ncid) )
        ! define dimensions (netcdf hands back a handle)
        call check( nf90_def_dim(io1%ncid, "times", NF90_UNLIMITED, io1%x_dimid) )
        call check( nf90_def_dim(io1%ncid, "nbins", n_bins, io1%bin_dimid) )
        call check( nf90_def_dim(io1%ncid, "nbinst", parcel1%n_bins1, io1%bin2_dimid) )
        call check( nf90_def_dim(io1%ncid, "nbinste", parcel1%n_bins1+1, io1%bin3_dimid) )
        call check( nf90_def_dim(io1%ncid, "nmodes", parcel1%n_modes, io1%mode_dimid) )
        call check( nf90_def_dim(io1%ncid, "ncomps", parcel1%n_comps, io1%comp_dimid) )
        

        ! close the file, freeing up any internal netCDF resources
        ! associated with the file, and flush any buffers
        call check( nf90_close(io1%ncid) )


        ! now define some variables, units, etc
        call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
        ! define mode
        call check( nf90_redef(io1%ncid) )


        ! define variable: time
        call check( nf90_def_var(io1%ncid, "time", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "time", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "seconds") )

        ! define variable: z
        call check( nf90_def_var(io1%ncid, "z", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "z", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "m") )

        ! define variable: p
        call check( nf90_def_var(io1%ncid, "p", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "p", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "Pa") )

        ! define variable: t
        call check( nf90_def_var(io1%ncid, "t", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "t", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "K") )

        ! define variable: rh
        call check( nf90_def_var(io1%ncid, "rh", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "rh", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "n/a") )

        ! define variable: w
        call check( nf90_def_var(io1%ncid, "w", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "w", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "m s-1") )
                   
        ! define variable: ql
        call check( nf90_def_var(io1%ncid, "ql", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "ql", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "kg kg-1") )
                   
        ! define variable: extinction
        call check( nf90_def_var(io1%ncid, "beta_ext", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "beta_ext", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "m-1") )
                   
        ! define variable: number > 2.5 microns (8.1812e-15 kg)
        call check( nf90_def_var(io1%ncid, "ndrop", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "ndrop", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "m-3") )
                   
        ! define variable: deff
        call check( nf90_def_var(io1%ncid, "deff", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "deff", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "m") )
                   
        ! define variable: mwat
        call check( nf90_def_var(io1%ncid, "mwat", NF90_DOUBLE, &
                    (/io1%bin_dimid,io1%mode_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "mwat", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "kg") )
                   
        ! define variable: nwat
        call check( nf90_def_var(io1%ncid, "nwat", NF90_DOUBLE, &
                    (/io1%bin2_dimid,io1%mode_dimid, io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "nwat", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "#/kg") )

        ! define variable: maer
        call check( nf90_def_var(io1%ncid, "maer", NF90_DOUBLE, &
                    (/io1%bin2_dimid,io1%mode_dimid,io1%comp_dimid, io1%x_dimid/), &
                        io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "maer", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "kg") )
                   

        ! define variable: mbinedges
        call check( nf90_def_var(io1%ncid, "mbinedges", NF90_DOUBLE, &
                    (/io1%bin3_dimid,io1%mode_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "mbinedges", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "kg") )

                   
                   
        if(ice_flag .eq. 1) then
            ! define variable: qi
            call check( nf90_def_var(io1%ncid, "qi", NF90_DOUBLE, &
                        (/io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "qi", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "kg kg-1") )
                   
                   
            ! define variable: number of ice crystals
            call check( nf90_def_var(io1%ncid, "nice", NF90_DOUBLE, &
                        (/io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "nice", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "m-3") )  
                        
            ! define variable: mice
            call check( nf90_def_var(io1%ncid, "mice", NF90_DOUBLE, &
                        (/io1%bin_dimid,io1%mode_dimid, io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "mice", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "kg") )
                   
            ! define variable: phi
            call check( nf90_def_var(io1%ncid, "phi", NF90_DOUBLE, &
                        (/io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "phi", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "") )  
                        
            ! define variable: nmon
            call check( nf90_def_var(io1%ncid, "nmon", NF90_DOUBLE, &
                        (/io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "nmon", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "") )  
                        
            ! define variable: rhoi
            call check( nf90_def_var(io1%ncid, "rhoi", NF90_DOUBLE, &
                        (/io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "rhoi", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "kg/m3") )  
                        
            ! define variable: nicem
            call check( nf90_def_var(io1%ncid, "nicem", NF90_DOUBLE, &
                        (/io1%bin2_dimid,io1%mode_dimid, io1%x_dimid/), io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "nicem", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "#/kg") )


            ! define variable: maeri
            call check( nf90_def_var(io1%ncid, "maeri", NF90_DOUBLE, &
                        (/io1%bin2_dimid,io1%mode_dimid,io1%comp_dimid, io1%x_dimid/), &
                            io1%varid) )
            ! get id to a_dimid
            call check( nf90_inq_varid(io1%ncid, "maeri", io1%a_dimid) )
            ! units
            call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                       "units", "kg") ) 
                   
                   
        endif
        
        call check( nf90_enddef(io1%ncid) )

        ! write variable: mbinedges
        call check( nf90_inq_varid(io1%ncid, "mbinedges", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, parcel1%mbinedges, &
                    start = (/io1%icur/)))

        call check( nf90_close(io1%ncid) )

        new_file=.false.
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write data to file                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call check( nf90_open(outputfile, NF90_WRITE, io1%ncid) )
    
    ! write variable: time
    call check( nf90_inq_varid(io1%ncid, "time", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, parcel1%tt, &
                start = (/io1%icur/)))
    ! write variable: z
    call check( nf90_inq_varid(io1%ncid, "z", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%iz), &
                start = (/io1%icur/)))

    ! write variable: p
    call check( nf90_inq_varid(io1%ncid, "p", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%ipr), &
                start = (/io1%icur/)))

    ! write variable: t
    call check( nf90_inq_varid(io1%ncid, "t", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%ite), &
                start = (/io1%icur/)))

    ! write variable: rh
    call check( nf90_inq_varid(io1%ncid, "rh", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%irh), &
                start = (/io1%icur/)))

    ! write variable: w
    call check( nf90_inq_varid(io1%ncid, "w", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%iw), &
                start = (/io1%icur/)))


    ! write variable: ql
    call check( nf90_inq_varid(io1%ncid, "ql", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        sum(parcel1%y(1:parcel1%n_bin_modew)*parcel1%npart(1:parcel1%n_bin_modew)), &
                start = (/io1%icur/)))

    ! write variable: beta_ext
    call check( nf90_inq_varid(io1%ncid, "beta_ext", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        2._wp*sum((parcel1%y(1:parcel1%n_bin_modew)* &
            6._wp/(rhow*pi))**(2._wp/3._wp)* pi/4._wp* &
            parcel1%npart(1:parcel1%n_bin_modew)), &
                start = (/io1%icur/)))

    ! write variable: number > 2.5 microns (8.1812e-15 kg)
    parcel1%ndrop=0._wp
    where (parcel1%y(1:parcel1%n_bin_modew) > 6.5450e-14_wp)
        parcel1%ndrop=parcel1%npart(:)
    end where
    
    call check( nf90_inq_varid(io1%ncid, "ndrop", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        sum(parcel1%ndrop), start = (/io1%icur/)))

    ! write variable: effective radius
    call check( nf90_inq_varid(io1%ncid, "deff", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        sum((parcel1%y(1:parcel1%n_bin_modew)* &
            6._wp/(rhow*pi))**(3._wp/3._wp)*  &
            parcel1%npart(1:parcel1%n_bin_modew)) / &
        sum((parcel1%y(1:parcel1%n_bin_modew)* &
            6._wp/(rhow*pi))**(2._wp/3._wp)*  &
            parcel1%npart(1:parcel1%n_bin_modew)), &
                start = (/io1%icur/)))

    call check( nf90_inq_varid(io1%ncid, "mwat", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        reshape(parcel1%y(1:parcel1%n_bin_modew),(/n_bins,n_mode/)), start = (/1,1,io1%icur/)))

    call check( nf90_inq_varid(io1%ncid, "nwat", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        reshape(parcel1%npart(1:parcel1%n_bin_modew), &
        (/parcel1%n_bins1,parcel1%n_modes/)), start = (/1,1,io1%icur/)))


    call check( nf90_inq_varid(io1%ncid, "maer", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        reshape(parcel1%mbin(1:parcel1%n_bin_modew,&
            1:parcel1%n_comps), &
        (/parcel1%n_bins1,parcel1%n_modes,parcel1%n_comps/)), start = (/1,1,1,io1%icur/)))

    if(ice_flag .eq. 1) then
        ! write variable: qi
        call check( nf90_inq_varid(io1%ncid, "qi", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            sum(parcel1%yice(1:parcel1%n_bin_modew)* &
                parcel1%npartice(1:parcel1%n_bin_modew)), &
                    start = (/io1%icur/)))

        ! write variable: number concentration of ice crystals
        parcel1%nice=parcel1%npartice
        call check( nf90_inq_varid(io1%ncid, "nice", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            sum(parcel1%nice), start = (/io1%icur/)))
    
        call check( nf90_inq_varid(io1%ncid, "mice", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            reshape(parcel1%yice(1:parcel1%n_bin_modew),(/n_bins,n_mode/)), &
             start = (/1,1,io1%icur/)))

        ! write variable: phi
        phi=sum(parcel1%moments(parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
            parcel1%n_comps+1)) / sum(parcel1%moments(parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
            parcel1%n_comps+2))
        call check( nf90_inq_varid(io1%ncid, "phi", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            phi, start = (/io1%icur/)))
    
        ! write variable: nmon
        phi=sum(parcel1%moments(parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
            parcel1%n_comps+2)) / sum(parcel1%npartice)
        call check( nf90_inq_varid(io1%ncid, "nmon", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            phi, start = (/io1%icur/)))
    
        ! write variable: rhoi
        phi=sum(parcel1%yice(1:parcel1%n_bin_modew)* &
                parcel1%npartice(1:parcel1%n_bin_modew) - &
                parcel1%moments(parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
            parcel1%n_comps+4)) / &
                sum(parcel1%moments(parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
                parcel1%n_comps+3))
        call check( nf90_inq_varid(io1%ncid, "rhoi", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            phi, start = (/io1%icur/)))
            
        call check( nf90_inq_varid(io1%ncid, "nicem", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            reshape(parcel1%npartice(1:parcel1%n_bin_modew), &
            (/parcel1%n_bins1,parcel1%n_modes/)), start = (/1,1,io1%icur/)))

        call check( nf90_inq_varid(io1%ncid, "maeri", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            reshape(parcel1%mbinice(:,&
                1:parcel1%n_comps), &
            (/parcel1%n_bins1,parcel1%n_modes,parcel1%n_comps/)), start = (/1,1,1,io1%icur/)))

            
            
    endif
    

    call check( nf90_close(io1%ncid) )


    io1%icur=io1%icur+1
    end subroutine output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! map to sce                                                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>map the BMM to the SCE variables
	!>@param[in] ice_flag - flag to say if we are computing ice
	subroutine map_to_sce(ice_flag)
	implicit none
	integer(i4b), intent(in) :: ice_flag
	
	
    ! map BMM to SCE
    parcel1%npartall(1:parcel1%n_bin_modew)=parcel1%npart
    parcel1%mbinall(1:parcel1%n_bin_modew,:)=parcel1%mbin
    parcel1%mbinall(1:parcel1%n_bin_modew,parcel1%n_comps+1)= &
        parcel1%y(1:parcel1%n_bin_modew)
    if(ice_flag.eq.1) then
        parcel1%npartall(1+parcel1%n_bin_modew:parcel1%n_bin_mode)= &
            parcel1%npartice           
        parcel1%mbinall(1+parcel1%n_bin_modew:parcel1%n_bin_mode,:)= &
            parcel1%mbinice
        parcel1%mbinall(1+parcel1%n_bin_modew:parcel1%n_bin_mode, &
                        parcel1%n_comps+1)= &
            parcel1%yice(1:parcel1%n_bin_modew)
    endif      
	
	
    end subroutine map_to_sce
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! map to bmm                                                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>map the SCE to the BMM variables
	!>@param[in] ice_flag - flag to say if we are computing ice
	subroutine map_to_bmm(ice_flag)
	implicit none
	integer(i4b), intent(in) :: ice_flag
	
	
    parcel1%npart=parcel1%npartall(1:parcel1%n_bin_modew)
    parcel1%mbin=parcel1%mbinall(1:parcel1%n_bin_modew,:)
    parcel1%y(1:parcel1%n_bin_modew)=parcel1%mbin(:,parcel1%n_comps+1)
    if(ice_flag.eq.1) then
        parcel1%npartice= &
            parcel1%npartall(1+parcel1%n_bin_modew:parcel1%n_bin_mode)         
        parcel1%mbinice= &
            parcel1%mbinall(1+parcel1%n_bin_modew:parcel1%n_bin_mode,:) 
        parcel1%yice(1:parcel1%n_bin_modew)=parcel1%mbinice(:,parcel1%n_comps+1)
                        
        
    endif 
	
	
    end subroutine map_to_bmm
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! adjust t and rh                                                              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>adjust the temperature and rh for latent heating
	!>@param[in] totmass,p
	!>@param[inout] t,rh
	subroutine adjust_t_rh(totmass,t,rh,p)
	implicit none
	real(wp), intent(inout) :: t,rh
	real(wp), intent(in) :: p,totmass
	
	real(wp) :: mr
	
	! total vapour m.r.
	mr=rh*eps1*svp_liq(t)/(p-svp_liq(t))
	! adjust t
	t=t+totmass*lf/cp
	! calculate rh
	rh = mr / (eps1*svp_liq(t)/(p-svp_liq(t)))

    end subroutine adjust_t_rh
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! driver for bmm                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>driver for the bin-microphysics module
	!>@param[in] sce_flag - flag to say if we are computing the SCE
    subroutine bmm_driver(sce_flag,hm_flag,break_flag,mode1_flag, mode2_flag)
    use numerics_type
    use sce, only : sce_microphysics, sce_sip_microphysics, qsmall
    implicit none
    integer(i4b), intent(in) :: sce_flag
    logical, intent(in) :: hm_flag, mode1_flag, mode2_flag
    integer(i4b), intent(in) :: break_flag
    integer(i4b) :: i, j, nt
    
    nt=ceiling(runtime / real(dt,kind=wp))
    do i=1,nt
        ! output to file
        call output(io1%new_file,outputfile)
        
        
        if ((updraft_type==2).and.(parcel1%TT>t_thresh)) parcel1%y(parcel1%iw)=0._wp
        ! one time-step of model
        call bin_microphysics(fparcelwarm, fparcelcold, & 
            icenucleation, noncollisional_iceformation)
        
        if(sce_flag.gt.0) then
            
            ! Map the BMM variables to the SCE variables
            call map_to_sce(ice_flag)
              
              
            ! one time-step of sce model
            if(sce_flag.eq.1) then
                call sce_microphysics(parcel1%n_bins1,parcel1%n_bin_mode,parcel1%n_comps+&
                                parcel1%imoms,&
                                parcel1%npartall,parcel1%moments,parcel1%momenttype, &
                                parcel1%ecoll,parcel1%indexc, &
                                parcel1%mbinall(:,n_comps+1),parcel1%dt, &
                                parcel1%y(parcel1%ite))
                ! latent heat of fusion
                if(ice_flag.eq.1) parcel1%yice(parcel1%itei)=parcel1%y(parcel1%ite)
            elseif(sce_flag.eq.2) then
                call sce_sip_microphysics(parcel1%n_bins1,parcel1%n_bin_mode,&
                                parcel1%n_comps+&
                                parcel1%imoms,&
                                parcel1%npartall,parcel1%moments,parcel1%momenttype, &
                                parcel1%ecoll,parcel1%indexc, &
                                parcel1%mbinall(:,n_comps+1),parcel1%vel,parcel1%dt, &
                                parcel1%y(parcel1%ite),parcel1%totaddto, &
                                mass_fragment1, mass_fragment2, mass_fragment3, &
                                hm_flag,break_flag,mode1_flag, mode2_flag )
                ! latent heat of fusion
                if(ice_flag.eq.1) then
                    call adjust_t_rh(parcel1%totaddto,parcel1%y(parcel1%ite), &
                                parcel1%y(parcel1%irh), parcel1%y(parcel1%ipr))
                    parcel1%yice(parcel1%irhi)=parcel1%y(parcel1%irh)
                    parcel1%yice(parcel1%itei)=parcel1%y(parcel1%ite)
                endif
            endif
                                        
            ! redefine the mass of each component of aerosol
            do j=1,parcel1%n_comps
                where (parcel1%npartall(:).gt.qsmall)
                    parcel1%mbinall(:,j)=parcel1%moments(:,j)/parcel1%npartall(:)
                end where
            enddo                             
                            
            ! map SCE to BMM
            call map_to_bmm(ice_flag)
            
        endif    
             
        ! break-out if flag has been set 
        if(parcel1%break_flag) exit
    enddo
    ! output to file
    call output(io1%new_file,outputfile)
    
    
    
    
    end subroutine bmm_driver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







	end module bmm	

