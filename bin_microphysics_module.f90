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
        						mass_fragment3=mass_fragment1, &
        						gam_fac_ent=1._wp/(1._wp+0.5_wp), & ! P+K, 12-25
        						gam_fac_ent2=1._wp+0.5_wp, &
        						onequarter=1._wp/4._wp, &
        						onethird=1._wp/3._wp, &
        						twothirds=2._wp/3._wp, fourthirds=4._wp/3._wp, &
                                dh2o=2.75e-10_wp, fhh_theta_min=1.e-12_wp, &
                                qsmall2=1.e-20_wp
		
		! the different bin schemes        						
		integer(i4b), parameter :: BIN_FULL_MOVING   = 0_i4b
		integer(i4b), parameter :: BIN_MOVING_CENTRE = 1_i4b
		integer(i4b), parameter :: BIN_CHEN_LAMB     = 2_i4b

        ! Conserved-moment semantics used by the collection solver.
        integer(i4b), parameter :: MOMENT_EXTENSIVE = 1_i4b
        integer(i4b), parameter :: MOMENT_NUMBER    = 2_i4b
        integer(i4b), parameter :: MOMENT_INHERIT   = 3_i4b

        ! Ice-nucleation mechanisms.  These are independent switches;
        ! component categories below determine heterogeneous eligibility.
        integer(i4b), parameter :: INUC_KOOP   = 1_i4b
        integer(i4b), parameter :: INUC_INAS   = 2_i4b
        integer(i4b), parameter :: INUC_DEMOTT = 3_i4b
        integer(i4b), parameter :: INUC_DAILY  = 4_i4b
        integer(i4b), parameter :: N_INUC_MECH = 4_i4b

        ! Component-level heterogeneous ice-nucleation categories.
        integer(i4b), parameter :: INP_NONE          = 0_i4b
        integer(i4b), parameter :: INP_DEMOTT        = 1_i4b
        integer(i4b), parameter :: INP_NIEMAND12     = 2_i4b
        integer(i4b), parameter :: INP_KAOLINITE_M11 = 3_i4b
        integer(i4b), parameter :: INP_KFELDSPAR_A13 = 4_i4b
        integer(i4b), parameter :: INP_DAILY25       = 5_i4b
        						
        ! l_inhom describes the mixing MODE used by the warm ODE on the
        ! current timestep.  l_inhom_event is true only when a discrete
        ! inhomogeneous entrainment event has actually been applied.
        logical :: l_inhom=.false., l_inhom_event=.false.

        type parcel
            ! variables for bin model
            integer(i4b) :: n_bins1,n_modes,n_comps, n_bin_mode, n_bin_modew, n_bin_mode1, &
                            n_sound, n_chamber, ice_flag, sce_flag,bin_scheme_flag, &
                            n_inp_classes, iinp_start, idemott
            real(wp) :: dt
			! Cumulative hydrometeor mass removed by fallout
			! kg hydrometeor / kg dry air
			real(wp) :: qfall_liq=0._wp
			real(wp) :: qfall_ice=0._wp
			
			! Hydrometeor mass removed during most recent timestep
			! kg hydrometeor / kg dry air
			real(wp) :: qfall_step_liq=0._wp
			real(wp) :: qfall_step_ice=0._wp
			
			! Cumulative particle number removed
			! particles / kg dry air
			real(wp) :: nfall_liq=0._wp
			real(wp) :: nfall_ice=0._wp

			! qchamber_bl is cumulative airborne water removed to the chamber
			! wall by the BL saturation cap. Positive values denote loss.
			! Fan and wall diagnostics are cumulative
            ! positive airborne hydrometeor/particle losses.
            real(wp) :: qchamber_bl=0._wp, qchamber_bl_step=0._wp
            real(wp) :: qfan_liq=0._wp, qfan_ice=0._wp
            real(wp) :: nfan_liq=0._wp, nfan_ice=0._wp
            real(wp) :: qwall_liq=0._wp, qwall_ice=0._wp
            real(wp) :: nwall_liq=0._wp, nwall_ice=0._wp
            real(wp), dimension(:), allocatable :: time_chamber, press_chamber, &
            										temp_chamber, dp_chamber, dt_chamber, &
            										qtot_chamber,dqtot_chamber
            real(wp), dimension(:,:), allocatable :: q_sound
            real(wp), dimension(:), allocatable :: t_sound, z_sound, rh_sound, &
                                                    p_sound, theta_q_sound
            real(wp) :: z,p,t,w,rh, rad,qinit, t_cbase, q_cbase, p_cbase, z_cbase, &
                        t_ctop, q_ctop, p_ctop, z_ctop, theta_q_cbase, theta_q_ctop, &
                        x_ent, theta_q, zlast, qtot
                        
                        
            ! liquid water
            real(wp), dimension(:), allocatable :: d, maer, npart, rho_core, &
                            rh_eq, rhoat, dw, da_dt, ndrop, npartall, npart_ent, &
                            npart_ent0, npart_temp
            real(wp), dimension(:,:), allocatable :: mbin, mbinall, rhobin, &
                                        nubin,molwbin,kappabin, &
                                        mbin_ent, mbin_ent0, &
                                        mbin_temp ! all bins x all comps
            ! variables for ODE:                    
            integer(i4b) :: neq, itol, ipr, ite, iz, iw, irh, ira, &
                            itask, istate, iopt, mf, lrw, liw
            integer(i4b), dimension(:), allocatable :: iwork
            integer(i4b), dimension(1) :: ipar
            real(wp) :: tt, tout
            real(wp), dimension(:), allocatable :: y, yold, atol, rwork
            real(wp), dimension(1) :: rpar
            real(wp), dimension(1) :: rtol
            
            ! ice water - dwice (Dmax); areaice (projected area normal to fall)
            real(wp), dimension(:), allocatable :: dice, maerice, npartice, rho_coreice, &
                            rh_eqice, rhoatice, dwice, areaice, da_dtice, nice, &
                            phi, rhoi, nump, rime
            real(wp), dimension(:,:), allocatable :: mbinice, rhobinice, &
                                        nubinice,molwbinice,kappabinice ! all bins x all comps  
                                        
                                        
            ! general
            integer(i4b) :: imoms
            real(wp), allocatable, dimension(:,:) :: moments, mbinedges,ecoll,ecoal, &
                                                    moments_ent, mbinedges_ent, &
                                                    moments_ent0, mbinedges_ent0, &
                                                    moments_temp, mbinedges_temp
            real(wp), allocatable, dimension(:) :: momtemp, vel, cd, nre
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
		logical :: fallout_flag=.false.
		logical  :: use_adt_optics=.false.
		real(wp) :: optics_wavelength=0.55e-6_wp		
		real(wp) :: residence_depth=100._wp        
        real(wp) :: zinit,tpert,winit,winit2, amplitude2, tau2, &
                    tinit,pinit,rhinit,radinit,z_ctop=-1._wp, alpha_therm, alpha_cond, &
                    alpha_therm_ice, alpha_dep, thresh_to_start_hom_mix

        ! Chamber forcing/options.  The measured time series remain in
        ! &chamber_spec; these switches say how they are used.
        logical :: chamber_force_pressure=.false., &
                   chamber_force_temperature=.false., &
                   chamber_force_qtot=.false.
        integer(i4b) :: chamber_bl_mix=0_i4b, chamber_bl_temp_mode=0_i4b, &
                        chamber_fan_loss=0_i4b, chamber_wall_loss=0_i4b
        real(wp) :: chamber_bl_tau=60._wp, chamber_bl_temp_offset=0._wp, &
                    chamber_fan_loss_kmax=0._wp, chamber_fan_loss_d50_ref=4.e-6_wp, &
                    chamber_fan_loss_exp=6._wp, chamber_fan_rpm=25000._wp, &
                    chamber_fan_rpm_ref=25000._wp, &
                    chamber_wall_ustar=0.02_wp, chamber_diameter=0._wp, &
                    chamber_height=0._wp

        ! Derived compatibility indicator: true whenever measured chamber
        ! P/T/qtot forcing is active.  It is no longer a namelist option and
        ! does not control the individual forcing terms.
        logical :: chamber_forcing_active=.false.

        integer(i4b) :: microphysics_flag=0, kappa_flag,updraft_type, vent_flag, &
                        sce_flag=0,ice_flag=0, bin_scheme_flag=1, entrain_period=0
        logical :: use_prof_for_tprh, hm_flag, mode1_flag, mode2_flag, &
        	bubble_flag, release_aerosol, entrain_aerosol
        integer(i4b) :: break_flag
        logical, dimension(N_INUC_MECH) :: ice_nucleation_mech = &
            [.true.,.false.,.true.,.false.]
        real(wp) :: dz,dt, runtime, t_thresh
        ! sounding spec
        real(wp) :: psurf, tsurf
        integer(i4b), parameter :: nlevels_r=1000
        integer(i4b), parameter :: nq=3
        integer(i4b) :: n_levels_s, n_levels_c=0_i4b, idum, n_sel
        real(wp) :: mult, rh_act
        real(wp), allocatable, dimension(:,:) :: q_read !nq x nlevels_r
        real(wp), allocatable, dimension(:) :: theta_read,rh_read,  z_read
        real(wp), allocatable, dimension(:) :: time_chamber, press_chamber, temp_chamber, &
        	wall_temp_chamber, qtot_chamber
        ! aerosol setup
        integer(i4b) :: n_intern, n_mode,n_sv,sv_flag,n_bins,n_comps, &
                        n_inp_classes=0_i4b
        ! aerosol_spec
        real(wp), allocatable, dimension(:,:) :: n_aer1,d_aer1,sig_aer1, mass_frac_aer1
        real(wp), allocatable, dimension(:) ::  molw_core1,density_core1,nu_core1, &
                                        kappa_core1, afhh_core1, bfhh_core1, inp_temp
        character(len=32), allocatable, dimension(:) :: inp_category
        integer(i4b), allocatable, dimension(:) :: inp_kind
        real(wp), allocatable, dimension(:) :: org_content1, molw_org1, kappa_org1, &
                                    density_org1, delta_h_vap1,nu_org1,log_c_star1
        
        

        ! variables for model
        real(wp) :: theta_surf,theta_init, &
            theta_q_sat,t1old, p111, w_cb, n_dummy, d_dummy, x2old=1.0_wp, &
            wvenv_send, tenv_send, dilute_send=1._wp
        logical :: set_theta_q_cb_flag=.true.
		integer(i4b) :: entrain_count=0
		real(wp) :: entrain_integral=0._wp

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
        real(wp), dimension(2), parameter :: gam_mu_cl=[260.163817050062335_wp, &
                                                    8.274747821396463_wp]


        character (len=200) :: outputfile='output', scefile='input'

	!private 
	!public :: read_in_bmm_namelist, initialise_bmm_arrays, bmm_driver, io1
    public
    private :: outputfile
	contains	
	
		
		
	! ============================================================================
	! allocate_arrays
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Allocates the sounding, chamber, aerosol-composition and volatility arrays used by BMM.
	!>@param[in] n_intern: number of internal lognormal modes in each aerosol mode
	!>@param[in] n_mode: number of external aerosol modes
	!>@param[in] n_sv: number of semivolatile/volatility bins
	!>@param[in] n_bins: number of aerosol size bins per mode
	!>@param[in] n_comps: number of aerosol composition components
	!>@param[in] nq: number of sounding water/species variables
	!>@param[in] n_levels_s: number of sounding levels
	!>@param[in] n_levels_c: number of chamber-data levels
	!>@param[inout] q_read,theta_read,rh_read,z_read: allocatable sounding arrays
	!>@param[inout] time_chamber,press_chamber,temp_chamber,wall_temp_chamber,qtot_chamber: allocatable chamber-forcing
	!>arrays
	!>@param[inout] n_aer1,d_aer1,sig_aer1: allocatable aerosol lognormal parameters
	!>@param[inout] mass_frac_aer1: allocatable component mass fractions for each aerosol mode
	!>@param[inout] molw_core1,density_core1,nu_core1,kappa_core1: allocatable aerosol-component
	!>thermodynamic properties
	!>@param[inout] org_content1,molw_org1,kappa_org1,density_org1: allocatable semivolatile-organic
	!>properties
	!>@param[inout] delta_h_vap1,nu_org1,log_c_star1: allocatable volatility-bin thermodynamic
	!>properties
	subroutine allocate_arrays(n_intern,n_mode,n_sv,n_bins,n_comps,nq,n_levels_s, &
		                    q_read,theta_read,rh_read,z_read, &
		                    n_levels_c, time_chamber, press_chamber, temp_chamber, &
		                    wall_temp_chamber, qtot_chamber, &
		                    n_aer1,d_aer1,sig_aer1,mass_frac_aer1, molw_core1, &
		                    density_core1, nu_core1, kappa_core1, &
		                    org_content1,molw_org1,kappa_org1,density_org1, &
		                    delta_h_vap1,nu_org1,log_c_star1)
		use numerics_type
		implicit none
		integer(i4b), intent(in) :: n_intern, n_mode, n_sv, n_bins,n_comps, nq, &
		                            n_levels_s, n_levels_c
		real(wp), dimension(:), allocatable, intent(inout) :: theta_read,rh_read,z_read, &
		                        org_content1,molw_org1,kappa_org1, &
		                        density_org1,delta_h_vap1,nu_org1,log_c_star1, &
		                        time_chamber, press_chamber, temp_chamber, wall_temp_chamber, qtot_chamber
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
		
		allocate( time_chamber(1:n_levels_c), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( press_chamber(1:n_levels_c), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( temp_chamber(1:n_levels_c), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( wall_temp_chamber(1:n_levels_c), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        wall_temp_chamber=-huge(1._wp)
		allocate( qtot_chamber(1:n_levels_c), STAT = AllocateStatus)
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
		allocate( afhh_core1(1:n_comps), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( bfhh_core1(1:n_comps), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( inp_temp(1:n_inp_classes), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( inp_category(1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( inp_kind(1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"

		! Backward-compatible defaults: no component uses FHH or INAS unless requested.
		afhh_core1=0._wp
		bfhh_core1=0._wp
        inp_temp=0._wp
        inp_category='demott'
        inp_kind=INP_DEMOTT

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
	
	
	
	! ============================================================================
	! read_in_bmm_namelist
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Reads the BMM run, sounding, chamber, aerosol and composition namelists and validates component
	!>properties, including FHH adsorption coefficients.
	!>@param[in] nmlfile: path to the BMM namelist file
	subroutine read_in_bmm_namelist(nmlfile)
		implicit none
        character(len=*), intent(in) :: nmlfile
        integer(i4b) :: i,j,k,ios
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        namelist /sounding_spec/ psurf, tsurf,  &
                    q_read, theta_read, rh_read, z_read
        namelist /chamber_spec/ time_chamber, press_chamber, temp_chamber, wall_temp_chamber, &
                    qtot_chamber
        namelist /chamber_options/ n_levels_c, &
                    chamber_force_pressure, chamber_force_temperature, chamber_force_qtot, &
                    chamber_bl_mix, chamber_bl_tau, chamber_bl_temp_mode, &
                    chamber_bl_temp_offset, &
                    chamber_fan_loss, chamber_fan_loss_kmax, chamber_fan_loss_d50_ref, &
                    chamber_fan_loss_exp, chamber_fan_rpm, chamber_fan_rpm_ref, &
                    chamber_wall_loss, chamber_wall_ustar, chamber_diameter, chamber_height
        ! define namelists for environment
        namelist /run_vars/ outputfile, scefile,runtime, dt, &
                    zinit,tpert,use_prof_for_tprh,winit,winit2,amplitude2, &
                    tinit,pinit,rhinit, radinit, bubble_flag, &
                    microphysics_flag, ice_flag, bin_scheme_flag, sce_flag, &
                    ice_nucleation_mech, &
                    hm_flag, break_flag, mode1_flag, mode2_flag, &
                    use_adt_optics, optics_wavelength, vent_flag, &
                    kappa_flag, updraft_type,t_thresh, adiabatic_prof, &
                    entrain_period, thresh_to_start_hom_mix, release_aerosol, &
                    entrain_aerosol, vert_ent, &
                    z_ctop, ent_rate,n_levels_s, &
                    fallout_flag,residence_depth, &
                    alpha_therm,alpha_cond,alpha_therm_ice,alpha_dep
        namelist /aerosol_setup/ n_intern,n_mode,n_sv,sv_flag, n_bins,n_comps,n_inp_classes
        namelist /aerosol_spec/ n_aer1,d_aer1,sig_aer1, dmina,dmaxa, &
                                mass_frac_aer1, molw_core1, &
                                density_core1,nu_core1,kappa_core1, &
                                afhh_core1,bfhh_core1, inp_temp,inp_category, &
                                org_content1, molw_org1,kappa_org1, &
                                density_org1, delta_h_vap1,nu_org1, log_c_star1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists	and allocate arrays								   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        rewind(8)
        read(8,nml=run_vars)

        ! Reset optional chamber controls on every namelist read so omission of
        ! &chamber_options is guaranteed to mean no chamber-specific forcing or
        ! process, even if this routine is called more than once in one process.
        n_levels_c=0_i4b
        chamber_force_pressure=.false.
        chamber_force_temperature=.false.
        chamber_force_qtot=.false.
        chamber_bl_mix=0_i4b
        chamber_bl_tau=60._wp
        chamber_bl_temp_mode=0_i4b
        chamber_bl_temp_offset=0._wp
        chamber_fan_loss=0_i4b
        chamber_fan_loss_kmax=0._wp
        chamber_fan_loss_d50_ref=4.e-6_wp
        chamber_fan_loss_exp=6._wp
        chamber_fan_rpm=25000._wp
        chamber_fan_rpm_ref=25000._wp
        chamber_wall_loss=0_i4b
        chamber_wall_ustar=0.02_wp
        chamber_diameter=0._wp
        chamber_height=0._wp

        ! Chamber options are deliberately optional.  With no &chamber_options
        ! group the defaults above leave every chamber-specific process off,
        ! which keeps future parcel/LES use independent of chamber machinery.
        rewind(8)
        read(8,nml=chamber_options,iostat=ios)
        if (ios.gt.0) error stop 'Malformed &chamber_options namelist'

        ! Derived indicator used internally; not a namelist control.
        chamber_forcing_active=chamber_force_pressure .or. chamber_force_temperature .or. &
                         chamber_force_qtot

        if (vert_ent) error stop &
            'vert_ent/Sanchez cloud-top entrainment has been removed; use lateral entrainment'
        if (n_levels_c.lt.0) error stop 'n_levels_c must be non-negative'
        if (chamber_forcing_active .and. n_levels_c.lt.2) error stop &
            'Measured chamber forcing requires n_levels_c >= 2'
        if (chamber_bl_mix.lt.0 .or. chamber_bl_mix.gt.2) error stop &
            'chamber_bl_mix must be 0 (off), 1 (homogeneous), or 2 (inhomogeneous)'
        if (chamber_bl_mix.gt.0 .and. chamber_bl_tau.le.0._wp) error stop &
            'chamber_bl_tau must be > 0 when chamber_bl_mix is enabled'
        if (chamber_bl_temp_mode.lt.0 .or. chamber_bl_temp_mode.gt.1) error stop &
            'chamber_bl_temp_mode must be 0 (fixed gas-temperature offset) or 1 (wall time series)'
        if (chamber_bl_mix.gt.0 .and. chamber_bl_temp_mode.eq.1 .and. n_levels_c.lt.2) error stop &
            'Wall-temperature BL forcing requires n_levels_c >= 2'
        if (chamber_fan_loss.lt.0 .or. chamber_fan_loss.gt.1) error stop &
            'chamber_fan_loss must be 0 (off) or 1 (sigmoid inertial collection)'
        if (chamber_fan_loss.eq.1) then
            if (chamber_fan_loss_kmax.lt.0._wp) error stop &
                'chamber_fan_loss_kmax must be non-negative'
            if (chamber_fan_loss_d50_ref.le.0._wp) error stop &
                'chamber_fan_loss_d50_ref must be > 0'
            if (chamber_fan_loss_exp.le.0._wp) error stop &
                'chamber_fan_loss_exp must be > 0'
            if (chamber_fan_rpm.le.0._wp .or. chamber_fan_rpm_ref.le.0._wp) error stop &
                'chamber_fan_rpm and chamber_fan_rpm_ref must be > 0'
        endif
        if (chamber_wall_loss.lt.0 .or. chamber_wall_loss.gt.1) error stop &
            'chamber_wall_loss must be 0 (off) or 1 (Lai-Nazaroff smooth-wall deposition)'
        if (chamber_wall_loss.eq.1) then
            if (chamber_wall_ustar.le.0._wp) error stop &
                'chamber_wall_ustar must be > 0 when wall deposition is enabled'
            if (chamber_diameter.le.0._wp .or. chamber_height.le.0._wp) error stop &
                'Positive chamber_diameter and chamber_height are required for wall deposition'
        endif
        if (chamber_force_qtot .and. (chamber_bl_mix.gt.0 .or. chamber_fan_loss.gt.0 .or. &
            chamber_wall_loss.gt.0)) then
            print *, 'Warning: measured qtot forcing plus chamber BL/fan/wall sinks can double-count water loss'
        endif

        rewind(8)
        read(8,nml=aerosol_setup)
        ! allocate memory / init
		call allocate_arrays(n_intern,n_mode,n_sv,n_bins,n_comps,nq,n_levels_s, &
		                    q_read,theta_read,rh_read,z_read, &
		                    n_levels_c, time_chamber, press_chamber, temp_chamber, &
		                    wall_temp_chamber, qtot_chamber, &
		                    n_aer1,d_aer1,sig_aer1,mass_frac_aer1, molw_core1, &
		                    density_core1, nu_core1, kappa_core1, &
		                    org_content1,molw_org1,kappa_org1,density_org1, &
		                    delta_h_vap1,nu_org1,log_c_star1)
        
        rewind(8)
        read(8,nml=sounding_spec)
        rewind(8)
        read(8,nml=aerosol_spec)

        ! Basic configuration checks.  Zero-number lognormal submodes may use
        ! zero placeholder diameters/widths for backward compatibility.
        if (n_intern.lt.1 .or. n_mode.lt.1 .or. n_bins.lt.1 .or. &
            n_comps.lt.1) error stop 'Aerosol dimensions must be positive'
        if (dmina.le.0._wp .or. dmaxa.le.dmina) error stop &
            'Require 0 < dmina < dmaxa'
        if (any(n_aer1.lt.0._wp)) error stop &
            'Aerosol number concentrations must be non-negative'
        do j=1,n_mode
            do i=1,n_intern
                if (n_aer1(i,j).gt.0._wp) then
                    if (d_aer1(i,j).le.0._wp) error stop &
                        'Populated aerosol submode requires positive median diameter'
                    if (sig_aer1(i,j).le.0._wp) error stop &
                        'Populated aerosol submode requires positive ln-sigma'
                endif
            enddo
            if (any(mass_frac_aer1(j,:).lt.0._wp)) error stop &
                'Aerosol component mass fractions must be non-negative'
            if (abs(sum(mass_frac_aer1(j,:))-1._wp).gt.1.e-6_wp) error stop &
                'Aerosol component mass fractions must sum to one in each mode'
        enddo

        ! A_FHH=0 means that a component has no adsorption contribution.
        ! If A_FHH>0, B_FHH must also be positive.  Components are otherwise
        ! treated like any other aerosol component; normally mineral dust has
        ! nu=0 and kappa=0, while soluble components have A_FHH=B_FHH=0.
        do k=1,n_comps
            if (afhh_core1(k).lt.0._wp .or. bfhh_core1(k).lt.0._wp) error stop &
                'A_FHH and B_FHH must be non-negative'
            if (afhh_core1(k).gt.0._wp .and. bfhh_core1(k).le.0._wp) error stop &
                'A component with A_FHH > 0 requires B_FHH > 0'
            if (density_core1(k).le.0._wp) error stop &
                'Aerosol component requires positive density'
            if (molw_core1(k).le.0._wp) error stop &
                'Aerosol component requires positive molecular weight'
        enddo

        ! Convert the namelist strings to compact integer categories once.
        do k=1,n_comps
            select case(trim(adjustl(inp_category(k))))
            case('none')
                inp_kind(k)=INP_NONE
            case('demott','demott2010')
                inp_kind(k)=INP_DEMOTT
            case('niemand','niemand12','desert_dust_niemand12')
                inp_kind(k)=INP_NIEMAND12
            case('kaolinite','kaolinite_murray11')
                inp_kind(k)=INP_KAOLINITE_M11
            case('kfeldspar','k-feldspar','kfeldspar_atkinson13')
                inp_kind(k)=INP_KFELDSPAR_A13
            case('daily25','daily2025','daly25','dcmex_daily25','dcmex_daly25')
                inp_kind(k)=INP_DAILY25
            case default
                print *,'Unknown inp_category(',k,'): ',trim(inp_category(k))
                error stop 'Unknown inp_category'
            end select
        enddo

        if (ice_nucleation_mech(INUC_INAS)) then
            if (n_inp_classes.lt.1) error stop &
                'INAS enabled but n_inp_classes < 1'
            do k=2,n_inp_classes
                if (inp_temp(k).ge.inp_temp(k-1)) error stop &
                    'inp_temp must be ordered from warm to cold'
            enddo
            if (.not.any(inp_kind.eq.INP_NIEMAND12 .or. &
                         inp_kind.eq.INP_KAOLINITE_M11 .or. &
                         inp_kind.eq.INP_KFELDSPAR_A13)) then
                print *,'Warning: INAS enabled but no component has an INAS category'
            endif
        endif

        if(n_levels_c.gt.0) then
            rewind(8)
            read(8,nml=chamber_spec)
            if (chamber_bl_mix.gt.0 .and. chamber_bl_temp_mode.eq.1) then
                if (any(wall_temp_chamber.lt.150._wp) .or. any(wall_temp_chamber.gt.400._wp)) &
                    error stop 'chamber_bl_temp_mode=1 requires valid wall_temp_chamber values [K]'
            endif
        endif
        if(chamber_force_pressure .or. chamber_force_temperature) then
            use_prof_for_tprh=.false.
        endif
        close(8)
        tau2 = 2._wp*pi/winit2*amplitude2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine read_in_bmm_namelist


    ! ============================================================================
    ! is_explicit_inas
    ! ============================================================================
    logical function is_explicit_inas(kind)
        implicit none
        integer(i4b), intent(in) :: kind

        is_explicit_inas = kind.eq.INP_NIEMAND12 .or. &
                           kind.eq.INP_KAOLINITE_M11 .or. &
                           kind.eq.INP_KFELDSPAR_A13
    end function is_explicit_inas


    ! ============================================================================
    ! ice_active_site_density
    ! ============================================================================
    ! Ice-active surface-site density, m-2 of geometric component area.
    ! The fits are evaluated only over their published temperature ranges.
    ! At warmer temperatures ns is zero; at colder temperatures the cold-end
    ! value is held fixed rather than extrapolating the exponential.
    function ice_active_site_density(tc,kind) result(ns_site)
        implicit none
        real(wp), intent(in) :: tc
        integer(i4b), intent(in) :: kind
        real(wp) :: ns_site, teval, tk

        ns_site=0._wp
        tk=tc+ttr

        select case(kind)
        case(INP_NIEMAND12)
            ! Niemand et al. (2012), natural desert dust.
            ! ns = exp(-0.517 Tc + 8.934) m-2; -12 >= Tc >= -36 C.
            if (tc.gt.-12._wp) return
            teval=max(tc,-36._wp)
            ns_site=exp(-0.517_wp*teval+8.934_wp)

        case(INP_KAOLINITE_M11)
            ! Murray et al. (2011), kaolinite singular fit.
            ! ln(ns_BET[cm-2]) = -0.8881 T[K] + 226.29.
            ! Convert cm-2 -> m-2 and BET -> geometric area using 3.49.
            if (tk.gt.245.5_wp) return
            teval=max(tk,236.1_wp)
            ns_site=3.49_wp*1.e4_wp*exp(-0.8881_wp*teval+226.29_wp)

        case(INP_KFELDSPAR_A13)
            ! Atkinson et al. (2013), K-feldspar.
            ! ln(ns_BET[cm-2]) = -1.038 T[K] + 275.26; 248--268 K.
            ! Convert cm-2 -> m-2 and BET -> geometric area using 3.5.
            if (tk.gt.268._wp) return
            teval=max(tk,248._wp)
            ns_site=3.5_wp*1.e4_wp*exp(-1.038_wp*teval+275.26_wp)

        case default
            ns_site=0._wp
        end select
    end function ice_active_site_density


    ! ============================================================================
    ! initialise_inp_moments
    ! ============================================================================
    ! Convert a fresh aerosol population into cumulative intrinsic freezing
    ! threshold moments.  This routine must be called only for newly created
    ! aerosol populations (initial aerosol, entrained aerosol, emissions).
    ! Once created, the IN moments are prognostic and must not be recalculated
    ! from composition after coagulation or phase changes.
    subroutine initialise_inp_moments(npart,mbin_comp,rhobin_comp,moments, &
                                      nbin,ncomps_local,iinp_start_local)
        implicit none
        integer(i4b), intent(in) :: nbin,ncomps_local,iinp_start_local
        real(wp), dimension(nbin), intent(in) :: npart
        real(wp), dimension(nbin,ncomps_local), intent(in) :: mbin_comp,rhobin_comp
        real(wp), dimension(:,:), intent(inout) :: moments

        integer(i4b) :: i,j,k
        real(wp) :: vcomp, acomp, lambda, ffrz, fprev, ns_site

        if (n_inp_classes.le.0) return

        moments(1:nbin,iinp_start_local: &
                iinp_start_local+n_inp_classes-1)=0._wp

        do i=1,nbin
            if (npart(i).le.qsmall2) cycle

            fprev=0._wp
            do k=1,n_inp_classes
                lambda=0._wp

                do j=1,ncomps_local
                    if (.not.is_explicit_inas(inp_kind(j))) cycle
                    if (mbin_comp(i,j).le.tiny(1._wp)) cycle
                    if (rhobin_comp(i,j).le.0._wp) cycle

                    ! Geometric area of the volume-equivalent sphere formed
                    ! by this component alone.  This is also consistent with
                    ! the component-area proxy used by the FHH treatment.
                    vcomp=mbin_comp(i,j)/rhobin_comp(i,j)
                    acomp=pi*(6._wp*vcomp/pi)**twothirds
                    ns_site=ice_active_site_density(inp_temp(k),inp_kind(j))
                    lambda=lambda+acomp*ns_site
                enddo

                ffrz=1._wp-exp(-max(lambda,0._wp))
                ffrz=max(fprev,min(ffrz,1._wp))

                moments(i,iinp_start_local+k-1)=npart(i)*ffrz
                fprev=ffrz
            enddo
        enddo
    end subroutine initialise_inp_moments


    ! ============================================================================
    ! get_inp_control
    ! ============================================================================
    ! Diagnose heterogeneous-nucleation control from the aerosol composition
    ! carried by a representative particle.  The runtime precedence among
    ! enabled heterogeneous mechanisms is INAS > DeMott > Daily/DCMEX.
    subroutine get_inp_control(mcomp,has_inas,has_demott,has_daily)
        implicit none
        real(wp), dimension(:), intent(in) :: mcomp
        logical, intent(out) :: has_inas,has_demott,has_daily
        integer(i4b) :: j

        has_inas=.false.
        has_demott=.false.
        has_daily=.false.

        do j=1,min(size(mcomp),n_comps)
            if (mcomp(j).le.tiny(1._wp)) cycle
            if (is_explicit_inas(inp_kind(j))) then
                has_inas=.true.
            elseif (inp_kind(j).eq.INP_DAILY25) then
                has_daily=.true.
            elseif (inp_kind(j).eq.INP_DEMOTT) then
                has_demott=.true.
            endif
        enddo
    end subroutine get_inp_control


    ! ============================================================================
    ! particle_is_activated
    ! ============================================================================
    !>Return true when a warm BMM particle lies beyond the maximum of its
    !>current Koehler/FHH equilibrium curve.  This is the same activation
    !>criterion used by the ndrop diagnostic: current water mass must exceed
    !>the critical water mass at the Koehler/FHH maximum.  For a fixed dry
    !>composition this is equivalent to wet diameter being greater than the
    !>critical wet diameter.
    !>
    !>The minimisation helper functions use the module scratch variables
    !>n_sel, rh_act and mult, and historically read parcel1%t.  Save/restore
    !>those values here and temporarily set parcel1%t to the current ODE
    !>temperature so the activation test is evaluated at the current state.
    logical function particle_is_activated(ibin,mwat_current,t_current) result(activated)
        use numerics, only : fmin
        implicit none
        integer(i4b), intent(in) :: ibin
        real(wp), intent(in) :: mwat_current,t_current
        integer(i4b) :: n_sel_save
        real(wp) :: rh_act_save,mult_save,t_save,nwcrit,mcrit

        activated=.false.
        if (ibin.lt.1 .or. ibin.gt.parcel1%n_bin_modew) return
        if (mwat_current.le.tiny(1._wp)) return

        n_sel_save=n_sel
        rh_act_save=rh_act
        mult_save=mult
        t_save=parcel1%t

        n_sel=ibin
        rh_act=0._wp
        mult=-1._wp
        parcel1%t=t_current

        select case(kappa_flag)
        case(0)
            nwcrit=fmin(1.e-50_wp,1.e1_wp,koehler02,1.e-30_wp)
        case(1)
            nwcrit=fmin(1.e-50_wp,1.e1_wp,kkoehler02,1.e-30_wp)
        case default
            parcel1%t=t_save
            mult=mult_save
            rh_act=rh_act_save
            n_sel=n_sel_save
            error stop 'Unknown kappa_flag in particle_is_activated'
        end select

        mcrit=max(nwcrit,0._wp)*molw_water
        activated=mwat_current.gt.mcrit

        parcel1%t=t_save
        mult=mult_save
        rh_act=rh_act_save
        n_sel=n_sel_save
    end function particle_is_activated



	! ============================================================================
	! initialise_bmm_arrays
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Initialises the BMM parcel state, environmental profiles, aerosol size/composition bins,
	!>equilibrium water contents and conserved moments.
	!>@param[in] psurf,tsurf: surface pressure and temperature
	!>@param[in] q_read,theta_read,rh_read,z_read: sounding water, potential-temperature,
	!>relative-humidity and height profiles
	!>@param[in] time_chamber,press_chamber,temp_chamber,qtot_chamber: chamber forcing data
	!>@param[in] runtime,dt: model run time and timestep
	!>@param[in] zinit,tpert,winit,tinit,pinit,rhinit,radinit: initial parcel height, perturbation,
	!>vertical velocity, temperature, pressure, relative humidity and radius
	!>@param[in] use_prof_for_tprh,bubble_flag: profile and parcel-interface switches; chamber forcing is controlled by chamber_force_*
	!>environmental/profile and parcel geometry setup
	!>@param[in] microphysics_flag,ice_flag,bin_scheme_flag,vent_flag,kappa_flag,updraft_type:
	!>microphysics configuration flags
	!>@param[in] adiabatic_prof,vert_ent,z_ctop,ent_rate: entrainment/profile controls
	!>@param[in] n_levels_s,n_levels_c: number of sounding and chamber levels
	!>@param[in] alpha_therm,alpha_cond,alpha_therm_ice,alpha_dep: heat/mass accommodation coefficients
	!>@param[in] n_intern,n_mode,n_sv,sv_flag,n_bins,n_comps: aerosol and bin dimensions/configuration
	!>@param[in] n_aer1,d_aer1,sig_aer1,dmina,dmaxa: aerosol lognormal parameters and dry-size limits
	!>@param[in] mass_frac_aer1: component mass fractions for each aerosol mode
	!>@param[in] molw_core1,density_core1,nu_core1,kappa_core1: aerosol-component molecular weights,
	!>densities, van't Hoff factors and kappa values
	!>@param[in] org_content1,molw_org1,kappa_org1,density_org1,delta_h_vap1,nu_org1,log_c_star1:
	!>semivolatile-organic properties
	!>@param[in] sce_flag: stochastic-collection-equation switch
    subroutine initialise_bmm_arrays(psurf, tsurf, q_read, theta_read, rh_read, z_read, &
    				time_chamber, press_chamber, temp_chamber, qtot_chamber, &
                    runtime, dt, zinit, tpert, use_prof_for_tprh, &
                    winit, tinit, pinit, &
                    rhinit, radinit, bubble_flag, &
                    microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                    kappa_flag, updraft_type, adiabatic_prof, vert_ent, z_ctop, &
                    ent_rate, n_levels_s, n_levels_c, &
                    alpha_therm, alpha_cond, alpha_therm_ice, &
                    alpha_dep, n_intern, n_mode, n_sv, sv_flag, n_bins, n_comps, &
                    n_aer1,d_aer1,sig_aer1,dmina,dmaxa,mass_frac_aer1,molw_core1, &
                    density_core1, nu_core1, kappa_core1, org_content1, molw_org1, &
                    kappa_org1, density_org1, delta_h_vap1,nu_org1, log_c_star1, &
                    sce_flag)
    use numerics_type
    use numerics, only : find_pos, poly_int, zeroin, fmin,vode_integrate

    implicit none
    logical, intent(in) :: use_prof_for_tprh, adiabatic_prof, vert_ent, bubble_flag
    integer(i4b), intent(in) :: microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                    kappa_flag, updraft_type, n_levels_s,n_levels_c, &
                    n_intern, n_mode, n_sv, &
                    sv_flag, n_bins, n_comps, sce_flag
    real(wp), intent(in) :: psurf, radinit, &
                    tsurf, runtime, dt, zinit, tpert, winit, tinit, &
                    pinit, rhinit, alpha_therm, alpha_cond, alpha_therm_ice, &
                    alpha_dep, dmina,dmaxa, z_ctop, ent_rate
    real(wp), dimension(1:n_levels_c), intent(in) :: time_chamber, press_chamber, &
    												temp_chamber, qtot_chamber
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
    parcel1%n_chamber=n_levels_c
    parcel1%n_bins1=n_bins    
    parcel1%n_modes=n_mode
    parcel1%n_comps=n_comps
    parcel1%n_inp_classes=n_inp_classes
    parcel1%iinp_start=n_comps+6
    parcel1%idemott=parcel1%iinp_start+n_inp_classes
    parcel1%n_bin_modew=n_bins*n_mode
    parcel1%n_bin_mode1=(n_bins+1)*n_mode
    parcel1%z=zinit
    parcel1%p=pinit
    parcel1%t=tinit
    parcel1%w=winit
    parcel1%rh=rhinit
    parcel1%rad=radinit
    parcel1%bin_scheme_flag=bin_scheme_flag
    parcel1%dt=dt
	! Fallout diagnostics
	parcel1%qfall_liq=0._wp
	parcel1%qfall_ice=0._wp
	
	parcel1%qfall_step_liq=0._wp
	parcel1%qfall_step_ice=0._wp
	
	parcel1%nfall_liq=0._wp
	parcel1%nfall_ice=0._wp
	! Chamber diagnostics
	parcel1%qchamber_bl=0._wp
	parcel1%qchamber_bl_step=0._wp
	parcel1%qfan_liq=0._wp
	parcel1%qfan_ice=0._wp
	parcel1%nfan_liq=0._wp
	parcel1%nfan_ice=0._wp
	parcel1%qwall_liq=0._wp
	parcel1%qwall_ice=0._wp
	parcel1%nwall_liq=0._wp
	parcel1%nwall_ice=0._wp
    

    parcel1%ice_flag=ice_flag
    parcel1%n_bin_mode=&
        parcel1%n_bins1*n_mode*(1+parcel1%ice_flag)     ! for all the liquid and ice    
    parcel1%imoms=ice_flag*(6+n_inp_classes)            ! phi, nmon, vol, rim, unf, cumulative Nin, n_demott
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! validate bin_scheme_flag					                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	select case(parcel1%bin_scheme_flag)
	case(BIN_FULL_MOVING,BIN_MOVING_CENTRE,BIN_CHEN_LAMB)
	
	case default
		error stop 'Unknown bin_scheme_flag'
	end select
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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
    allocate( parcel1%cd(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%nre(1:parcel1%n_bin_mode), STAT = AllocateStatus)
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

    allocate( parcel1%time_chamber(1:parcel1%n_chamber), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%press_chamber(1:parcel1%n_chamber), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%temp_chamber(1:parcel1%n_chamber), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%qtot_chamber(1:parcel1%n_chamber), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%dp_chamber(1:parcel1%n_chamber-1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%dt_chamber(1:parcel1%n_chamber-1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%dqtot_chamber(1:parcel1%n_chamber-1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! Initialise BMM's diagnostic arrays
	parcel1%ecoll=0._wp
	parcel1%vel=0._wp
	parcel1%nre=0._wp
	parcel1%cd=0._wp


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
        ! for this make sure it isnt zero
        call lognormal_n_between_limits(max(n_aer1(:,k),1.e-6_wp),d_aer1(:,k),sig_aer1(:,k), &
                                        n_intern,dmina,dmaxa, num)
        !print *,num
        ! set up variables for parcel model
        ntot=num
        number_per_bin=ntot/real(n_bins,wp)
!         ! make sure it is zero if needed
!         parcel1%npart(1+(k-1)*n_bins:(k)*n_bins)=min(number_per_bin, &
!         		sum(n_aer1(:,k))/real(n_bins,wp))
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
    ! aerosol number and dry mass                                                  !
    !                                                                              !
    ! Calculate the actual number and third diameter moment between the            !
    ! bin edges. The representative aerosol mass is then the mean aerosol          !
    ! mass in the bin. This conserves both aerosol number and dry aerosol mass.     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=1,parcel1%n_modes
        do j=1,parcel1%n_bins1
            ! index of diameter edges
            i=j+(k-1)*(n_bins+1)
            ! actual number represented by this bin
            call lognormal_n_between_limits( &
                n_aer1(:,k), &
                d_aer1(:,k), &
                sig_aer1(:,k), &
                n_intern, &
                parcel1%d(i), &
                parcel1%d(i+1), &
                num)

            parcel1%npart(j+(k-1)*n_bins)=num
            ! third diameter moment represented by this bin
            call lognormal_m3_between_limits( &
                n_aer1(:,k), &
                d_aer1(:,k), &
                sig_aer1(:,k), &
                n_intern, &
                parcel1%d(i), &
                parcel1%d(i+1), &
                var)
                
            if (num > tiny(1._wp)) then
                ! mean dry aerosol mass per particle
                parcel1%maer(j+(k-1)*n_bins)= &
                    pi/6._wp * parcel1%rho_core(k) * var/num
            else
                ! Empty bin: mass is irrelevant because npart=0.
                ! Retain a sensible positive representative mass to avoid
                ! problems in routines that evaluate every bin.
                parcel1%maer(j+(k-1)*n_bins)= &
                    pi/6._wp * &
                    (0.5_wp*(parcel1%d(i+1)+parcel1%d(i)))**3 * &
                    parcel1%rho_core(k)
            endif
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
        parcel1%rh=rhinit
        !parcel1%t=parcel1%t+tpert
        print *,'t,p,rh from sounding: ', parcel1%t, parcel1%p, parcel1%rh
    endif
    
    ! Initialise/copy chamber data.  The measured time series can be retained
    ! even when only a subset of P/T/qtot is used as forcing.
    if (n_levels_c.ge.2) then
        parcel1%time_chamber=time_chamber
        parcel1%press_chamber=press_chamber
        parcel1%temp_chamber=temp_chamber
        parcel1%qtot_chamber=qtot_chamber

        ! Chamber model time starts at t=0.  Do not use parcel1%TT here: the
        ! ODE clock is initialised later in this routine.
        iloc=find_pos(parcel1%time_chamber(1:n_levels_c),0._wp)
        iloc=min(n_levels_c-1,iloc)
        iloc=max(1,iloc)

        if (chamber_force_temperature) then
            call poly_int(parcel1%time_chamber(iloc:iloc+1), &
                parcel1%temp_chamber(iloc:iloc+1), &
                min(max(0._wp,parcel1%time_chamber(1)), &
                    parcel1%time_chamber(n_levels_c)),var,dummy)
            parcel1%t=var+tpert
        endif

        if (chamber_force_pressure) then
            call poly_int(parcel1%time_chamber(iloc:iloc+1), &
                parcel1%press_chamber(iloc:iloc+1), &
                min(max(0._wp,parcel1%time_chamber(1)), &
                    parcel1%time_chamber(n_levels_c)),var,dummy)
            parcel1%p=var
        endif

        if (chamber_force_qtot) then
            call poly_int(parcel1%time_chamber(iloc:iloc+1), &
                parcel1%qtot_chamber(iloc:iloc+1), &
                min(max(0._wp,parcel1%time_chamber(1)), &
                    parcel1%time_chamber(n_levels_c)),var,dummy)
            parcel1%qtot=var
            ! At chamber initialisation the supplied qtot is interpreted as
            ! vapour for setting RH, matching the historical chamber setup.
            parcel1%rh=parcel1%qtot / &
                (eps1*svp_liq(parcel1%t)/(parcel1%p-svp_liq(parcel1%t)))
        else
            parcel1%rh=rhinit
        endif

        ! Piecewise-linear forcing tendencies.  All are retained so the user
        ! can switch individual forcing variables without rebuilding the data.
        do i=1,n_levels_c-1
            if (parcel1%time_chamber(i+1).le.parcel1%time_chamber(i)) error stop &
                'time_chamber must be strictly increasing'
            parcel1%dp_chamber(i)= &
                (parcel1%press_chamber(i+1)-parcel1%press_chamber(i))/ &
                (parcel1%time_chamber(i+1)-parcel1%time_chamber(i))
            parcel1%dt_chamber(i)= &
                (parcel1%temp_chamber(i+1)-parcel1%temp_chamber(i))/ &
                (parcel1%time_chamber(i+1)-parcel1%time_chamber(i))
            parcel1%dqtot_chamber(i)= &
                (parcel1%qtot_chamber(i+1)-parcel1%qtot_chamber(i))/ &
                (parcel1%time_chamber(i+1)-parcel1%time_chamber(i))
        enddo

        if (chamber_forcing_active) then
            print *,'t,p,rh from chamber forcing: ',parcel1%t,parcel1%p,parcel1%rh
        endif
    endif
    
    parcel1%zlast=parcel1%z
    entrain_count=0
    entrain_integral=0._wp
    parcel1%qinit=parcel1%rh*eps1*svp_liq(parcel1%t)/(parcel1%p-svp_liq(parcel1%t))
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
    if(adiabatic_prof) then
        parcel1%neq=parcel1%n_bin_modew+5 ! p,t,rh,z,w
    else
        parcel1%neq=parcel1%n_bin_modew+6 ! p,t,rh,z,w,ra
    endif
    parcel1%tt=0._wp

    ! A positive entrain_period selects the discrete extreme-inhomogeneous
    ! pathway before thresh_to_start_hom_mix.  At and after the threshold the
    ! model switches permanently to continuous homogeneous P&K entrainment.
    ! entrain_period=0 therefore means homogeneous mixing from t=0.
    l_inhom=(.not.adiabatic_prof) .and. (entrain_period.gt.0) .and. &
        (parcel1%tt.lt.thresh_to_start_hom_mix)
    l_inhom_event=.false.
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
    ! extra variables for entrainment
    if(.not.adiabatic_prof) then
        parcel1%ira =parcel1%n_bin_modew+6 ! radius of bubble or jet
        parcel1%atol(parcel1%ira) =2.e-2_wp
        if(parcel1%ira .ne. parcel1%neq) stop "*** problem with array lengths ***"
    else
        if(parcel1%iw .ne. parcel1%neq) stop "*** problem with array lengths ***"
    endif
    
    parcel1%atol(parcel1%ipr)=10._wp
    parcel1%atol(parcel1%ite)=1.e-4_wp
    parcel1%atol(parcel1%irh)=1.e-8_wp
    parcel1%atol(parcel1%iz) =2.e-2_wp
    parcel1%atol(parcel1%iw) =2.e-2_wp
    
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
    if(.not.adiabatic_prof) then
        parcel1%y(parcel1%ira) =parcel1%rad   
        ! note, this allocates some space for the aerosol to be entrained into the parcel
        allocate( parcel1%mbinedges_ent(1:parcel1%n_bins1+1,1:parcel1%n_modes), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%mbinedges_ent0(1:parcel1%n_bins1+1,1:parcel1%n_modes), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( parcel1%npart_ent(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%npart_ent0(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( parcel1%mbin_ent(1:parcel1%n_bin_modew,1:n_comps+1), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%mbin_ent0(1:parcel1%n_bin_modew,1:n_comps+1), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
		allocate( parcel1%moments_ent(1:parcel1%n_bin_mode,1:n_comps+parcel1%imoms), &
			STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( parcel1%moments_ent0(1:parcel1%n_bin_mode,1:n_comps+parcel1%imoms), &
			STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    endif

    ! Temporary warm population used for aerosol residuals released by either
    ! atmospheric extreme-inhomogeneous mixing or chamber BL mode 2.
    if ((.not.adiabatic_prof) .or. chamber_bl_mix.eq.2) then
        allocate( parcel1%mbinedges_temp(1:parcel1%n_bins1+1,1:parcel1%n_modes), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( parcel1%npart_temp(1:parcel1%n_bin_modew), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( parcel1%mbin_temp(1:parcel1%n_bin_modew,1:n_comps+1), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        allocate( parcel1%moments_temp(1:parcel1%n_bin_mode,1:n_comps+parcel1%imoms), &
            STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    endif
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
		allocate( parcel1%areaice(1:parcel1%n_bin_modew), STAT = AllocateStatus)
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
		parcel1%dwice=0._wp
		parcel1%areaice=0._wp        
        
                
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
    parcel1%momenttype(1:parcel1%n_comps)=MOMENT_EXTENSIVE

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
			parcel1%moments(i,parcel1%n_comps+1)= parcel1%npart(i)
			parcel1%moments(i,parcel1%n_comps+2)= parcel1%npart(i)
		
            ! rim: mass
            parcel1%moments(i,parcel1%n_comps+4)=parcel1%npart(i)* &
                parcel1%mbin(i,parcel1%n_comps+1)
            ! unf: mass
            parcel1%moments(i,parcel1%n_comps+5)=parcel1%npart(i)* &
                parcel1%mbin(i,parcel1%n_comps+1)

        enddo        
        parcel1%momenttype(parcel1%n_comps+1:parcel1%n_comps+5)= &
            [MOMENT_NUMBER,MOMENT_NUMBER,MOMENT_EXTENSIVE, &
             MOMENT_EXTENSIVE,MOMENT_EXTENSIVE]
        if (n_inp_classes.gt.0) then
            parcel1%momenttype(parcel1%iinp_start: &
                parcel1%iinp_start+n_inp_classes-1)=MOMENT_INHERIT

            ! Fresh aerosol receives its intrinsic IN threshold spectrum once.
            call initialise_inp_moments(parcel1%npart, &
                parcel1%mbin(:,1:n_comps),parcel1%rhobin,parcel1%moments, &
                parcel1%n_bin_modew,n_comps,parcel1%iinp_start)
        endif
        ! Number of currently frozen DeMott-origin primary monomers.
        ! This is extensive through aggregation and is zero in fresh liquid.
        parcel1%momenttype(parcel1%idemott)=MOMENT_EXTENSIVE
        parcel1%moments(:,parcel1%idemott)=0._wp
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end subroutine initialise_bmm_arrays
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! ============================================================================
	! write_sce_grid_to_bmm
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Copies the fixed SCE mass-bin edges and collision-product mode map into the BMM parcel state.
	!>@param[in] n_bin_mode: total number of bins represented by the SCE state
	!>@param[in] n_bin_modew: number of liquid/warm bins
	!>@param[in] n_binst: number of fixed bins per external mode
	!>@param[in] n_modes: number of external aerosol modes
	!>@param[in] n_comps: number of aerosol composition components
	!>@param[in] indexc: receiving-mode lookup table for collision products
	!>@param[in] mbinedges: fixed particle-mass bin edges for each mode
	subroutine write_sce_grid_to_bmm( &
		n_bin_mode,n_bin_modew,n_binst,n_modes,n_comps, &
		indexc,mbinedges)
		implicit none
	
		integer(i4b), intent(in) :: n_bin_mode,n_bin_modew,n_binst,n_modes,n_comps
		integer(i4b), dimension(:,:), intent(in) :: indexc	
		real(wp), dimension(:,:), intent(in) :: mbinedges

		! Copy only fixed-grid information
		parcel1%indexc = indexc
		parcel1%mbinedges = mbinedges
		! If these are intended to share the same fixed grid:
		parcel1%mbinedges_ent = mbinedges	
	end subroutine write_sce_grid_to_bmm
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! project_initial_bmm_to_fixed_grid
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Conservatively remaps the initial BMM liquid population onto the fixed SCE mass grid and
	!>synchronises the solution-vector water masses.
	subroutine project_initial_bmm_to_fixed_grid()
		implicit none
	
		call moving_centre( &
			parcel1%n_bin_mode, &
			parcel1%n_bin_modew, &
			parcel1%n_bins1, &
			parcel1%n_modes, &
			parcel1%n_comps, &
			parcel1%n_comps+parcel1%imoms, &
			parcel1%npart, &
			parcel1%y(1:parcel1%n_bin_modew), &
			parcel1%moments(1:parcel1%n_bin_modew,:), &
			parcel1%mbin, &
			parcel1%mbinedges)
	
		! Keep solution vector consistent with remapped wet masses
		parcel1%y(1:parcel1%n_bin_modew) = parcel1%mbin(:,parcel1%n_comps+1)
	end subroutine project_initial_bmm_to_fixed_grid
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	! ============================================================================
	! write_sce_to_bmm
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Copies an SCE particle population, composition and conserved moments into the BMM state and
	!>initialises associated liquid, ice and entraining populations.
	!>@param[in] n_bin_mode: total number of liquid and ice SCE bins
	!>@param[in] n_bin_modew: number of liquid/warm bins
	!>@param[in] n_binst: number of bins per mode
	!>@param[in] n_mode: number of external modes
	!>@param[in] n_comps: number of aerosol composition components
	!>@param[in] n_moments: number of conserved moments per bin
	!>@param[in] ice_flag: switch indicating whether ice bins are present
	!>@param[in] npart: particle number concentrations on the SCE grid
	!>@param[in] moments: conserved SCE moments
	!>@param[in] mbin: per-particle component and water/ice masses
	!>@param[in] indexc: collision-product receiving-mode lookup table
	!>@param[in] mbinedges: fixed particle-mass bin edges
	!>@param[in] adiabatic_prof: true for an adiabatic profile; false when an entraining environmental
	!>aerosol population is required
	
    subroutine write_sce_to_bmm(n_bin_mode,n_bin_modew,n_binst,n_mode, n_comps, n_moments, &
                    ice_flag, &
                    npart, moments, mbin, indexc,mbinedges, &
                    adiabatic_prof)
    implicit none
    integer(i4b), intent(in) :: n_bin_mode, n_bin_modew, &
        n_binst, n_mode, n_comps, n_moments, ice_flag
    real(wp), dimension(n_bin_mode), intent(in) :: npart
    real(wp), dimension(n_bin_mode,n_moments), intent(in) :: moments
    real(wp), dimension(n_bin_mode,n_comps+1), intent(in) :: mbin
    real(wp), dimension(n_binst+1,n_mode), intent(in) :: mbinedges
    integer(i4b), dimension(n_bin_mode,n_bin_mode), intent(in) :: indexc
    logical, intent(in) :: adiabatic_prof
    integer(i4b) :: i,j


    
    parcel1%moments=moments
    parcel1%npart=npart(1:n_bin_modew)
    parcel1%npartall=npart
    parcel1%y(1:n_bin_modew)=mbin(1:n_bin_modew,n_comps+1)
    parcel1%maer=sum(mbin(1:n_bin_modew,1:n_comps),2)
    parcel1%mbin=mbin(1:n_bin_modew,1:n_comps+1)
    parcel1%mbinedges=mbinedges
    parcel1%mbinall(:,:)=mbin
    
    ! entraining aerosol
    if(.not.adiabatic_prof) then
        parcel1%npart_ent=npart(1:n_bin_modew)
        ! note, 1:n_comps are the masses of the aerosol componenets in each bin
        ! n_comps + 1 is the water mass
        parcel1%mbin_ent=mbin(1:n_bin_modew,1:n_comps+1)
        parcel1%mbinedges_ent=mbinedges
		parcel1%moments_ent=0._wp
		do j=1,parcel1%n_comps
			! above the aerosol
			do i=1,parcel1%n_bin_modew
				! aerosol moments
				parcel1%moments_ent(i,j)=parcel1%npart_ent(i)*parcel1%mbin_ent(i,j)
			enddo
		enddo
		! Ice auxiliary moments only exist when ice microphysics is enabled.
		if (parcel1%ice_flag.eq.1) then
            parcel1%moments_ent(1:parcel1%n_bin_modew,parcel1%n_comps+1)= &
                parcel1%npart_ent
            parcel1%moments_ent(1:parcel1%n_bin_modew,parcel1%n_comps+2)= &
                parcel1%npart_ent
			do i=1,parcel1%n_bin_modew
				! other moments: phi, nmon, vol, rim, unf
				! rim: mass
				parcel1%moments_ent(i,parcel1%n_comps+4)=parcel1%npart_ent(i)* &
					parcel1%mbin_ent(i,parcel1%n_comps+1)
				! unf: mass
				parcel1%moments_ent(i,parcel1%n_comps+5)=parcel1%npart_ent(i)* &
					parcel1%mbin_ent(i,parcel1%n_comps+1)
			enddo
            if (parcel1%n_inp_classes.gt.0) then
                call initialise_inp_moments(parcel1%npart_ent, &
                    parcel1%mbin_ent(:,1:parcel1%n_comps),parcel1%rhobin, &
                    parcel1%moments_ent,parcel1%n_bin_modew,parcel1%n_comps, &
                    parcel1%iinp_start)
            endif
		endif

        ! Keep an immutable copy of the environmental aerosol distribution.
        ! npart_ent/mbin_ent/moments_ent are working arrays which may be
        ! equilibrated and remapped each timestep before they are mixed into
        ! the parcel.  The *_ent0 arrays always retain the t=0 aerosol PSD and
        ! composition so every entrainment event samples the same environment.
        parcel1%npart_ent0=parcel1%npart_ent
        parcel1%mbin_ent0=parcel1%mbin_ent
        parcel1%moments_ent0=parcel1%moments_ent
        parcel1%mbinedges_ent0=parcel1%mbinedges_ent
    endif
    
    
    parcel1%indexc=indexc

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

	! ============================================================================
	! lognormal_n_between_limits
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Integrates the number concentration of a sum of lognormal aerosol modes between two dry diameters.
	!>@param[in] n_aer1,d_aer1,sig_aer1: number concentration, median diameter and logarithmic width
	!>of each internal mode
	!>@param[in] n_intern: number of internal lognormal modes
	!>@param[in] dmin,dmax: lower and upper dry-diameter limits
	!>@param[inout] num: integrated number concentration between the limits
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
    
    
	! ============================================================================
	! lognormal_m3_between_limits
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Integrates the third diameter moment of a sum of lognormal aerosol modes between two dry
	!>diameters.
	!>@param[in] n_aer1,d_aer1,sig_aer1: number concentration, median diameter and logarithmic width
	!>of each internal mode
	!>@param[in] n_intern: number of internal lognormal modes
	!>@param[in] dmin,dmax: lower and upper dry-diameter limits
	!>@param[out] m3: integrated third diameter moment between the limits
    subroutine lognormal_m3_between_limits(n_aer1,d_aer1,sig_aer1, &
                                           n_intern,dmin,dmax,m3)
        implicit none

        integer(i4b), intent(in) :: n_intern
        real(wp), intent(in) :: dmin,dmax
        real(wp), dimension(n_intern), intent(in) :: n_aer1,d_aer1,sig_aer1
        real(wp), intent(out) :: m3
        integer(i4b) :: i
        real(wp) :: zmin,zmax,sig

        m3=0._wp
        do i=1,n_intern
            if (n_aer1(i) <= 0._wp) cycle

            sig=sig_aer1(i)
            zmin=log(dmin/d_aer1(i))/sig
            zmax=log(dmax/d_aer1(i))/sig

            m3=m3 + n_aer1(i)*d_aer1(i)**3 * &
                exp(4.5_wp*sig**2) * &
                ( &
                0.5_wp*erfc(-(zmax-3._wp*sig)/sqrt(2._wp)) - &
                0.5_wp*erfc(-(zmin-3._wp*sig)/sqrt(2._wp)) &
                )
        enddo
    end subroutine lognormal_m3_between_limits
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    


	! ============================================================================
	! find_upper_diameter
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Returns the root-finder residual used to locate a size-bin upper diameter containing the target
	!>aerosol number concentration.
	!>@param[in] x: trial upper dry diameter
	!>@return find_upper_diameter: integrated number between the current lower edge and x minus the
	!>target bin number
    function find_upper_diameter(x)
        use numerics_type
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: find_upper_diameter, num
        
        call lognormal_n_between_limits(max(n_aer1(:,idum),1e-6_wp), &
        			d_aer1(:,idum),sig_aer1(:,idum), &
                                    n_intern,d_dummy,x, num)
        find_upper_diameter=num-n_dummy
        
    end function find_upper_diameter
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
	! ============================================================================
	! cloud_base
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Returns the saturation-mixing-ratio residual used to locate parcel cloud-base pressure.
	!>@param[in] p: trial pressure
	!>@return cloud_base: saturation water mixing ratio at the trial pressure minus the parcel
	!>cloud-base total-water mixing ratio
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



	! ============================================================================
	! hydrostatic1
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates dz/dp for a dry hydrostatic profile using the prescribed surface potential temperature.
	!>@param[in] p: pressure
	!>@param[in] z: integration-state array containing height
	!>@param[out] dzdp: hydrostatic height derivative with respect to pressure
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
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! hydrostatic1a
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates dp/dz for a dry hydrostatic profile using the prescribed surface potential temperature.
	!>@param[in] z: height
	!>@param[in] p: integration-state array containing pressure
	!>@param[out] dpdz: hydrostatic pressure derivative with respect to height
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
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! hydrostatic1b
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates dp/dz hydrostatically using potential temperature interpolated from the input sounding.
	!>@param[in] z: height at which the derivative is required
	!>@param[in] p: integration-state array containing pressure
	!>@param[out] dpdz: hydrostatic pressure derivative with respect to height
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
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! hydrostatic2
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates dz/dp for a saturated/moist hydrostatic profile, solving temperature from the
	!>moist-potential-temperature constraint.
	!>@param[in] p: pressure
	!>@param[in] z: integration-state array containing height
	!>@param[out] dzdp: hydrostatic height derivative with respect to pressure
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
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! ============================================================================
	! hydrostatic2a
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates dp/dz for a saturated/moist hydrostatic profile, solving temperature from the
	!>moist-potential-temperature constraint.
	!>@param[in] z: height
	!>@param[in] p: integration-state array containing pressure
	!>@param[out] dpdz: hydrostatic pressure derivative with respect to height
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
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

	! ============================================================================
	! calc_theta_q
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Returns the moist-potential-temperature residual used by a root finder to recover saturated
	!>temperature at the module pressure p111.
	!>@param[in] t111: trial temperature
	!>@return calc_theta_q: calculated saturated moist potential temperature minus theta_q_sat
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
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! ============================================================================
	! calc_theta_q2
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Returns the moist-potential-temperature residual used by a root finder to recover pressure at
	!>the module temperature t1old.
	!>@param[in] p: trial pressure
	!>@return calc_theta_q2: calculated saturated moist potential temperature minus theta_q_sat
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


	! ============================================================================
	! calc_theta_q3
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates the parcel moist potential temperature from temperature, pressure and total-water
	!>mixing ratio.
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[in] q: total-water mixing ratio
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


	! ============================================================================
	! svp_liq
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates saturation vapour pressure over liquid water using the Buck formulation.
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
	
	
	! ============================================================================
	! svp_ice
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates saturation vapour pressure over ice using the Buck formulation.
	!>@param[in] t: temperature
	!>@return svp_ice: saturation vapour pressure over ice
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


	! ============================================================================
	! dd
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates the molecular diffusivity of water vapour in air.
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@return dd: water-vapour diffusivity in air
    function dd(t,p)
      use numerics_type
      implicit none
      real(wp), intent(in) :: t, p
      real(wp) :: dd, t1
      t1=max(t,200._wp)
      dd=2.11e-5_wp*(t1/ttr)**1.94_wp*(101325._wp/p)
    end function dd
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! ============================================================================
	! ka
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates the thermal conductivity of air.
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

	! ============================================================================
	! viscosity_air
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates the dynamic viscosity of air - page 417 pruppacher and klett.
	!>@param[in] t: temperature
	!>@return viscosity_air: dynamic viscosity of air
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
    

	! ============================================================================
	! surface_tension
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates the surface tension of liquid water - pruppacher and klett	.
	!>@param[in] t: temperature
	!>@return surface_tension: surface tension of water
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
    
    
    
	! ============================================================================
	! wetdiam
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates wet volume-equivalent particle diameter and, optionally, mean particle density from
	!>component and water masses.
	!>@param[in] mwat: water mass per representative particle
	!>@param[in] mbin: aerosol-component masses per representative particle
	!>@param[in] rhobin: densities of the aerosol components
	!>@param[in] sz: number of particles/bins to process
	!>@param[out] dw: wet volume-equivalent diameter
	!>@param[out] rhoat: optional mean wet-particle density
 	subroutine wetdiam(mwat,mbin,rhobin,sz,dw,rhoat)
		implicit none
		real(wp), dimension(:), intent(in) :: mwat
		real(wp), dimension(:,:), intent(in) :: mbin,rhobin
		integer(i4b), intent(in) :: sz
		real(wp), dimension(:), intent(out) :: dw
		real(wp), dimension(:), intent(out), optional :: rhoat	
		real(wp), dimension(sz) :: rhoat1

        ! calculate the diameter and radius
		rhoat1(:)=mwat(:)/rhow + &
			sum(mbin(:,1:n_comps)/rhobin(:,:),2)
		rhoat1(:)=(mwat(:)+sum(mbin(:,1:n_comps),2))/rhoat1(:)
		dw(:)=((mwat(:)+sum(mbin(:,1:n_comps),2))*6._wp / &
			(pi*rhoat1(:)))**onethird
	
        ! wet diameter:
		if (present(rhoat)) rhoat(:)=rhoat1(:)
		
	end subroutine wetdiam
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
	! ============================================================================
	! fhh_adsorption_factor
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates the FHH adsorption activity factor for a mixed particle using surface-area-weighted
	!> effective adsorption coefficients from components with A_FHH greater than zero.
	!> Returns the adsorption contribution to water activity.  The insoluble FHH
    !> component is treated as a spherical core.  Any soluble dry material plus water
    !> lies outside that core.  When no FHH component (or no FHH mass) is present the
    !> factor is exactly unity, so the existing Koehler equations are recovered.
	!>@param[in] mwat: water mass per representative particle
	!>@param[in] dw: wet particle diameter
	!>@param[in] mbin1: component masses per representative particle
	!>@param[in] rhobin1: component densities
	!>@return fads: multiplicative FHH adsorption factor; unity when no adsorbing component is present
    function fhh_adsorption_factor(mwat,dw,mbin1,rhobin1) result(fads)
        implicit none
        real(wp), intent(in) :: mwat,dw
        real(wp), dimension(:), intent(in) :: mbin1,rhobin1
        integer(i4b) :: j
        real(wp) :: fads,vads,vj,voutside,dcore,delta_d,theta,expo,denom
        real(wp) :: area_j,area_tot,afhh_eff,bfhh_eff

        fads=1._wp
        vads=0._wp
        area_tot=0._wp

        ! A_FHH=0 means that this component is non-adsorbing.
        ! For each adsorbing component, use the surface area of its
        ! volume-equivalent sphere as the weighting proxy: S_j ~ V_j^(2/3).
        do j=1,min(size(mbin1),n_comps)
            if (afhh_core1(j).le.0._wp) cycle
            if (mbin1(j).le.tiny(1._wp)) cycle
            if (rhobin1(j).le.0._wp) cycle

            vj=mbin1(j)/rhobin1(j)
            if (vj.le.tiny(1._wp)) cycle

            vads=vads+vj
            area_j=vj**twothirds
            area_tot=area_tot+area_j
        enddo

        ! No adsorbing material -> exactly the existing Koehler/K-kohler result.
        if (vads.le.tiny(1._wp) .or. area_tot.le.tiny(1._wp)) return

        ! Surface-area-weighted effective adsorption coefficients.
        afhh_eff=0._wp
        bfhh_eff=0._wp
        do j=1,min(size(mbin1),n_comps)
            if (afhh_core1(j).le.0._wp) cycle
            if (mbin1(j).le.tiny(1._wp)) cycle
            if (rhobin1(j).le.0._wp) cycle

            vj=mbin1(j)/rhobin1(j)
            if (vj.le.tiny(1._wp)) cycle

            area_j=vj**twothirds
            afhh_eff=afhh_eff+(area_j/area_tot)*afhh_core1(j)
            bfhh_eff=bfhh_eff+(area_j/area_tot)*bfhh_core1(j)
        enddo

        ! One equivalent adsorbing core containing all material with A_FHH>0.
        dcore=(6._wp*vads/pi)**onethird

        ! Volume outside that core = non-adsorbing dry material + water.
        voutside=max(max(mwat,0._wp)/rhow + sum(mbin1/rhobin1) - vads,0._wp)

        ! Evaluate Dw-Dcore without subtracting two nearly equal diameters.
        denom=dw*dw+dw*dcore+dcore*dcore
        if (denom.gt.tiny(1._wp)) then
            delta_d=(6._wp*voutside/pi)/denom
        else
            delta_d=0._wp
        endif

        theta=max(delta_d/(2._wp*dh2o),fhh_theta_min)
        expo=-afhh_eff*theta**(-bfhh_eff)

        ! Avoid signalling underflow for an essentially dry adsorbing surface.
        if (expo.le.log(tiny(1._wp))) then
            fads=0._wp
        else
            fads=exp(expo)
        endif
    end function fhh_adsorption_factor
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! koehler01
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates equilibrium relative humidity for particle arrays using classical Kohler solute
	!>activity with the mixed FHH adsorption correction.
	!>@param[in] t: temperature
	!>@param[in] mwat: water mass per representative particle
	!>@param[in] mbin: aerosol-component masses per representative particle
	!>@param[in] rhobin: component densities
	!>@param[in] nubin: van't Hoff factors for the components
	!>@param[in] molwbin: component molecular weights
	!>@param[in] sz: number of particles/bins to process
	!>@param[inout] rh_eq: equilibrium relative humidity
	!>@param[inout] rhoat: mean wet-particle density
	!>@param[inout] dw: wet particle diameter
    subroutine koehler01(t,mwat,mbin,rhobin,nubin,molwbin,sz,rh_eq,rhoat,dw) 
      use numerics_type
      implicit none
      real(wp), dimension(:), intent(in) :: mwat
      real(wp), dimension(:,:), intent(in) :: mbin,rhobin,nubin,molwbin
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz) :: nw,nsolute,aw,fads
      real(wp), dimension(:),intent(inout) :: rh_eq,rhoat, dw
      real(wp), intent(in) :: t
      real(wp) :: sigma
      integer(i4b) :: i

      nw(:)=mwat(:)/molw_water
      call wetdiam(mwat,mbin,rhobin,sz,dw,rhoat)

      nsolute(:)=sum(mbin(:,1:n_comps)/molwbin(:,:)*nubin(:,:),2)
      do i=1,sz
          if (nsolute(i).le.tiny(1._wp)) then
              aw(i)=1._wp
          else
              aw(i)=nw(i)/(nw(i)+nsolute(i))
          endif
          fads(i)=fhh_adsorption_factor(mwat(i),dw(i), &
              mbin(i,1:n_comps),rhobin(i,1:n_comps))
      enddo

      sigma=surface_tension(t)

      rh_eq(:)=exp(4._wp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))*aw(:)*fads(:)

    end subroutine koehler01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

	! ============================================================================
	! kkoehler01
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates equilibrium relative humidity for particle arrays using volume-mixed kappa-Kohler
	!>activity with the mixed FHH adsorption correction.
	!>@param[in] t: temperature
	!>@param[in] mwat: water mass per representative particle
	!>@param[in] mbin: aerosol-component masses per representative particle
	!>@param[in] rhobin: component densities
	!>@param[in] kappabin: component kappa values
	!>@param[in] molwbin: component molecular weights; retained for interface consistency
	!>@param[in] sz: number of particles/bins to process
	!>@param[inout] rh_eq: equilibrium relative humidity
	!>@param[inout] rhoat: mean wet-particle density
	!>@param[inout] dw: wet particle diameter
    subroutine kkoehler01(t,mwat,mbin,rhobin,kappabin,molwbin,sz,rh_eq,rhoat,dw) 
      use numerics_type
      implicit none
      real(wp), dimension(:), intent(in) :: mwat
      real(wp), dimension(:,:), intent(in) :: mbin,rhobin,kappabin,molwbin
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz) :: vdry,vwat,kappa,aw,fads
      real(wp), dimension(:),intent(inout) :: rh_eq,rhoat, dw
      real(wp), intent(in) :: t
      real(wp) :: sigma
      integer(i4b) :: i

      call wetdiam(mwat,mbin,rhobin,sz,dw,rhoat)

      do i=1,sz
          vdry(i)=sum(mbin(i,1:n_comps)/rhobin(i,1:n_comps))
          vwat(i)=max(mwat(i),0._wp)/rhow
          if (vdry(i).gt.tiny(1._wp)) then
              kappa(i)=sum(mbin(i,1:n_comps)/rhobin(i,1:n_comps)* &
                  kappabin(i,1:n_comps))/vdry(i)
          else
              kappa(i)=0._wp
          endif

          ! Algebraically identical to the existing P&K expression, but stable
          ! for a pure insoluble particle with kappa=0.
          if (kappa(i).le.tiny(1._wp)) then
              aw(i)=1._wp
          else
              aw(i)=vwat(i)/(vwat(i)+kappa(i)*vdry(i))
          endif

          fads(i)=fhh_adsorption_factor(mwat(i),dw(i), &
              mbin(i,1:n_comps),rhobin(i,1:n_comps))
      enddo

      sigma=surface_tension(t)

      rh_eq(:)=exp(4._wp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))*aw(:)*fads(:)
    end subroutine kkoehler01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	! ============================================================================
	! koehler02
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates the classical-Kohler plus FHH equilibrium-humidity residual for the currently selected
	!>BMM particle; intended for minimisation/root finding in water content.
	!>@param[in] nw: trial number of moles of particle water
	!>@return koehler02: signed equilibrium-humidity residual relative to rh_act, including the module
	!>multiplier mult
    function koehler02(nw) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: nw
      real(wp) :: massw
      real(wp) :: rhoat,dw,koehler02
      real(wp) :: sigma,maer1,nsolute,aw,fads

      maer1=sum(parcel1%mbin(n_sel,1:n_comps))
      massw=nw*molw_water
      rhoat=massw/rhow+sum(parcel1%mbin(n_sel,1:n_comps)/ &
            parcel1%rhobin(n_sel,1:n_comps))
      rhoat=(massw+maer1)/rhoat
      dw=((massw+maer1)*6._wp/(pi*rhoat))**onethird

      sigma=surface_tension(parcel1%t)
      nsolute=sum(parcel1%mbin(n_sel,1:n_comps)/ &
          parcel1%molwbin(n_sel,1:n_comps)*parcel1%nubin(n_sel,1:n_comps))
      if (nsolute.le.tiny(1._wp)) then
          aw=1._wp
      else
          aw=nw/(nw+nsolute)
      endif
      fads=fhh_adsorption_factor(massw,dw, &
          parcel1%mbin(n_sel,1:n_comps),parcel1%rhobin(n_sel,1:n_comps))

      koehler02=mult*(exp(4._wp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           aw*fads)-rh_act

    end function koehler02
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	! ============================================================================
	! kkoehler02
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates the kappa-Kohler plus FHH equilibrium-humidity residual for the currently selected BMM
	!>particle; intended for minimisation/root finding in water content.
	!>@param[in] nw: trial number of moles of particle water
	!>@return kkoehler02: signed equilibrium-humidity residual relative to rh_act, including the
	!>module multiplier mult
    function kkoehler02(nw) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: nw
      real(wp) :: massw
      real(wp) :: rhoat,dw,kappa,kkoehler02
      real(wp) :: sigma,maer1,vdry,vwat,aw,fads

      maer1=sum(parcel1%mbin(n_sel,1:n_comps))
      massw=nw*molw_water
      rhoat=massw/rhow+sum(parcel1%mbin(n_sel,1:n_comps)/ &
            parcel1%rhobin(n_sel,1:n_comps))
      rhoat=(massw+maer1)/rhoat
      dw=((massw+maer1)*6._wp/(pi*rhoat))**onethird

      sigma=surface_tension(parcel1%t)

      vdry=sum(parcel1%mbin(n_sel,1:n_comps)/parcel1%rhobin(n_sel,1:n_comps))
      vwat=max(massw,0._wp)/rhow
      if (vdry.gt.tiny(1._wp)) then
          kappa=sum(parcel1%mbin(n_sel,1:n_comps)/parcel1%rhobin(n_sel,1:n_comps)* &
               parcel1%kappabin(n_sel,1:n_comps))/vdry
      else
          kappa=0._wp
      endif
      if (kappa.le.tiny(1._wp)) then
          aw=1._wp
      else
          aw=vwat/(vwat+kappa*vdry)
      endif

      fads=fhh_adsorption_factor(massw,dw, &
          parcel1%mbin(n_sel,1:n_comps),parcel1%rhobin(n_sel,1:n_comps))

      kkoehler02=mult*(exp(4._wp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           aw*fads)-rh_act
    end function kkoehler02
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
	! ============================================================================
	! koehler02_ent
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates the classical-Kohler plus FHH equilibrium-humidity residual for the currently selected
	!>entraining aerosol particle.
	!>@param[in] nw: trial number of moles of particle water
	!>@return koehler02_ent: signed equilibrium-humidity residual relative to rh_act, including the
	!>module multiplier mult
    function koehler02_ent(nw) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: nw
      real(wp) :: massw
      real(wp) :: rhoat,dw,koehler02_ent
      real(wp) :: sigma,maer1,nsolute,aw,fads

      maer1=sum(parcel1%mbin_ent(n_sel,1:n_comps))
      massw=nw*molw_water
      rhoat=massw/rhow+sum(parcel1%mbin_ent(n_sel,1:n_comps)/ &
            parcel1%rhobin(n_sel,1:n_comps))
      rhoat=(massw+maer1)/rhoat
      dw=((massw+maer1)*6._wp/(pi*rhoat))**onethird

      sigma=surface_tension(tenv_send)
      nsolute=sum(parcel1%mbin_ent(n_sel,1:n_comps)/ &
          parcel1%molwbin(n_sel,1:n_comps)*parcel1%nubin(n_sel,1:n_comps))
      if (nsolute.le.tiny(1._wp)) then
          aw=1._wp
      else
          aw=nw/(nw+nsolute)
      endif
      fads=fhh_adsorption_factor(massw,dw, &
          parcel1%mbin_ent(n_sel,1:n_comps),parcel1%rhobin(n_sel,1:n_comps))

      koehler02_ent=mult*(exp(4._wp*molw_water*sigma/r_gas/tenv_send/rhoat/dw)* &
           aw*fads)-rh_act

    end function koehler02_ent
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	! ============================================================================
	! kkoehler02_ent
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates the kappa-Kohler plus FHH equilibrium-humidity residual for the currently selected
	!>entraining aerosol particle.
	!>@param[in] nw: trial number of moles of particle water
	!>@return kkoehler02_ent: signed equilibrium-humidity residual relative to rh_act, including the
	!>module multiplier mult
    function kkoehler02_ent(nw) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: nw
      real(wp) :: massw
      real(wp) :: rhoat,dw,kappa,kkoehler02_ent
      real(wp) :: sigma,maer1,vdry,vwat,aw,fads

      maer1=sum(parcel1%mbin_ent(n_sel,1:n_comps))
      massw=nw*molw_water
      rhoat=massw/rhow+sum(parcel1%mbin_ent(n_sel,1:n_comps)/ &
            parcel1%rhobin(n_sel,1:n_comps))
      rhoat=(massw+maer1)/rhoat
      dw=((massw+maer1)*6._wp/(pi*rhoat))**onethird

      sigma=surface_tension(tenv_send)

      vdry=sum(parcel1%mbin_ent(n_sel,1:n_comps)/parcel1%rhobin(n_sel,1:n_comps))
      vwat=max(massw,0._wp)/rhow
      if (vdry.gt.tiny(1._wp)) then
          kappa=sum(parcel1%mbin_ent(n_sel,1:n_comps)/parcel1%rhobin(n_sel,1:n_comps)* &
               parcel1%kappabin(n_sel,1:n_comps))/vdry
      else
          kappa=0._wp
      endif
      if (kappa.le.tiny(1._wp)) then
          aw=1._wp
      else
          aw=vwat/(vwat+kappa*vdry)
      endif

      fads=fhh_adsorption_factor(massw,dw, &
          parcel1%mbin_ent(n_sel,1:n_comps),parcel1%rhobin(n_sel,1:n_comps))

      kkoehler02_ent=mult*(exp(4._wp*molw_water*sigma/r_gas/tenv_send/rhoat/dw)* &
           aw*fads)-rh_act
    end function kkoehler02_ent
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
	! ============================================================================
	! koehler03
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates the classical-Kohler plus FHH residual as a function of total dry aerosol mass for the
	!>selected prescribed aerosol mode.
	!> n_sel must be set to the mode, d_dummy the water in the bin to return the 
	!> aerosol dry mass
	!>@param[in] mbin: trial total dry aerosol mass per particle
	!>@return koehler03: signed equilibrium-humidity residual relative to rh_act, including the module
	!>multiplier mult
    function koehler03(mbin) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: mbin ! mass of aerosol particle
      real(wp) :: rhoat,dw,nw,koehler03
      real(wp) :: sigma,nsolute,aw,fads
      real(wp), dimension(n_comps) :: mbin_comp

      nw=d_dummy/molw_water
      rhoat=d_dummy/rhow+mbin*sum(mass_frac_aer1(n_sel,1:n_comps)/ &
                           density_core1(1:n_comps))
      rhoat=(d_dummy+mbin)/rhoat
      dw=((d_dummy+mbin)*6._wp/(pi*rhoat))**onethird

      sigma=surface_tension(parcel1%t)
      mbin_comp=mbin*mass_frac_aer1(n_sel,1:n_comps)
      nsolute=mbin*sum(mass_frac_aer1(n_sel,1:n_comps)/ &
          molw_core1(1:n_comps)*nu_core1(1:n_comps))
      if (nsolute.le.tiny(1._wp)) then
          aw=1._wp
      else
          aw=nw/(nw+nsolute)
      endif
      fads=fhh_adsorption_factor(d_dummy,dw,mbin_comp,density_core1(1:n_comps))

      koehler03=mult*(exp(4._wp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           aw*fads)-rh_act

    end function koehler03
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    

	! ============================================================================
	! kkoehler03
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates the kappa-Kohler plus FHH residual as a function of total dry aerosol mass for the
	!>selected prescribed aerosol mode.
	!> n_sel must be set to the mode, d_dummy the water in the bin to return the 
	!> aerosol dry mass
	!>@param[in] mbin: trial total dry aerosol mass per particle
	!>@return kkoehler03: signed equilibrium-humidity residual relative to rh_act, including the
	!>module multiplier mult
    function kkoehler03(mbin) ! only pass one variable so can use root-finders
      use numerics_type
      implicit none
      real(wp), intent(in) :: mbin ! mass of aerosol particle
      real(wp) :: rhoat,dw,kappa,kkoehler03
      real(wp) :: sigma,vdry,vwat,aw,fads
      real(wp), dimension(n_comps) :: mbin_comp

      rhoat=d_dummy/rhow+mbin*sum(mass_frac_aer1(n_sel,1:n_comps)/ &
                           density_core1(1:n_comps))
      rhoat=(d_dummy+mbin)/rhoat
      dw=((d_dummy+mbin)*6._wp/(pi*rhoat))**onethird

      vdry=mbin*sum(mass_frac_aer1(n_sel,1:n_comps)/density_core1(1:n_comps))
      vwat=max(d_dummy,0._wp)/rhow
      if (vdry.gt.tiny(1._wp)) then
          kappa=sum(mbin*mass_frac_aer1(n_sel,1:n_comps)/density_core1(1:n_comps)* &
               kappa_core1(1:n_comps))/vdry
      else
          kappa=0._wp
      endif
      if (kappa.le.tiny(1._wp)) then
          aw=1._wp
      else
          aw=vwat/(vwat+kappa*vdry)
      endif

      sigma=surface_tension(parcel1%t)
      mbin_comp=mbin*mass_frac_aer1(n_sel,1:n_comps)
      fads=fhh_adsorption_factor(d_dummy,dw,mbin_comp,density_core1(1:n_comps))

      kkoehler03=mult*(exp(4._wp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           aw*fads)-rh_act
    end function kkoehler03
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
	! ============================================================================
	! terminal01
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates liquid-drop terminal velocity, Reynolds number and drag coefficient over the model
	!>size range. see pruppacher and klett    
	!>@param[inout] vel: terminal fall speed
	!>@param[in] diam: wet drop diameter
	!>@param[in] rhoat: mean drop density
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[inout] nre: drop Reynolds number
	!>@param[inout] cd: drag coefficient
	!>@param[in] sz: number of drops/bins to process
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


	! ============================================================================
	! ventilation01
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates vapour and heat ventilation factors for liquid drops from their current terminal
	!>velocities and transport properties.
    !> see  pruppacher and klett - page 538-541                                     !
    !> original reference: pruppacher and rasmussen 1979                            !
	!>@param[in] diam: wet drop diameter
	!>@param[in] rhoat: mean drop density
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[inout] fv: ventilation factor for vapour diffusion
	!>@param[inout] fh: ventilation factor for heat transfer
	!>@param[in] sz: number of drops/bins to process
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
      calc = (nsc1**(onethird)) * sqrt(nre)
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
      calc = (nsc2**(onethird)) * sqrt(nre)
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


	! ============================================================================
	! maxdimension01
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates maximum dimension of unrimed ice aggregates and the equivalent spherical diameter of
	!>their rime mass.
	!>@param[in] mice: unrimed/depositional ice mass per representative particle
	!>@param[in] rhoi: mean depositional-ice density
	!>@param[in] phi: mean monomer aspect ratio
	!>@param[in] nump: mean number of monomers per aggregate
	!>@param[in] rime: rime mass per representative particle
	!>@param[out] dmax: maximum particle dimension including the provisional rime envelope
	!>@param[out] drime: equivalent solid-ice diameter of the rime mass
	!>@param[in] sz: number of ice particles/bins to process
	subroutine maxdimension01(mice,rhoi,phi,nump,rime,dmax,drime,sz)
		implicit none
		integer(i4b), intent(in) :: sz
		real(wp), dimension(:), intent(in) :: &
			mice,rhoi,phi,nump,rime
		real(wp), dimension(sz), intent(out) :: &
			dmax,drime
		integer(i4b) :: i
		real(wp) :: nmon
		real(wp) :: phi1,rho1
		real(wp) :: vmon
		real(wp) :: a,c,dmon,dagg
		real(wp), parameter :: small=1.e-60_wp

		dmax=0._wp
		drime=0._wp
		
		do i=1,sz
			! -------------------------------------------------------------
			! Equivalent solid-rime diameter.
			!
			! This is only used as a provisional rime-envelope treatment.
			! -------------------------------------------------------------
			if (rime(i) > small) then
				drime(i)= &
					(6._wp*rime(i)/(pi*rhoice))**onethird
			endif
			! -------------------------------------------------------------
			! Unrimed ice aggregate
			! -------------------------------------------------------------
			if (mice(i) > small) then
				nmon=max(nump(i),1._wp)
				phi1=max(phi(i),1.e-6_wp)
				rho1=max(rhoi(i),1._wp)
				! Do not allow mean deposited-ice density above solid ice
				rho1=min(rho1,rhoice)
				! Mean material volume of one monomer
				vmon=mice(i)/(rho1*nmon)
				! Chen-Lamb spheroidal monomer:
				!
				! phi = c/a
				!
				! V = 4/3*pi*a^2*c
				!   = 4/3*pi*a^3*phi
				!
				a=(3._wp*vmon/(4._wp*pi*phi1))**onethird
				c=a*phi1
				! Maximum dimension of one monomer
				dmon=2._wp*max(a,c)
				! Connolly et al. aggregate scaling:
				! fractal mass dimension ~= 2
				dagg=dmon*sqrt(nmon)
			else
				dagg=0._wp
			endif
			! -------------------------------------------------------------
			! PROVISIONAL RIME TREATMENT
			!
			! The Connolly et al. (2012) experiments effectively had
			! riming switched off, so their aggregate scaling does not
			! prescribe partially-rimed geometry.
			!
			! For now do not allow the particle to be smaller than the
			! equivalent sphere containing the rime itself.
			! -------------------------------------------------------------
			dmax(i)=max(dagg,drime(i))
		enddo
	end subroutine maxdimension01
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! areaaggregates01
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates projected area of ice monomers/aggregates using the current habit, aggregate scaling
	!>and provisional rime envelope.
	!>@param[out] area: projected particle area
	!>@param[in] mice: unrimed/depositional ice mass per representative particle
	!>@param[in] rhoi: mean depositional-ice density
	!>@param[in] phi: mean monomer aspect ratio
	!>@param[in] nump: mean number of monomers per aggregate
	!>@param[in] dmax: maximum particle dimension
	!>@param[in] drime: equivalent spherical diameter of rime mass
	!>@param[in] sz: number of ice particles/bins to process
	subroutine areaaggregates01( &
		area,mice,rhoi,phi,nump,dmax,drime,sz)
		implicit none
		integer(i4b), intent(in) :: sz
		real(wp), dimension(:), intent(in) :: mice,rhoi,phi,nump,dmax,drime
		real(wp), dimension(sz), intent(out) :: area
		integer(i4b) :: i
		real(wp), parameter :: dfract=1.33_wp
		real(wp), parameter :: small=1.e-60_wp
		real(wp) :: nmon
		real(wp) :: phi1,rho1
		real(wp) :: vmon
		real(wp) :: a,c,dmon,dagg
		real(wp) :: amon
		real(wp) :: aagg
		real(wp) :: arime
		real(wp) :: acirc
	
	
		area=0._wp
		do i=1,sz
			if (mice(i) > small) then
				nmon=max(nump(i),1._wp)
				phi1=max(phi(i),1.e-6_wp)
				rho1=min(max(rhoi(i),1._wp),rhoice)
				! Mean monomer material volume
				vmon=mice(i)/(rho1*nmon)

				! Monomer spheroid
				a=(3._wp*vmon/(4._wp*pi*phi1))**onethird	
				c=a*phi1
				dmon=2._wp*max(a,c)
	
				! Unrimed aggregate maximum dimension
				dagg=dmon*sqrt(nmon)
				! ---------------------------------------------------------
				! Projected monomer area normal to fall
				! ---------------------------------------------------------
				if (phi1 <= 1._wp) then
					! Oblate / plate
					amon=pi*a**2
				else
					! Prolate / column
					amon=pi*a*c
				endif
				! ---------------------------------------------------------
				! Aggregate projected area:
				!
				! A = c Dmax^1.33
				!
				! Set prefactor from monomer area.
				! ---------------------------------------------------------
				if (dmon > small) then
					aagg=amon*(dagg/dmon)**dfract
				else
					aagg=0._wp
				endif
			else
				aagg=0._wp
			endif
			! -------------------------------------------------------------
			! PROVISIONAL RIME ENVELOPE
			!
			! Do not let projected area be smaller than the area of the
			! equivalent solid-rime sphere.
			! -------------------------------------------------------------
			arime=pi*onequarter*drime(i)**2
			area(i)=max(aagg,arime)
			! -------------------------------------------------------------
			! Physical constraint:
			!
			! projected area cannot exceed a circle with diameter Dmax
			! -------------------------------------------------------------
			acirc=pi*onequarter*dmax(i)**2
			area(i)=min(area(i),acirc)
			area(i)=max(area(i),0._wp)
		enddo
	
	end subroutine areaaggregates01
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! terminal02
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates ice-particle terminal velocity and Reynolds number using the Heymsfield-Westbrook
	!>area-ratio formulation and current aggregate geometry.
    !> calculate the terminal velocity of ice particles                             
    !> see heymsfield and westbrook (2010, jas)                                     
    !> who corrected a bias in the fall speeds derived by mitchell and others       
    !> this method also works reasonably well for sub-100 micron crystals           
    !> see westbrook qj paper for latter point                                      
	!>@param[inout] vel: terminal fall speed
	!>@param[in] mwat: total ice-particle mass per representative particle
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[in] phi: mean monomer aspect ratio
	!>@param[in] rhoi: mean depositional-ice density
	!>@param[in] nump: mean number of monomers per aggregate
	!>@param[in] rime: rime mass per representative particle
	!>@param[inout] nre: Reynolds number
	!>@param[in] sz: number of ice particles/bins to process
	!>@param[out] dmax_out: optional maximum-dimension diagnostic
	!>@param[out] area_out: optional projected-area diagnostic
	subroutine terminal02( &
		vel,mwat,t,p,phi,rhoi,nump,rime,nre,sz,dmax_out,area_out)
      use numerics_type
      implicit none
      real(wp), intent(in) :: t, p
      real(wp), dimension(:), intent(in) :: mwat, phi, rhoi, nump, rime
      integer(i4b), intent(in) :: sz
      real(wp), dimension(sz), intent(inout) :: nre,vel 
      real(wp), dimension(:), intent(out), optional :: dmax_out
	  real(wp), dimension(:), intent(out), optional :: area_out
      real(wp) :: eta, rhoa
      real(wp), dimension(sz) :: dmax,drime,area,ar, x
      integer(i4b) :: i	
      ! air properties
      eta = viscosity_air(t)
      rhoa = p/ra/t
    
      ! calculate the maximum dimension of the particle
	  ! ----------------------------------------------------------------------
   	  ! Ice particle geometry
	  ! ----------------------------------------------------------------------	
	  call maxdimension01( &
			max(mwat-rime,0._wp), &
			rhoi,phi,nump,rime, &
			dmax,drime,sz)
	  call areaaggregates01( &
			area, &
			max(mwat-rime,0._wp), &
			rhoi,phi,nump, &
			dmax,drime,sz)
	  ! ----------------------------------------------------------------------
	  ! Area ratio
	  !
	  ! Ar = projected area /
	  !      area of circumscribed circle based on Dmax
	  ! ----------------------------------------------------------------------
	  ar=1._wp
	  where (dmax > tiny(1._wp))
	  	ar=area/(pi*onequarter*dmax**2)
	  end where
	  ! The Heymsfield-Westbrook formulation is specifically designed to
	  ! retain the effects of low area ratio ice particles. Only impose a
	  ! very small numerical floor.
	  ar=min(max(ar,1.e-6_wp),1._wp)
	  ! Optional diagnostic outputs
	  if (present(dmax_out)) then
	  	dmax_out(1:sz)=dmax
	  endif
	  if (present(area_out)) then
		area_out(1:sz)=area
	  endif
  
      ! heymsfield and westbrook
	  x=0._wp
	  nre=0._wp
	  vel=0._wp
	  where ((mwat > tiny(1._wp)).and. (dmax > tiny(1._wp)))	
		! Modified Best number
		x=rhoa*8._wp*mwat*grav / &
			((eta**2)*pi*sqrt(ar))
		! Heymsfield-Westbrook Re-X relation
		nre=(8._wp**2)*onequarter * &
				( sqrt(1._wp + 4._wp*sqrt(x) / &
				((8._wp**2)*sqrt(0.35_wp))) - 1._wp )**2
		vel=eta*nre/(rhoa*dmax)	
	  end where
	  ! Viscous/Stokes-like correction
	  where ((mwat > tiny(1._wp)).and. (dmax > tiny(1._wp)).and. (nre < 1._wp))
	    vel = grav*mwat / (6._wp*pi*eta*0.465_wp*dmax*sqrt(ar))	
	  end where
	  where (isnan(vel))
			vel=0._wp
	  end where	
	  where (isnan(nre))
			nre=0._wp
	  end where
    end subroutine terminal02
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
	! ============================================================================
	! ventilation02
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates vapour and heat ventilation factors for ice particles using current shape, density,
	!>aggregation state and terminal velocity.
    !> calculate ventilation coefficients for ice crystals                          
    !> see page 553 p+k                                                             
    !> original reference: ji 1991 and wang and ji 1992                             
	!>@param[in] mwat: total ice-particle mass per representative particle
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[in] phi: mean monomer aspect ratio
	!>@param[in] rhoi: mean depositional-ice density
	!>@param[in] nump: mean number of monomers per aggregate
	!>@param[in] rime: rime mass per representative particle
	!>@param[inout] fv: ventilation factor for vapour diffusion
	!>@param[inout] fh: ventilation factor for heat transfer
	!>@param[in] sz: number of ice particles/bins to process
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
      calc1 = nsc1**(onethird)
      calc2 = nsc2**(onethird)
  
      ! columns
      nre2=min(nre,20._wp)
      where(phi.gt.1.0_wp)  
        x = calc1*sqrt(nre2)	
        fv = 1.0_wp - 0.000668_wp*x*onequarter + 2.39402_wp*((x*onequarter)**2._wp) + &
             0.73409_wp*((x*onequarter)**3._wp)-0.73911_wp*((x*onequarter)**4._wp)
        x = calc2*sqrt(nre2);	
        fh = 1.0_wp - 0.000668_wp*x*onequarter + 2.39402_wp*((x*onequarter)**2._wp) + &
             0.73409_wp*((x*onequarter)**3._wp)-0.73911_wp*((x*onequarter)**4._wp)
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
    
    

	! ============================================================================
	! dropgrowthrate01
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates diffusional liquid-drop radial growth rate including kinetic, thermal and optional
	!>ventilation corrections.
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[in] rh: ambient liquid-water saturation ratio/relative humidity used by the growth equation
	!>@param[in] rh_eq: equilibrium relative humidity at each particle surface
	!>@param[in] rhoat: mean wet-particle density
	!>@param[in] diam: wet particle diameter
	!>@param[in] sz: number of particles/bins to process
	!>@return dropgrowthrate01: radial growth rate of each liquid particle
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



	! ============================================================================
	! icegrowthrate01
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates diffusional ice mass-growth rate using particle capacitance, deposition kinetics,
	!>thermal resistance and optional ventilation.
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[in] rh_ice: ambient relative humidity with respect to ice
	!>@param[in] rh_eq: equilibrium relative humidity with respect to ice for each particle
	!>@param[in] mwat: total ice-particle mass per representative particle
	!>@param[in] mbin: aerosol-component masses carried by the ice particles
	!>@param[in] rhobin: aerosol-component densities
	!>@param[in] phi: mean monomer aspect ratio
	!>@param[in] rhoi: mean depositional-ice density
	!>@param[in] nump: mean monomer number per aggregate
	!>@param[in] rime: rime mass per representative particle
	!>@param[in] sz: number of ice particles/bins to process
	!>@return icegrowthrate01: vapour-deposition/sublimation mass-growth rate for each ice particle
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


	! ============================================================================
	! capacitance01
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates electrostatic/diffusional capacitance of oblate and prolate spheroidal ice particles
	!>from their mass, density and aspect ratio.
    !> calculates the capacitance of oblate (longer a) and prolate (longer c) 
    !> spheroids. a is the length of the prism axis (half of) and c the basal 
    !> (half of). see page 1214 of chen and lamb, jas, 1994.
	!>@param[in] mwat: ice-particle mass
	!>@param[in] phi: aspect ratio c/a
	!>@param[in] rhoice: particle/depositional ice density
	!>@param[in] numi: mean monomer number; retained for aggregate-capacitance extensions
	!>@param[in] rimemass: rime mass; retained for aggregate-capacitance extensions
	!>@param[in] sz: number of particles/bins to process
	!>@return capacitance01: ice-particle capacitance length
    function capacitance01(mwat,phi,rhoice,numi,rimemass,sz)
      use numerics_type
      implicit none
      real(wp), dimension(:), intent(in) :: mwat, phi, rhoice, numi, rimemass
      real(wp), dimension(sz) :: capacitance01,vol,a,c,ecc,dmax,drime
  
      integer(i4b), intent(in) :: sz
  
      vol=mwat/rhoice
  
      a=( 3._wp*vol/(4._wp*pi*phi) )**(onethird)
      c=a*phi
  
      ! Handle nearly spherical particles before evaluating either
      ! spheroidal expression; both analytic forms are 0/0 at phi=1.
      where(abs(phi-1._wp).lt.1.e-4_wp)
        ecc=0._wp
        capacitance01=a
      elsewhere(phi.lt.1._wp)
        ecc=sqrt(max(1._wp-phi**2._wp,0._wp))
        capacitance01=a*ecc/asin(ecc)
      elsewhere
        ecc=sqrt(max(1._wp-phi**(-2._wp),0._wp))
        capacitance01=c*ecc/log((1._wp+ecc)*phi)
      end where
      ! westbrook et al. (2008, jas): capacitance of aggregates is 0.25 times the 
      ! maximum dimension of the ice particle
!       call maxdimension01(mwat-rimemass,rhoice,phi,numi,rimemass,dmax,drime,sz)
!       where(numi.ge.2._wp)
!          capacitance01=0.25_wp*dmax
!       end where
  
    end function capacitance01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	! ============================================================================
	! koopnucrate
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates homogeneous ice nucleation rate in supercooled liquid water using 
	!> the Koop et al. (2000)
	!> water-activity formulation with pressure correction.
	!>@param[in] aw: liquid-water activity
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[in] sz: number of elements in the activity array
	!>@return koopnucrate: homogeneous nucleation rate for each element
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
    
    
    
	! ============================================================================
	! fparcelwarm
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates the warm-parcel ODE tendencies for condensation/evaporation, pressure, temperature,
	!>relative humidity, height, velocity and parcel geometry.
	!>@param[inout] neq: number of ODE equations
	!>@param[inout] tt: ODE integration time
	!>@param[inout] y: parcel state vector
	!>@param[inout] ydot: parcel-state time derivatives
	!>@param[inout] rpar: real workspace required by the ODE interface
	!>@param[inout] ipar: integer workspace required by the ODE interface
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
        real(wp) :: wv=0._wp, wl=0._wp, wi=0._wp, rm, rme, cpm, &
                  drv=0._wp, dri=0._wp,dri2=0._wp, dwl_micro=0._wp, &
                  rh,t,p,err,sl, w, &
                  te, qve, pe, var, dummy, rhoe, rhop, b, mu, w_e,dlnrho, wv_old, &
                  rm_old, ratio=1.0_wp, flux_old, flux_new, svp1, svp2, winv, &
                  mu_old, mu_now, dz_ent

        integer(i4b) :: i, j,iloc, ipart, ipr, ite, irh, iz,iw, ira

        
        ipart=parcel1%n_bin_modew
        ipr=parcel1%ipr
        ite=parcel1%ite
        irh=parcel1%irh
        iz =parcel1%iz
        iw =parcel1%iw
        ira = parcel1%ira
        if(.not.adiabatic_prof) ydot(ira)=0._wp

        if ((updraft_type==3).and.(tt>t_thresh)) then
            y(iw)=winit2*cos(2._wp*pi*(tt-t_thresh)/tau2)
        endif

        rh=y(irh)
        t=y(ite)
        p=y(ipr)
    

        ! check there are no negative values
        where(y(1:ipart).le.0.e1_wp)
            y(1:ipart)=1.e-22_wp
        end where


        ! calculate mixing ratios from rh, etc
        svp1=svp_liq(t)
        sl=svp1*rh/(p-svp1) ! saturation ratio
        sl=(sl*p/(1._wp+sl))/svp1
        wv=eps1*rh*svp1 / (p-svp1) ! vapour mixing ratio
		! calculate the moist gas constant
        rm=ra+wv*rv

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! lateral entrainment reducing drop number conc.                       !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if((.not. adiabatic_prof) .and. (.not.l_inhom)) then
			svp2=svp_liq(parcel1%yold(parcel1%ite))
			wv_old=eps1*parcel1%yold(parcel1%irh)* &
				svp2 / (parcel1%yold(parcel1%ipr)- svp2) ! wv mixing ratio
			rm_old = ra+wv_old*rv
			if(bubble_flag) then
				! actually the mass of a bubble
				flux_old = fourthirds*pi*parcel1%yold(parcel1%ira)**3 * &
					parcel1%yold(parcel1%ipr)/(parcel1%yold(parcel1%ite)*rm_old)
				flux_new = fourthirds*pi*y(parcel1%ira)**3 * &
					y(parcel1%ipr)/(y(parcel1%ite)*rm)
			else
				! actually the mass flux for a jet (conserve number flux)
				flux_old = pi*parcel1%yold(parcel1%ira)**2 * &
					parcel1%yold(parcel1%ipr)/(parcel1%yold(parcel1%ite)*rm_old) * &
					parcel1%yold(parcel1%iw)
				flux_new = pi*y(parcel1%ira)**2 * &
					y(parcel1%ipr)/(y(parcel1%ite)*rm) * y(parcel1%iw)
			endif
			mu_old=ent_rate/max(parcel1%yold(ira),tiny(1._wp))
			mu_now=ent_rate/max(y(ira),tiny(1._wp))
			
			dz_ent=max(y(iz)-parcel1%yold(iz),0._wp)
			
			ratio=exp(-0.5_wp*(mu_old+mu_now)*dz_ent)
		endif
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		! calculate lw mixing ratio
        wl=sum(parcel1%npart*y(1:ipart))*ratio             ! liquid mixing ratio
        if(ice_flag==1) wi=sum(parcel1%npartice*parcel1%yice(1:ipart))*ratio ! ice mr

        ! calculate specific heats
        cpm=cp+wv*cpv+wl*cpw+wi*cpi

        ! now calculate derivatives
        ! adiabatic parcel model
        ydot(iz )=y(iw)                         ! vertical wind
        if(chamber_force_pressure .or. chamber_force_temperature .or. chamber_force_qtot) then
			iloc=find_pos(parcel1%time_chamber(1:n_levels_c),TT)
			iloc=min(n_levels_c-1,iloc)
			iloc=max(1,iloc)
            if (chamber_force_pressure) then
                ydot(ipr)=parcel1%dp_chamber(iloc)
            else
                ydot(ipr)=-p/rm/t*grav*ydot(iz)
            endif
        else
	        ydot(ipr)=-p/rm/t*grav*ydot(iz)      ! hydrostatic equation
		endif
		
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
        ! Microphysical change of liquid-water mixing ratio.  The ratio
        ! accounts for the homogeneous increase in parcel mass/flux during the
        ! current ODE step; the stored particle concentrations are diluted after
        ! the accepted ODE step.
        dwl_micro=sum(ydot(1:ipart)*parcel1%npart)*ratio
        drv=-dwl_micro
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in temperature of parcel                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (chamber_force_qtot) then
            ! The measured total-water trajectory is an explicit external
            ! source/sink.  It is independent of BL mixing and particle loss.
            drv=drv+parcel1%dqtot_chamber(iloc)
        endif

        if (chamber_force_temperature) then
            ydot(ite)=parcel1%dt_chamber(iloc)
        else
            ydot(ite)=rm/p*ydot(ipr)*t/cpm  ! temperature change: expansion
            ydot(ite)=ydot(ite)-lv/cpm*drv ! temp change: condensation/source
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        if(.not.adiabatic_prof) then ! homogeneous lateral entrainment
			! parcel density:
			rhop=p/(rm*t)
			w_e=y(iw)
			! mu
			if(l_inhom) then
				mu=0.0_wp
			else
				mu=ent_rate/y(ira)
			endif		
			!calculate the environmental p, qv, te, density
			! parcel p, density
			! buoyancy...
			! locate position
			iloc=find_pos(parcel1%z_sound(1:n_levels_s),y(iz))
			iloc=min(n_levels_s-1,iloc)
			iloc=max(1,iloc)
			! linear interp p
			call poly_int(parcel1%z_sound(iloc:iloc+1), &
					parcel1%p_sound(iloc:iloc+1), &
						min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
			pe=p
			! linear interp qv
			call poly_int(parcel1%p_sound(iloc:iloc+1), &
					parcel1%q_sound(1,iloc:iloc+1), &
						max(y(ipr),parcel1%p_sound(n_levels_s)), var,dummy)        
			qve=var
			! linear interp te
			call poly_int(parcel1%p_sound(iloc:iloc+1), &
					parcel1%t_sound(iloc:iloc+1), &
						max(y(ipr),parcel1%p_sound(n_levels_s)), var,dummy)        
			te=var
			! env density:
			rme=ra+qve*rv
			rhoe=pe/(rme*te)
			!buoyancy
			if((parcel1%z_sound(n_levels_s) .lt. y(iz)) .or. &
				(parcel1%z_sound(1) .gt. y(iz))) then
				b=0._wp
			else
!                 b=grav*(rhoe-rhop)/rhop
				b=grav*(t-te)/te
			endif
			if(updraft_type==4) then
				if(y(iw) > 0.001_wp) then
					ydot(iw)=gam_fac_ent*(b-grav*wl-grav*wi)-mu*gam_fac_ent*w_e**2 
				else
					ydot(iw)=(0.001_wp-y(iw))/10._wp
				endif
			endif        
            ! Pruppacher and Klett (1997), Eqs. 12-29 and 12-26.
            !
            ! Eq. 12-29 is written in terms of the TOTAL derivative of liquid
            ! water mixing ratio:
            !
            !   d wv/dt = -d wl/dt - mu*W*(wv-wv_env+wl).
            !
            ! In this bin implementation the homogeneous dilution part of wl
            ! is already represented by ratio=old_mass/new_mass, i.e.
            !
            !   d wl/dt = d wl_micro/dt - mu*W*wl.
            !
            ! Substitution cancels the explicit +wl term.  Therefore only the
            ! microphysical liquid tendency and vapour-vapour mixing belong in
            ! drv here.  This avoids applying condensate dilution twice.
            drv=drv-w_e*mu*(wv-qve)

            ! With the same substitution Eq. 12-26 reduces to the existing
            ! expansion/latent-heating terms above plus direct sensible mixing
            ! of parcel and environmental temperature.
            ydot(ite)=ydot(ite)-w_e*mu*(t-te)

            ! Equation 12-32 or 12-34 P+K
            !dlnrho=1._wp/rhop*(1._wp/(ra*t))*(ydot(ipr)- p/(t)*ydot(ite))
            ! includes contribution from rv
            dlnrho = ydot(ipr)/p - ydot(ite)/t - rv*drv/rm
            if(bubble_flag) then
                ydot(ira) = y(ira)*onethird*(mu*w_e-dlnrho)
            else
            	winv = w_e/(w_e**2 + 1.e-3_wp**2)
                ydot(ira) = y(ira)*0.5_wp*(mu*w_e- dlnrho - &
                    winv*ydot(iw))
            endif

        endif
        ! jobs: 
        !       add 2nd aerosol structure, which can be entrained in
        !       entrainment of aerosol / drops outside of solver 
        !       move entrainment outside of solver?     

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in rh of parcel                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ydot(irh)=(p-svp1)*svp1*drv
        ydot(irh)=ydot(irh)+svp1*wv*ydot(ipr)
        ydot(irh)=ydot(irh)-wv*p*dfsid1(svp_liq,t,1.e0_wp,1.e-8_wp,err)*ydot(ite)
        ydot(irh)=ydot(irh) / (eps1*svp1**2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
    end subroutine fparcelwarm
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! fparcelcold
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Evaluates the ice-parcel ODE tendencies for vapour deposition/sublimation and the associated
	!>thermodynamic changes.
	!>@param[inout] neq: number of ODE equations
	!>@param[inout] tt: ODE integration time
	!>@param[inout] y: ice/parcel state vector
	!>@param[inout] ydot: ice/parcel-state time derivatives
	!>@param[inout] rpar: real workspace required by the ODE interface
	!>@param[inout] ipar: integer workspace required by the ODE interface
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
        								! the pressure change is zero, so 0
        if (.not.chamber_force_temperature) &
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



	! ============================================================================
	! jparcelwarm
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Dummy Jacobian callback retained for compatibility with the parcel ODE solver interface; no
	!>Jacobian entries are currently calculated.
	!>@param[in] neq: number of ODE equations
	!>@param[in] t: ODE integration time
	!>@param[in] y: current solution vector
	!>@param[in] ml,mu: lower and upper Jacobian bandwidth parameters
	!>@param[out] pd: Jacobian work array; currently left unchanged
	!>@param[in] nrpd: leading dimension of pd
	!>@param[inout] rpar: real solver workspace
	!>@param[inout] ipar: integer solver workspace
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

    

    ! ============================================================================
    ! daily_dcmex_2025
    ! ============================================================================
    !>@author
    !>Paul J. Connolly, The University of Manchester
    !>@brief
    !>Temperature-only DCMEX INP concentration spectrum from the supplied
    !>Daily/Daly et al. example.  The returned concentration is in m-3.
    !>
    !>The spectrum is the sum of (1) a fixed mineral-dust/K-feldspar term using
    !>the Harrison et al. polynomial and (2) a warm biological/empirical term.
    !>It is zero at and above -4 deg C.  The supplied example evaluates the
    !>curve down to -40 deg C, so colder temperatures are clamped at -40 deg C.
    !>
    !>This function deliberately has no aerosol-size argument.  The caller caps
    !>the requested INP number by the number of eligible aerosol particles.
    !>@param[in] t: temperature [K]
    !>@return daily_dcmex_2025: target INP number concentration [m-3]
    function daily_dcmex_2025(t) result(ninp)
        implicit none
        real(wp), intent(in) :: t
        real(wp) :: ninp
        real(wp) :: tc,teval,log_ns,ns_site,ninp1,ninp2
        real(wp), parameter :: adust=7.5e-6_wp ! 7.5 um2 cm-3 -> m2 m-3
        real(wp), parameter :: fkfsp=0.05_wp

        tc=t-ttr
        ninp=0._wp
        if(tc.ge.-4._wp) return

        teval=max(tc,-40._wp)

        ! Harrison et al. K-feldspar polynomial as used in the supplied
        ! DCMEX example.  log_ns is log10(sites cm-2), hence the 1e4 conversion
        ! to sites m-2.
        log_ns = -3.25_wp - 0.793_wp*teval - &
                 6.91e-2_wp*teval**2 - 4.17e-3_wp*teval**3 - &
                 1.05e-4_wp*teval**4 - 9.08e-7_wp*teval**5
        ns_site=10._wp**log_ns*1.e4_wp
        ninp1=adust*fkfsp*ns_site

        ! a*exp(b*(Tmax-T)^c), written exactly as in the supplied example:
        ! 1000 converts L-1 to m-3 and exp(-50)=1.9287e-22.
        ninp2=1000._wp*exp(-50._wp + &
            45.25_wp*(-4._wp-teval)**0.046_wp)

        ninp=max(ninp1+ninp2,0._wp)
    end function daily_dcmex_2025

	! ============================================================================
	! demott_2010
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates primary ice-nucleating-particle number concentration using the DeMott et al. (2010)
	!>parameterisation.
	!>  (Predicting global atmospheric ice...
	!>                    https://doi.org/10.1073/pnas.0910818107)
	!>@param[in] t: temperature
	!>@param[in] naer05: aerosol number concentration with dry diameter greater than 0.5 micrometres
	!>@return demott_2010: predicted INP number concentration
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
		
! 		demott_2010=0._wp
! 		if(t<270.15_wp) then
! 			demott_2010=(10.e3-1.e3)/(3.0-25.0)*(t-273.15)-227.2727272727273
! 		endif
		
	end function demott_2010
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! ============================================================================
	! noncollisional_iceformation
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Forms ice without particle-particle collisions, transfers liquid particles and conserved
	!>component moments to the fixed ice grid, and applies latent heating.
	!>@param[inout] npart: liquid-particle number concentration
	!>@param[inout] npartice: ice-particle number concentration
	!>@param[in] mwat: liquid-water mass per representative drop
	!>@param[in] mbin2: aerosol-component masses in liquid particles
	!>@param[inout] mbin2_ice: aerosol-component and ice masses in ice particles
	!>@param[in] rhobin,nubin,kappabin,molwbin: aerosol-component thermodynamic properties
	!>@param[inout] moments: conserved liquid and ice moments
	!>@param[in] medges: fixed mass-bin edges for each mode
	!>@param[inout] t: parcel temperature, updated for latent heat of freezing
	!>@param[in] p: pressure
	!>@param[in] nbins1: number of fixed bins per mode
	!>@param[in] ncomps: number of aerosol composition components
	!>@param[in] nbinw: number of liquid/warm bins
	!>@param[in] nmoms: number of conserved moments
	!>@param[in] nmodes: number of external modes
	!>@param[inout] yice: ice mass per representative particle
	!>@param[inout] rh: ambient relative humidity
	!>@param[in] dt: timestep
	!>@param[in] sce_flag: SCE switch controlling fixed-grid treatment
	!>@param[in] mode1_flag: switch for mode-1 secondary-ice treatment
	!>@param[in] ice_nucleation_mech_in: logical switches for Koop, INAS, DeMott and Daily/DCMEX
    subroutine noncollisional_iceformation(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin, &
                         moments,medges, &
                         t,p,nbins1,ncomps,nbinw,nmoms,nmodes,yice,rh,dt,&
                         sce_flag,mode1_flag, ice_nucleation_mech_in)
    use numerics_type
    use sce, only : calculate_mode1, sce_receiving_bin
    implicit none
    real(wp), intent(inout) :: t
    real(wp), intent(in) :: p,dt
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
    logical, dimension(N_INUC_MECH), intent(in) :: ice_nucleation_mech_in

    real(wp), dimension(nbinw) :: nw,aw,jw,dn01,dn_inas,dn_daily,dn_demott,dn_koop, &
                                  nremain,m01,dw,dd,kappa,rhoat
    real(wp), dimension(nbinw,n_inp_classes) :: nin_freeze
    real(wp), dimension(nmoms) :: momtemp

    real(wp), intent(inout), dimension(nbinw) :: yice
    real(wp), intent(inout) :: rh

    integer(i4b) :: i,j,k,kk,im,inew,it,ib,jl,jh,imode,ibin,iedge
    real(wp) :: fracinliq,naer05,naer_daily,nprimary,nprimary_existing, &
                ndaily_target,ndaily_existing,avail, &
                n,nt,nb,mt,mb,mnew,nleft,mttot,mbtot,mleft,mall, &
                dprimary,dn,frac,dcin,tc,dlo,dhi,frac05
    logical :: has_inas,has_demott,has_daily,full_moving
    logical, dimension(nbinw) :: activated_mask
    integer(i4b) :: idemott

    ! Non-collisional freezing is only relevant below the melting point.
    if(t.gt.ttr) return

    tc=t-ttr
    m01=yice*npartice
    dn_inas=0._wp
    dn_daily=0._wp
    dn_demott=0._wp
    dn_koop=0._wp
    dn01=0._wp
    nin_freeze=0._wp
    idemott=ncomps+6+n_inp_classes
    full_moving=(parcel1%bin_scheme_flag.eq.BIN_FULL_MOVING)

    ! Evaluate the expensive Koehler/FHH activation test once per warm bin and
    ! reuse it for all heterogeneous immersion-nucleation mechanisms.
    activated_mask=.false.
    if (ice_nucleation_mech_in(INUC_INAS) .or. &
        ice_nucleation_mech_in(INUC_DEMOTT) .or. &
        ice_nucleation_mech_in(INUC_DAILY)) then
        do i=1,nbinw
            if (npart(i).le.qsmall2 .or. mwat(i).le.tiny(1._wp)) cycle
            activated_mask(i)=particle_is_activated(i,mwat(i),t)
        enddo
    endif

    ! -------------------------------------------------------------------------
    ! (1) Immersion freezing by the prognostic discrete INAS spectrum.
    ! Only activated liquid particles (current water mass above the
    ! Koehler/FHH critical water mass) are eligible.
    !
    ! The IN moments are cumulative and ordered from warm to cold threshold.
    ! When a threshold is crossed, the remaining cumulative value at that
    ! threshold is the number in that cohort after warmer cohorts have already
    ! been removed.  The intrinsic IN spectrum is transferred with the frozen
    ! particles; ns(T) is not evaluated here.
    ! -------------------------------------------------------------------------
    if(ice_nucleation_mech_in(INUC_INAS) .and. n_inp_classes.gt.0) then

        do i=1,nbinw
            if(npart(i).le.qsmall2) cycle
            if(mwat(i).le.tiny(1._wp)) cycle
            if(.not.activated_mask(i)) cycle

            call get_inp_control(mbin2(i,:),has_inas,has_demott,has_daily)
            if(.not.has_inas) cycle

            do kk=1,n_inp_classes
                if(tc.gt.inp_temp(kk)) cycle

                im=ncomps+5+kk
                dn=max(moments(i,im),0._wp)
                dn=min(dn,max(npart(i)-dn_inas(i),0._wp))
                if(dn.le.qsmall2) cycle

                ! Remove one differential threshold cohort from all cumulative
                ! moments that contain it.
                moments(i,im:ncomps+5+n_inp_classes)= &
                    max(moments(i,im:ncomps+5+n_inp_classes)-dn,0._wp)

                ! The same intrinsic threshold identity follows the particles
                ! into the ice phase.
                nin_freeze(i,kk:n_inp_classes)= &
                    nin_freeze(i,kk:n_inp_classes)+dn
                dn_inas(i)=dn_inas(i)+dn
            enddo
        enddo
    endif

    ! -------------------------------------------------------------------------
    ! (2) DeMott et al. (2010) for particles assigned to the DeMott category.
    ! Explicit INAS components take precedence after internal mixing.
    ! Existing DeMott primary ice is tracked by the extensive n_demott moment,
    ! so aggregation with other ice types does not erase DeMott provenance.
    ! -------------------------------------------------------------------------
    if(ice_nucleation_mech_in(INUC_DEMOTT)) then
        dd(:)=((sum(mbin2(:,:)/rhobin(:,:),2))*6._wp/pi)**onethird
        naer05=0._wp

        do i=1,nbinw
            if(mwat(i).le.tiny(1._wp)) cycle
            if(.not.activated_mask(i)) cycle
            call get_inp_control(mbin2(i,:),has_inas,has_demott,has_daily)
            if((has_inas .and. ice_nucleation_mech_in(INUC_INAS)) .or. .not.has_demott) cycle

			if (full_moving) then
				! Full-moving bins retain their original dry-aerosol-bin identity.
				! Integrate a flat number distribution across the original dry bin
				! when it straddles the DeMott 0.5 micron threshold.
				imode=(i-1)/nbins1+1
				ibin=i-(imode-1)*nbins1
				iedge=ibin+(imode-1)*(nbins1+1)
			
				dlo=parcel1%d(iedge)
				dhi=parcel1%d(iedge+1)
			
				if (dhi.le.0.5e-6_wp) then
					frac05=0._wp
				elseif (dlo.ge.0.5e-6_wp) then
					frac05=1._wp
				elseif (dhi.gt.dlo) then
					frac05=(dhi-0.5e-6_wp)/(dhi-dlo)
				else
					frac05=0._wp
				endif
			else
				! Moving-centre and Chen-Lamb remap particles between numerical bins,
				! so the original parcel1%d interval is no longer associated with bin i.
				! Classify using the representative dry-equivalent diameter reconstructed
				! from the aerosol component masses carried by the remapped population.
				frac05=merge(1._wp,0._wp,dd(i).gt.0.5e-6_wp)
			endif
           naer05=naer05+frac05*max(npart(i)-dn_inas(i),0._wp)
        enddo

        nprimary=min(naer05,demott_2010(t,naer05))

        ! Existing DeMott primary monomers are tracked explicitly.
        nprimary_existing=sum(max(moments(nbinw+1:2*nbinw,idemott),0._wp))
        nprimary=max(nprimary-nprimary_existing,0._wp)

        ! Deplete eligible aerosol from the large end, retaining the existing
        ! BMM ordering used by the DeMott implementation.
        do i=nbinw,1,-1
            if(nprimary.le.qsmall2) exit
            if(mwat(i).le.tiny(1._wp)) cycle
            if(.not.activated_mask(i)) cycle

            call get_inp_control(mbin2(i,:),has_inas,has_demott,has_daily)
            if((has_inas .and. ice_nucleation_mech_in(INUC_INAS)) .or. .not.has_demott) cycle

			if (full_moving) then
				! Full-moving bins retain their original dry-aerosol-bin identity.
				! Integrate a flat number distribution across the original dry bin
				! when it straddles the DeMott 0.5 micron threshold.
				imode=(i-1)/nbins1+1
				ibin=i-(imode-1)*nbins1
				iedge=ibin+(imode-1)*(nbins1+1)
			
				dlo=parcel1%d(iedge)
				dhi=parcel1%d(iedge+1)
			
				if (dhi.le.0.5e-6_wp) then
					frac05=0._wp
				elseif (dlo.ge.0.5e-6_wp) then
					frac05=1._wp
				elseif (dhi.gt.dlo) then
					frac05=(dhi-0.5e-6_wp)/(dhi-dlo)
				else
					frac05=0._wp
				endif
			else
				! Moving-centre and Chen-Lamb remap particles between numerical bins,
				! so the original parcel1%d interval is no longer associated with bin i.
				! Classify using the representative dry-equivalent diameter reconstructed
				! from the aerosol component masses carried by the remapped population.
				frac05=merge(1._wp,0._wp,dd(i).gt.0.5e-6_wp)
			endif

            avail=max(frac05*max(npart(i)-dn_inas(i),0._wp)-dn_demott(i),0._wp)
            dprimary=min(nprimary,avail)
            dn_demott(i)=dn_demott(i)+dprimary
            nprimary=nprimary-dprimary
        enddo
    endif

    ! -------------------------------------------------------------------------
    ! (3) Daily/Daly et al. DCMEX concentration spectrum.
    !
    ! This is a temperature-only target number concentration: there is no dry
    ! diameter criterion.  It is therefore applied proportionally across every
    ! eligible liquid aerosol bin.  The new ice requested at this timestep is
    ! min(max(N_target-N_existing,0),N_available), which explicitly limits the
    ! parameterisation by the available aerosol reservoir.
    !
    ! Runtime precedence among enabled mechanisms is INAS > DeMott > Daily.
    ! Hence Daily is used only when no enabled higher-priority treatment applies.
    ! -------------------------------------------------------------------------
    if(ice_nucleation_mech_in(INUC_DAILY)) then
        naer_daily=0._wp

        do i=1,nbinw
            if(npart(i).le.qsmall2) cycle
            if(mwat(i).le.tiny(1._wp)) cycle
            if(.not.activated_mask(i)) cycle
            call get_inp_control(mbin2(i,:),has_inas,has_demott,has_daily)
            if((has_inas .and. ice_nucleation_mech_in(INUC_INAS)) .or. &
               (has_demott .and. ice_nucleation_mech_in(INUC_DEMOTT)) .or. &
               .not.has_daily) cycle
            naer_daily=naer_daily+max(npart(i)-dn_inas(i)-dn_demott(i),0._wp)
        enddo

        ! Daily has no dedicated provenance moment.  Use carried aerosol
        ! composition plus ice monomer number to diagnose already nucleated
        ! Daily/DCMEX primary ice.  After mixed-mechanism aggregation this is
        ! an intentional approximation; n_demott does not have this ambiguity.
        ndaily_existing=0._wp
        do i=1,nbinw
            if(moments(nbinw+i,ncomps+2).le.qsmall2) cycle
            call get_inp_control(mbin2_ice(i,1:ncomps), &
                                 has_inas,has_demott,has_daily)
            if((has_inas .and. ice_nucleation_mech_in(INUC_INAS)) .or. &
               (has_demott .and. ice_nucleation_mech_in(INUC_DEMOTT)) .or. &
               .not.has_daily) cycle
            ndaily_existing=ndaily_existing+moments(nbinw+i,ncomps+2)
        enddo

        ndaily_target=max(daily_dcmex_2025(t)-ndaily_existing,0._wp)
        ndaily_target=min(ndaily_target,naer_daily)

        ! No aerosol-size dependence: sample all eligible liquid bins in
        ! proportion to their currently available aerosol number.
        if(ndaily_target.gt.qsmall2 .and. naer_daily.gt.qsmall2) then
            do i=1,nbinw
                if(mwat(i).le.tiny(1._wp)) cycle
                if(.not.activated_mask(i)) cycle
                call get_inp_control(mbin2(i,:), &
                                     has_inas,has_demott,has_daily)
                if((has_inas .and. ice_nucleation_mech_in(INUC_INAS)) .or. &
                   (has_demott .and. ice_nucleation_mech_in(INUC_DEMOTT)) .or. &
                   .not.has_daily) cycle

                avail=max(npart(i)-dn_inas(i)-dn_demott(i),0._wp)
                dn_daily(i)=min(avail,ndaily_target*avail/naer_daily)
            enddo
        endif
    endif

    ! -------------------------------------------------------------------------
    ! (4) Koop homogeneous freezing of whatever liquid remains.
    ! Homogeneous freezing samples the remaining intrinsic IN spectrum
    ! proportionally rather than preferentially selecting a threshold class.
    ! -------------------------------------------------------------------------
    nremain=max(npart-dn_inas-dn_daily-dn_demott,0._wp)

    if(ice_nucleation_mech_in(INUC_KOOP)) then
        nw=mwat/molw_water

        select case(kappa_flag)
        case(0)
            aw(:)=nw(:)/(nw(:)+sum(mbin2(:,:)/molwbin(:,:)*nubin(:,:),2))

        case(1)
            rhoat(:)=mwat(:)/rhow+sum(mbin2(:,:)/rhobin(:,:),2)
            rhoat(:)=(mwat(:)+sum(mbin2(:,:),2))/max(rhoat(:),tiny(1._wp))
            dw(:)=((mwat(:)+sum(mbin2(:,:),2))*6._wp/ &
                    (pi*max(rhoat(:),tiny(1._wp))))**onethird
            dd(:)=((sum(mbin2(:,:)/rhobin(:,:),2))*6._wp/pi)**onethird
            kappa(:)=sum((mbin2(:,:)+1.e-60_wp)/rhobin(:,:)*kappabin(:,:),2) / &
                     sum((mbin2(:,:)+1.e-60_wp)/rhobin(:,:),2)
            aw=(dw**3-dd**3)/(dw**3-dd**3*(1._wp-kappa))

        case default
            error stop 'Unknown kappa_flag in noncollisional_iceformation'
        end select

        jw=koopnucrate(aw,t,p,nbinw)
        dn_koop=abs(nremain*(1._wp-exp(-jw*mwat/rhow*dt)))
        dn_koop=min(dn_koop,nremain)

        if(n_inp_classes.gt.0) then
            do i=1,nbinw
                if(nremain(i).le.qsmall2 .or. dn_koop(i).le.qsmall2) cycle
                frac=min(max(dn_koop(i)/nremain(i),0._wp),1._wp)
                do kk=1,n_inp_classes
                    im=ncomps+5+kk
                    dcin=frac*moments(i,im)
                    moments(i,im)=max(moments(i,im)-dcin,0._wp)
                    nin_freeze(i,kk)=nin_freeze(i,kk)+dcin
                enddo
            enddo
        endif
    endif

    dn01=min(dn_inas+dn_daily+dn_demott+dn_koop,npart)

    ! -------------------------------------------------------------------------
    ! (5) Freezing fragmentation and transfer to the ice grid.
    ! -------------------------------------------------------------------------
    do i=1,nmodes
        do j=1,nbins1
            k=j+(i-1)*nbins1
            if(dn01(k).le.qsmall2) cycle

            mall=mwat(k)*dn01(k)
            mnew=mwat(k)
            nleft=dn01(k)
            n=0._wp
            nt=0._wp
            nb=0._wp
            mt=0._wp
            mb=0._wp

            if(mode1_flag) call calculate_mode1(mwat(k),0._wp,t,n,nt,nb,mb,mt)

            mnew=mnew-nt*mt-nb*mb
            n=n*dn01(k)
            nt=nt*dn01(k)
            nb=nb*dn01(k)
            mttot=mt*nt
            mbtot=mb*nb
            mleft=nleft*mnew

            ! Number concentration remaining in the liquid bin.
            npart(k)=npart(k)-dn01(k)
            fracinliq=npart(k)/max(npart(k)+dn01(k),1.e-30_wp)

            ! Aerosol component moments are sampled in proportion to particle
            ! number for freezing.  IN moments were already removed selectively
            ! above and are therefore not included in this scaling.
            momtemp=0._wp
            momtemp(1:ncomps)=(1._wp-fracinliq)*moments(k,1:ncomps)
            moments(k,1:ncomps)=moments(k,1:ncomps)*fracinliq
            moments(k,ncomps+1)=moments(k,ncomps+1)*fracinliq
            moments(k,ncomps+2)=moments(k,ncomps+2)*fracinliq
            moments(k,ncomps+4)=npart(k)*mwat(k)
            moments(k,ncomps+5)=npart(k)*mwat(k)

            if (full_moving) then
                ! Primary freezing must remain on the current moving ice
                ! representation.  Search only within the parent aerosol mode;
                ! do not project the new crystal back onto the fixed SCE grid.
                jl=(i-1)*nbins1+1
                jh=i*nbins1

                inew=sce_receiving_bin(mnew,jl,jh,yice,.true.)

                ! Mode-1 fragments are optional.  If a fragment class is absent,
                ! leave its destination equal to the primary destination; its
                ! zero number/mass contribution below then has no effect.
                it=inew
                ib=inew
                if (nt.gt.qsmall2 .and. mt.gt.qsmall2) &
                    it=sce_receiving_bin(mt,jl,jh,yice,.true.)
                if (nb.gt.qsmall2 .and. mb.gt.qsmall2) &
                    ib=sce_receiving_bin(mb,jl,jh,yice,.true.)
            else
                ! Existing moving-centre / Chen-Lamb fixed-grid placement.
                inew=find_medge(medges,mnew,nbins1,nmodes,i)
                it=max(find_medge(medges,mt,nbins1,nmodes,i),1)
                ib=max(find_medge(medges,mb,nbins1,nmodes,i),1)
            endif

            if(mall.gt.0._wp) then
                moments(inew+nbinw,1:ncomps)=moments(inew+nbinw,1:ncomps)+ &
                    momtemp(1:ncomps)*mleft/mall
                moments(it+nbinw,1:ncomps)=moments(it+nbinw,1:ncomps)+ &
                    momtemp(1:ncomps)*mttot/mall
                moments(ib+nbinw,1:ncomps)=moments(ib+nbinw,1:ncomps)+ &
                    momtemp(1:ncomps)*mbtot/mall

                ! The aerosol residual/IN identity is not cloned into every
                ! freezing fragment.  Its ensemble concentration is partitioned
                ! over fragment categories by the same mass fractions.
                if(n_inp_classes.gt.0) then
                    do kk=1,n_inp_classes
                        im=ncomps+5+kk
                        moments(inew+nbinw,im)=moments(inew+nbinw,im)+ &
                            nin_freeze(k,kk)*mleft/mall
                        moments(it+nbinw,im)=moments(it+nbinw,im)+ &
                            nin_freeze(k,kk)*mttot/mall
                        moments(ib+nbinw,im)=moments(ib+nbinw,im)+ &
                            nin_freeze(k,kk)*mbtot/mall
                    enddo
                endif

                ! DeMott provenance follows the original primary residual but
                ! is not cloned into every freezing fragment.  Partition its
                ! ensemble concentration with the same fragment mass fractions.
                moments(inew+nbinw,idemott)=moments(inew+nbinw,idemott)+ &
                    dn_demott(k)*mleft/mall
                moments(it+nbinw,idemott)=moments(it+nbinw,idemott)+ &
                    dn_demott(k)*mttot/mall
                moments(ib+nbinw,idemott)=moments(ib+nbinw,idemott)+ &
                    dn_demott(k)*mbtot/mall
            endif

            m01(inew)=m01(inew)+mleft
            m01(it)=m01(it)+mttot
            m01(ib)=m01(ib)+mbtot

            npartice(inew)=npartice(inew)+nleft
            npartice(it)=npartice(it)+nt
            npartice(ib)=npartice(ib)+nb

            ! Keep the moving pivots current during the freezing loop so the
            ! next source bin is placed relative to the already-updated ice PSD.
            ! For an empty receiving bin this sets the pivot exactly to the new
            ! product mass; for an occupied bin it gives the exact number-weighted
            ! mean mass.
            if (full_moving) then
                if (npartice(inew).gt.qsmall2) yice(inew)=m01(inew)/npartice(inew)
                if (npartice(it).gt.qsmall2)    yice(it)=m01(it)/npartice(it)
                if (npartice(ib).gt.qsmall2)    yice(ib)=m01(ib)/npartice(ib)
            endif

            ! Ice moments: phi, nmon and volume.
            moments(inew+nbinw,ncomps+1)=moments(inew+nbinw,ncomps+1)+nleft
            moments(it+nbinw,ncomps+1)=moments(it+nbinw,ncomps+1)+nt
            moments(ib+nbinw,ncomps+1)=moments(ib+nbinw,ncomps+1)+nb

            moments(inew+nbinw,ncomps+2)=moments(inew+nbinw,ncomps+2)+nleft
            moments(it+nbinw,ncomps+2)=moments(it+nbinw,ncomps+2)+nt
            moments(ib+nbinw,ncomps+2)=moments(ib+nbinw,ncomps+2)+nb

            moments(inew+nbinw,ncomps+3)=moments(inew+nbinw,ncomps+3)+mleft/rhoice
            moments(it+nbinw,ncomps+3)=moments(it+nbinw,ncomps+3)+mttot/rhoice
            moments(ib+nbinw,ncomps+3)=moments(ib+nbinw,ncomps+3)+mbtot/rhoice
        enddo
    enddo

    where((m01.gt.0._wp).and.(npartice.gt.0._wp))
        yice=m01/npartice
    end where

    mbin2_ice(:,1:ncomps)=moments(nbinw+1:2*nbinw,1:ncomps) / &
        (1.e-50_wp+spread(npartice,2,ncomps))
    mbin2_ice(:,ncomps+1)=yice

    ! Freezing releases latent heat.  Vapour mass is unchanged, so recompute RH
    ! at the updated temperature.
    if(.not.chamber_force_temperature) then
        call adjust_t_rh(sum(mwat(:)*dn01(:)),t,rh,p)
    endif

    end subroutine noncollisional_iceformation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
	! ============================================================================
	! find_medge
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Finds the flattened fixed-grid bin index immediately below a specified particle mass in a
	!>selected external mode.
	!>@param[in] medges: mass-bin edges for all modes
	!>@param[in] m: particle mass to locate
	!>@param[in] nbins1: number of bins per mode
	!>@param[in] nmodes: number of external modes
	!>@param[in] imode: external mode to search
	!>@return find_medge: flattened bin index below the requested mass
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
    
    
	! ============================================================================
	! icenucleation
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Legacy no-fragmentation ice-nucleation callback. Uses the same Koop, INAS, DeMott and
	!>Daily/DCMEX treatment, bin placement, moment transfer and latent heating as
	!>noncollisional_iceformation, but with Mode-1 freezing fragmentation disabled.
	!>@param[inout] npart: liquid-particle number concentration
	!>@param[inout] npartice: ice-particle number concentration
	!>@param[in] mwat: liquid-water mass per representative particle
	!>@param[in] mbin2: aerosol-component masses in liquid particles
	!>@param[inout] mbin2_ice: aerosol-component and ice masses in ice particles
	!>@param[in] rhobin,nubin,kappabin,molwbin: aerosol-component thermodynamic properties
	!>@param[inout] moments: conserved liquid and ice moments
	!>@param[inout] t: parcel temperature, updated for latent heat of freezing
	!>@param[in] p: pressure
	!>@param[in] sz: number of aerosol composition components
	!>@param[in] sz2: number of liquid/warm bins
	!>@param[in] sz3: number of conserved moments
	!>@param[inout] yice: ice mass per representative particle
	!>@param[inout] rh: ambient relative humidity
	!>@param[in] dt: timestep
    subroutine icenucleation(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin, &
                         moments,t,p,sz,sz2,sz3,yice,rh,dt)
        use numerics_type
        implicit none

        real(wp), intent(inout) :: t
        real(wp), intent(in) :: p,dt
        real(wp), dimension(sz2), intent(inout) :: npart,npartice
        real(wp), dimension(sz2), intent(in) :: mwat
        real(wp), dimension(sz2,sz), intent(in) :: mbin2, &
                                                  rhobin,nubin,kappabin,molwbin
        real(wp), dimension(2*sz2,sz3), intent(inout) :: moments
        integer(i4b), intent(in) :: sz,sz2,sz3
        real(wp), dimension(sz2,sz+1), intent(inout) :: mbin2_ice
        real(wp), dimension(sz2), intent(inout) :: yice
        real(wp), intent(inout) :: rh

        ! Keep this legacy callback numerically identical to
        ! noncollisional_iceformation, except that freezing fragmentation
        ! (Mode 1) is deliberately disabled here.
        !
        ! This ensures that both callbacks use the same:
        !   * Koehler/FHH activated-particle criterion;
        !   * INAS, DeMott, Daily/DCMEX and Koop mechanisms;
        !   * flat-within-dry-bin treatment of the DeMott >0.5 um reservoir;
        !   * IN/provenance-moment depletion and transfer;
        !   * full-moving versus fixed-grid ice receiving-bin treatment;
        !   * aerosol, number, ice-mass, morphology and latent-heat updates.
        !
        ! Using a wrapper rather than duplicating the implementation also
        ! prevents the two nucleation paths drifting apart in future.
        call noncollisional_iceformation( &
            npart,npartice,mwat,mbin2,mbin2_ice, &
            rhobin,nubin,kappabin,molwbin, &
            moments,parcel1%mbinedges, &
            t,p,parcel1%n_bins1,sz,sz2,sz3,parcel1%n_modes, &
            yice,rh,dt,parcel1%sce_flag,.false.,ice_nucleation_mech)

    end subroutine icenucleation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

	! ============================================================================
	! chen_and_lamb_prop
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Updates deposited ice volume and spheroidal aspect ratio for a vapour-deposition mass increment
	!>following Chen and Lamb (1994).
	!>@param[in] dm: deposited/sublimated ice-mass increment
	!>@param[in] gamma_t: Chen-Lamb temperature-dependent growth-ratio parameter
	!>@param[inout] v: deposited ice volume
	!>@param[inout] phi: particle aspect ratio c/a
	!>@param[in] dep_density: density assigned to newly deposited ice
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
        rgamma_tp2=1._wp/(gamma_t+2._wp)
        ln_vn_vo=log(v/v_old)
        phi=phi*exp((gamma_t-1._wp)*rgamma_tp2*ln_vn_vo)       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    end subroutine chen_and_lamb_prop
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	! ============================================================================
	! moving_centre
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Conservatively reassigns particles and extensive moments to the fixed bin whose mass interval
	!>contains each current representative mass using the moving-centre scheme.
	!> see Jacobson's Book, Fundamentals of Atmospheric Modelling
	!>@param[in] n_bin_mode: total number of bins represented by the state
	!>@param[in] n_bin_modew: number of bins being remapped
	!>@param[in] n_binst: number of bins per mode
	!>@param[in] n_mode: number of external modes
	!>@param[in] n_comps: number of aerosol composition components
	!>@param[in] n_moments: number of conserved moments
	!>@param[inout] npart: particle number concentration
	!>@param[inout] masses: representative particle mass in each bin
	!>@param[inout] moments: extensive conserved moments
	!>@param[inout] mbin: per-particle component and water/ice masses
	!>@param[in] mbinedges: fixed mass-bin edges for each mode
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


	! ============================================================================
	! apply_growth_bin_scheme
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Dispatches diffusional-growth remapping to the configured full-moving, moving-centre or
	!>Chen-Lamb bin scheme.
	!>@param[inout] npart: particle number concentration
	!>@param[in] mass_old: representative masses before diffusional growth
	!>@param[inout] mass_new: representative masses after growth and, when required, after remapping
	!>@param[inout] moments: extensive conserved moments
	!>@param[inout] mbin: per-particle component and water/ice masses
	subroutine apply_growth_bin_scheme(npart,mass_old,mass_new,moments,mbin)
		implicit none
		real(wp), dimension(parcel1%n_bin_modew), intent(inout) :: npart,mass_new
		real(wp), dimension(parcel1%n_bin_modew),intent(in) :: mass_old
		real(wp), dimension(parcel1%n_bin_modew,parcel1%n_comps+parcel1%imoms), &
			intent(inout) :: moments
		real(wp), dimension(parcel1%n_bin_modew, parcel1%n_comps+1), &
			intent(inout) :: mbin
	
		select case(parcel1%bin_scheme_flag)
		case(BIN_FULL_MOVING)
			! No remapping after diffusional growth.  Keep the duplicated
			! per-particle hydrometeor-mass column synchronized with the
			! accepted DVODE representative mass.
			mbin(:,n_comps+1)=mass_new
		case(BIN_MOVING_CENTRE)
			call moving_centre( parcel1%n_bin_mode, &
				parcel1%n_bin_modew, parcel1%n_bins1, &
				parcel1%n_modes, parcel1%n_comps, &
				parcel1%n_comps+parcel1%imoms, npart, &
				mass_new, moments, mbin, parcel1%mbinedges)
		case(BIN_CHEN_LAMB)
			call chen_lamb_growth_remap( parcel1%n_bin_modew, &
				parcel1%n_bins1, parcel1%n_modes, &
				parcel1%n_comps, parcel1%n_comps+parcel1%imoms, &
				npart, mass_old, mass_new, &
				moments, mbin, parcel1%mbinedges)
		case default
			error stop 'Unknown bin_scheme_flag'
		end select
	end subroutine apply_growth_bin_scheme
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! chen_lamb_int_number
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Analytically integrates the Chen-Lamb linear within-bin number distribution over a specified
	!>mass interval using a numerically stable local coordinate.
	!> Chen and Lamb (1994) binning														 
	!> Eqs. (6)–(8), with analytical transfer integrals in Eq. (15) for Chen and Lamb	 
	!> binning number integral															 
	!>@param[in] nref: number-distribution value at reference mass xref
	!>@param[in] xref: reference mass coordinate
	!>@param[in] slope: slope of the linear number distribution
	!>@param[in] a,b: lower and upper integration limits
	!>@return val: integrated particle number over [a,b]
	pure function chen_lamb_int_number(nref,xref,slope,a,b) result(val)
		implicit none
		real(wp), intent(in) :: nref,xref,slope,a,b
		real(wp) :: val
		real(wp) :: width, na
	
		if (b <= a) then
			val=0._wp
			return
		endif
		! ----------------------------------------------------------
		! Write the linear distribution relative to the LEFT
		! integration boundary:
		!
		!   n(x) = na + slope*(x-a)
		!
		! This avoids subtracting nearly equal powers of a-xref
		! and b-xref for very narrow overlaps.
		! ----------------------------------------------------------
		width=b-a
		na=nref+slope*(a-xref)
		! Integral_a^b n(x) dx
		val=width*(na + 0.5_wp*slope*width)
	
	end function chen_lamb_int_number
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! chen_lamb_int_mass
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Analytically integrates the first mass moment of the Chen-Lamb linear within-bin number
	!>distribution over a specified mass interval.
	!> Chen and Lamb (1994) binning														 
	!> Eqs. (6)–(8), with analytical transfer integrals in Eq. (15) for Chen and Lamb	 
	!> binning mass integral				    											 
	!>@param[in] nref: number-distribution value at reference mass xref
	!>@param[in] xref: reference mass coordinate
	!>@param[in] slope: slope of the linear number distribution
	!>@param[in] a,b: lower and upper integration limits
	!>@return val: integrated particle mass moment over [a,b]
	pure function chen_lamb_int_mass(nref,xref,slope,a,b) result(val)
		implicit none
		real(wp), intent(in) :: nref,xref,slope,a,b
		real(wp) :: val
		real(wp) :: width,na,dN,dM_about_a
	
		if (b <= a) then
			val=0._wp
			return
		endif
		! ----------------------------------------------------------
		! Local coordinate:
		!
		!   t = x-a
		!
		!   n(x) = na + slope*t
		!
		! with 0 <= t <= width.
		!
		! Then
		!
		!   M = integral x*n(x) dx
		!
		!     = a*N
		!       + integral t*n(t) dt
		!
		! This form is substantially more stable for tiny overlaps.
		! ----------------------------------------------------------
		width=b-a
		na=nref+slope*(a-xref)
		! Zeroth moment over this overlap
		dN=width*(na + 0.5_wp*slope*width)
		! First moment measured relative to x=a
		dM_about_a=width**2 * ( 0.5_wp*na + slope*width/3._wp)
		! Absolute first mass moment
		val=a*dN+dM_about_a
	end function chen_lamb_int_mass
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! ============================================================================
	! chen_lamb_distribution
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Constructs the positive linear Chen-Lamb within-bin number distribution that reproduces
	!>specified bin number and mass moments.
	!> Chen and Lamb (1994) binning														 
	!> Eqs. (8)–(10)  Chen and Lamb binning, calculate the density and slope, with        
	!> special cases	for when they go negative - see paper                                
	!>@param[in] x1,x2: lower and upper source-bin mass edges
	!>@param[in] N: total particle number in the source bin
	!>@param[in] M: total first mass moment in the source bin
	!>@param[out] nref: reconstructed distribution value at xref
	!>@param[out] xref: reference mass for the linear representation
	!>@param[out] slope: reconstructed distribution slope
	!>@param[out] xlo,xhi: lower and upper support over which the reconstructed distribution is positive
	subroutine chen_lamb_distribution(x1,x2,N,M,nref,xref,slope,xlo,xhi)
		implicit none
		real(wp), intent(in) :: x1,x2,N,M
		real(wp), intent(out) :: nref,xref,slope,xlo,xhi
		real(wp) :: width,x0,n0,xmean,nleft,nright, xstar,tolx,toln
	
		if (N <= tiny(1._wp)) then
			nref=0._wp
			xref=0._wp
			slope=0._wp
			xlo=x1
			xhi=x2
			return
		endif
		width=x2-x1
		if (width <= 0._wp) then
			error stop 'Chen-Lamb: non-positive bin width'
		endif
		xmean=M/N
		! The representative mean of a valid fixed bin must lie inside
		! its own boundaries.
		tolx=100._wp*epsilon(1._wp)* &
			max(abs(x1),abs(x2),abs(xmean),tiny(1._wp))
		if ((xmean < x1-tolx).or. &
			(xmean > x2+tolx)) then
			print *,'Chen-Lamb source mean outside bin'
			print *,'x1, mean, x2 = ',x1,xmean,x2
			error stop
		endif
		! Protect only against tiny floating-point excursions.
		xmean=min(max(xmean,x1),x2)
		x0=0.5_wp*(x1+x2)
		n0=N/width
		slope=12._wp*(M-x0*N)/width**3
		nleft =n0+slope*(x1-x0)
		nright=n0+slope*(x2-x0)
		
		toln=100._wp*epsilon(1._wp)* &
			max(abs(n0),abs(slope*width),tiny(1._wp))
	
		! ================================================================
		! Normal positive linear distribution
		! ================================================================
		if ((nleft >= -toln).and. &
			(nright >= -toln)) then
			nref=n0
			xref=x0
			xlo=x1
			xhi=x2
			return
		endif
		! ================================================================
		! Distribution would be negative at lower edge.
		!
		! Chen & Lamb:
		!
		!     x* = 3 M/N - 2 x2
		!
		! Represent only x* <= x <= x2 with
		!
		!     n(x)=k*(x-x*)
		!
		! The slope below is obtained directly from conservation of N.
		! ================================================================
		if (nleft < -toln) then
			xstar=3._wp*xmean-2._wp*x2
			xstar=min(max(xstar,x1),x2)
			if ((x2-xstar) <= tiny(1._wp)) then
				error stop 'Chen-Lamb lower positivity correction failed'
			endif
			nref=0._wp
			xref=xstar
			slope=2._wp*N/(x2-xstar)**2
			xlo=xstar
			xhi=x2
			return
		endif
		! ================================================================
		! Distribution would be negative at upper edge.
		!
		!     x* = 3 M/N - 2 x1
		!
		! Represent only x1 <= x <= x*.
		! ================================================================
		if (nright < -toln) then
			xstar=3._wp*xmean-2._wp*x1
			xstar=min(max(xstar,x1),x2)
			if ((xstar-x1) <= tiny(1._wp)) then
				error stop 'Chen-Lamb upper positivity correction failed'
			endif
			nref=0._wp
			xref=xstar
			slope=-2._wp*N/(xstar-x1)**2
			xlo=x1
			xhi=xstar
			return
		endif
		error stop 'Chen-Lamb positivity logic failure'
	end subroutine chen_lamb_distribution
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	! ============================================================================
	! chen_lamb_growth_remap
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Remaps a population after diffusional growth using the Chen-Lamb linear within-bin distribution
	!>and analytical number/mass transfer integrals while conserving extensive moments.
	!> The accepted mean mass change is supplied by the BMM/DVODE solution.
	!> The lower and upper sub-bin limits are shifted independently assuming
	!> dm/dt proportional to m^(1/3), following the Chen & Lamb condensation
	!> example.  A new postgrowth linear distribution is then reconstructed
	!> from number and total postgrowth mass.
	!>
	!> For ice, use of the m^(1/3) edge scaling is a BMM approximation.
	!>@param[in] n_bin_modew: number of bins being remapped
	!>@param[in] n_binst: number of fixed bins per mode
	!>@param[in] n_mode: number of external modes
	!>@param[in] n_comps: number of aerosol composition components
	!>@param[in] n_moments: number of conserved moments
	!>@param[inout] npart: particle number concentration
	!>@param[in] mass_old: representative mass before growth
	!>@param[inout] mass_new: representative mass after growth/remapping
	!>@param[inout] moments: extensive conserved moments
	!>@param[inout] mbin: per-particle component and water/ice masses
	!>@param[in] mbinedges: fixed mass-bin edges for each mode
	subroutine chen_lamb_growth_remap( n_bin_modew,n_binst,&
		n_mode,n_comps,n_moments, &
		npart,mass_old,mass_new,moments,mbin,mbinedges)
		use sce, only : qsmall2
		implicit none
		integer(i4b), intent(in) :: n_bin_modew,n_binst,&
			n_mode,n_comps,n_moments
		real(wp), dimension(n_binst+1,n_mode), intent(in) :: mbinedges
		real(wp), dimension(n_bin_modew), intent(in) :: mass_old
		real(wp), dimension(n_bin_modew), intent(inout) :: &
			npart,mass_new
		real(wp), dimension(n_bin_modew,n_moments), &
			intent(inout) :: moments
		real(wp), dimension(n_bin_modew,n_comps+1), &
			intent(inout) :: mbin
		integer(i4b) :: imode,isbin,idbin,isrc,idest,k,idlo,idhi,last_idbin
		real(wp) :: Nsrc,Mtarget, x1,x2,xmean_old,xmean_new, &
			dxmean,dx1,dx2, x1_new,x2_new, nref,xref,slope,xlo_new,xhi_new, &
			a,b,dN,dM,frac,frac_sum, Nraw,Mraw,Nmapped,Mmapped,Nfull,Mfull, &
			corrN,corrM,relN,relM, Nbefore,Nafter,Mbefore,Mafter, tolN,tolM
		real(wp), dimension(n_bin_modew) :: npart_tmp,mass_tmp
		real(wp), dimension(n_bin_modew,n_moments) :: moments_tmp
		real(wp), dimension(n_binst) :: dN_work,dM_work
		real(wp), dimension(n_moments) :: moments_before,moments_after
		! --------------------------------------------------------------
		! Scratch destination state
		! --------------------------------------------------------------
		npart_tmp=0._wp
		mass_tmp=0._wp
		moments_tmp=0._wp
		! --------------------------------------------------------------
		! Remove populations below the model's numerical number floor.
		!
		! SCE itself does not treat populations below qsmall2 as
		! physically active.  Do this before calculating the remapping
		! conservation targets so that numerical ghosts are not included
		! in Nbefore/Mbefore.
		! --------------------------------------------------------------
		do isrc=1,n_bin_modew
			if (npart(isrc) <= qsmall2) then
				npart(isrc)=0._wp
				moments(isrc,:)=0._wp
			endif	
		enddo
		! --------------------------------------------------------------
		! Conservation targets AFTER physical growth, BEFORE remapping
		! --------------------------------------------------------------
		Nbefore=sum(npart)
		Mbefore=sum(npart*mass_new)
		moments_before=sum(moments,dim=1)
		! ==============================================================
		! Treat each aerosol mode independently.
		!
		! There is no numerical transfer from one aerosol mode to another.
		! ==============================================================
		do imode=1,n_mode
			do isbin=1,n_binst
				isrc=(imode-1)*n_binst+isbin
				Nsrc=npart(isrc)
				if (Nsrc <= qsmall2) cycle	

				x1=mbinedges(isbin,imode)
				x2=mbinedges(isbin+1,imode)
				! ==============================================================
				! Chen-Lamb POSTGROWTH linear method
				!
				! The actual mean mass change has already been calculated by
				! BMM/DVODE.  The m^(1/3) law is used only to estimate the
				! relative motion of the unresolved lower and upper bin limits.
				! ==============================================================		
				xmean_old=mass_old(isrc)
				xmean_new=mass_new(isrc)
				! --------------------------------------------------------------
				! An occupied category must have a positive old mean mass.
				! --------------------------------------------------------------
				if (xmean_old <= 0._wp) then
					print *,'Chen-Lamb invalid old mean mass'
					print *,'mode/bin = ',imode,isbin
					print *,'mass_old = ',xmean_old
					error stop
				endif
				! --------------------------------------------------------------
				! Before diffusional growth, the representative mean should lie
				! within the fixed source bin.
				! --------------------------------------------------------------
				tolM=100._wp*epsilon(1._wp)* &
					max(abs(x1),abs(x2),abs(xmean_old),tiny(1._wp))
				if ((xmean_old < x1-tolM).or. (xmean_old > x2+tolM)) then
					print *,'Chen-Lamb old mean outside source bin'
					print *,'mode/bin = ',imode,isbin
					print *,'x1,mean,x2 = ',x1,xmean_old,x2
					error stop
				endif
				! --------------------------------------------------------------
				! Accepted mean growth from the BMM microphysics.
				! --------------------------------------------------------------
				dxmean=xmean_new-xmean_old
				! --------------------------------------------------------------
				! Chen & Lamb nonlinear edge-growth approximation:
				!
				!       dx/dt proportional to x^(1/3)
				!
				! therefore:
				!
				!       dx_edge =
				!           dx_mean * (x_edge/x_mean)^(1/3)
				!
				! This does NOT replace the physical BMM growth law.
				! --------------------------------------------------------------
				dx1 = dxmean * (max(x1,0._wp)/xmean_old)**onethird
				dx2 = dxmean * (max(x2,0._wp)/xmean_old)**onethird
				x1_new=x1+dx1
				x2_new=x2+dx2
				! --------------------------------------------------------------
				! Water/ice mass cannot be negative.
				!
				! In particular, for x1=0 the m^(1/3) relation gives dx1=0,
				! so the lowest boundary naturally remains at zero.
				! --------------------------------------------------------------
				x1_new=max(x1_new,0._wp)
				x2_new=max(x2_new,0._wp)
				! --------------------------------------------------------------
				! Closed upper boundary for the finite Chen-Lamb mass grid.
				!
				! If an occupied final source category undergoes positive growth,
				! the Chen-Lamb upper sub-bin edge will necessarily move beyond
				! the finite grid even when the accepted mean mass remains well
				! inside the final category.
				!
				! In that case constrain only the unresolved upper support edge
				! to the grid maximum.  The postgrowth distribution is subsequently
				! reconstructed from Nsrc and Mtarget, so zeroth- and first-moment
				! conservation are retained.
				!
				! If the accepted mean itself lies beyond the grid, however, the
				! grid is genuinely too small and execution must stop.
				! --------------------------------------------------------------
				tolM=1000._wp*epsilon(1._wp)* &
					max(abs(x1_new),abs(x2_new),abs(xmean_new), &
						abs(mbinedges(n_binst+1,imode)),tiny(1._wp))
				if (x2_new > mbinedges(n_binst+1,imode)) then
					if (isbin == n_binst .and. &
						xmean_new <= mbinedges(n_binst+1,imode)+tolM) then
						x2_new=mbinedges(n_binst+1,imode)
					else
						print *,'Chen-Lamb genuine upper-grid overflow'
						print *,'mode/source bin = ',imode,isbin
						print *,'Nsrc = ',Nsrc
						print *,'old/new mean = ',xmean_old,xmean_new
						print *,'new limits before boundary = ',x1+dx1,x2+dx2
						print *,'grid maximum = ',mbinedges(n_binst+1,imode)
						error stop
					endif
				endif
				! --------------------------------------------------------------
				! The postgrowth interval must remain ordered.
				! --------------------------------------------------------------
				if (x2_new <= x1_new) then
					print *,'Chen-Lamb postgrowth bin collapse'
					print *,'mode/bin = ',imode,isbin
					print *,'old limits = ',x1,x2
					print *,'new limits = ',x1_new,x2_new
					print *,'old/new means = ',xmean_old,xmean_new
					print *,'dxmean = ',dxmean
					print *,'dx1,dx2 = ',dx1,dx2
					error stop
				endif
				! --------------------------------------------------------------
				! The accepted postgrowth mean must lie within the postgrowth
				! limits.  If this fails, the m^(1/3) edge approximation has
				! become inconsistent with the DVODE mean growth over this dt.
				! --------------------------------------------------------------
				tolM=100._wp*epsilon(1._wp)* max(abs(x1_new),abs(x2_new), &
						abs(xmean_new),tiny(1._wp))
				if ((xmean_new < x1_new-tolM).or. (xmean_new > x2_new+tolM)) then
					print *,'Chen-Lamb new mean outside postgrowth limits'
					print *,'mode/bin = ',imode,isbin
					print *,'new limits = ',x1_new,x2_new
					print *,'old/new means = ',xmean_old,xmean_new
					print *,'dxmean = ',dxmean
					print *,'dx1,dx2 = ',dx1,dx2
					error stop
				endif
				! --------------------------------------------------------------
				! Postgrowth first moment supplied by DVODE.
				! --------------------------------------------------------------
				Mtarget=Nsrc*xmean_new
				! --------------------------------------------------------------
				! Reconstruct the POSTGROWTH linear distribution from:
				!
				!       N = Nsrc
				!       M = Nsrc*xmean_new
				!       limits = x1_new,x2_new
				!
				! chen_lamb_distribution also applies the positivity correction
				! if the unconstrained line becomes negative at an edge.
				! --------------------------------------------------------------
				call chen_lamb_distribution(x1_new,x2_new,Nsrc,Mtarget, &
					nref,xref,slope,xlo_new,xhi_new)
					
				! --------------------------------------------------------------
				! Verify that the reconstructed POSTGROWTH distribution contains
				! exactly the intended number and first mass moment before
				! remapping it onto the fixed grid.
				! --------------------------------------------------------------
				Nfull=chen_lamb_int_number(nref,xref,slope,xlo_new,xhi_new)
				Mfull=chen_lamb_int_mass(nref,xref,slope,xlo_new,xhi_new)
				tolN=1.e-10_wp*max(abs(Nsrc),1.e-300_wp)
				if (abs(Nfull-Nsrc) > tolN) then
					print *,'Chen-Lamb postgrowth reconstruction number error'
					print *,'mode/bin = ',imode,isbin
					print *,'Nsrc,Nfull = ',Nsrc,Nfull
					print *,'postgrowth limits = ',x1_new,x2_new
					print *,'positive support = ',xlo_new,xhi_new
					error stop
				endif
				tolM=1.e-10_wp* max(abs(Mtarget),tiny(1._wp))
				if (abs(Mfull-Mtarget) > tolM) then
					print *,'Chen-Lamb postgrowth reconstruction mass error'
					print *,'mode/bin = ',imode,isbin
					print *,'Mtarget,Mfull = ',Mtarget,Mfull
					print *,'postgrowth limits = ',x1_new,x2_new
					print *,'positive support = ',xlo_new,xhi_new
					error stop
				endif
				
    			! ======================================================
				! Local mass-grid boundary protection
				! ======================================================
				tolM=1000._wp*epsilon(1._wp)* &
					max(abs(x1_new),abs(x2_new), abs(xlo_new),abs(xhi_new), &
						abs(xmean_new), abs(x2_new-x1_new), tiny(1._wp))
				if (xlo_new < mbinedges(1,imode)-tolM) then
					print *,'Chen-Lamb lower-grid overflow'
					print *,'mode/source bin = ',imode,isbin
					print *,'support = ',xlo_new,xhi_new
					print *,'grid = ', mbinedges(1,imode), &
						mbinedges(n_binst+1,imode)
					error stop
				endif
				if (xhi_new > mbinedges(n_binst+1,imode)+tolM) then
					print *,'Chen-Lamb upper-grid overflow'
					print *,'mode/source bin = ',imode,isbin
					print *,'Nsrc = ',Nsrc
					print *,'mode N = ', &
						sum(npart((imode-1)*n_binst+1:imode*n_binst))
					print *,'Nsrc/modeN = ', &
						Nsrc/max(sum(npart((imode-1)*n_binst+1:imode*n_binst)), &
								 tiny(1._wp))
					print *,'old limits = ',x1,x2
					print *,'old/new mean = ',xmean_old,xmean_new
					print *,'dxmean = ',dxmean
					print *,'dx1,dx2 = ',dx1,dx2
					print *,'new limits = ',x1_new,x2_new
					print *,'positive support = ',xlo_new,xhi_new
					print *,'grid = ', &
						mbinedges(1,imode), &
						mbinedges(n_binst+1,imode)
					print *,'overflow amount = ', &
						xhi_new-mbinedges(n_binst+1,imode)
					print *,'component masses = ', &
						mbin(isrc,1:n_comps)
					print *,'current water/ice mass = ', &
						mbin(isrc,n_comps+1)
					print *,'moments = ',moments(isrc,:)
					error stop
				endif
				! Tiny floating-point excursions only
				xlo_new=max(xlo_new,mbinedges(1,imode))
				xhi_new=min(xhi_new,mbinedges(n_binst+1,imode))
				! ==============================================================
				! Find the first and last fixed bins overlapped by the
				! postgrowth distribution.
				!
				! Start from the source bin.  Diffusional growth normally moves
				! the distribution by only a small number of fixed categories,
				! so this avoids scanning the entire grid.
				! ==============================================================
				idlo=isbin
				! Lower support may have moved to smaller bins.
				do while (idlo > 1)
					if (xlo_new >= mbinedges(idlo,imode)) exit
					idlo=idlo-1
				enddo
				! Or the entire lower support may have moved to larger bins.
				do while (idlo < n_binst)
					if (xlo_new < mbinedges(idlo+1,imode)) exit
					idlo=idlo+1
				enddo
				idhi=isbin
				! Upper support may have moved to larger bins.
				do while (idhi < n_binst)
					if (xhi_new <= mbinedges(idhi+1,imode)) exit
					idhi=idhi+1
				enddo
				! Or the entire upper support may have moved to smaller bins.
				do while (idhi > 1)
					if (xhi_new > mbinedges(idhi,imode)) exit
					idhi=idhi-1
				enddo
				! This should never occur for a valid ordered support.
				if (idhi < idlo) then
					print *,'Chen-Lamb invalid destination range'
					print *,'mode/source bin = ',imode,isbin
					print *,'idlo,idhi = ',idlo,idhi
					print *,'support = ',xlo_new,xhi_new
					error stop
				endif
				
				! ==============================================================
				! First pass:
				!
				! Calculate the analytical Chen-Lamb transfer into each
				! overlapping destination category.
				!
				! Keep these locally so that tiny floating-point integration
				! errors can be normalized source-by-source before modifying
				! the destination spectrum.
				! ==============================================================
				Nraw=0._wp
				Mraw=0._wp
				dN_work(idlo:idhi)=0._wp
				dM_work(idlo:idhi)=0._wp
				
				do idbin=idlo,idhi
					a=max(xlo_new,mbinedges(idbin,imode))
					b=min(xhi_new,mbinedges(idbin+1,imode))
					if (b <= a) cycle
				
					dN=chen_lamb_int_number(nref,xref,slope,a,b)
					dM=chen_lamb_int_mass(nref,xref,slope,a,b)
					! ----------------------------------------------------------
					! Allow only genuine floating-point negative noise.
					! Scale tolerance to THIS source category, not to 1.
					! ----------------------------------------------------------
					tolN=1.e-12_wp*max(abs(Nsrc),1.e-300_wp)
					if (dN < -tolN) then
						print *,'Chen-Lamb negative transferred number'
						print *,'mode/source/destination = ',imode,isbin,idbin
						print *,'dN = ',dN
						error stop
					endif
					if (dN < 0._wp) dN=0._wp
					tolM=1.e-12_wp* max(abs(Mtarget),1.e-300_wp)
					if (dM < -tolM) then
						print *,'Chen-Lamb negative transferred mass'
						print *,'mode/source/destination = ',imode,isbin,idbin
						print *,'dM = ',dM
						error stop
					endif
					if (dM < 0._wp) dM=0._wp
					dN_work(idbin)=dN
					dM_work(idbin)=dM
					Nraw=Nraw+dN
					Mraw=Mraw+dM
				enddo
				if (Nraw <= 0._wp) then				
					print *,'Chen-Lamb zero mapped particle number'
					print *,'mode/source bin = ',imode,isbin
					print *,'Nsrc = ',Nsrc
					print *,'support = ',xlo_new,xhi_new
					error stop
				endif
				relN=abs(Nraw-Nsrc) / max(abs(Nsrc),1.e-300_wp)
				if (abs(Mtarget) > 1.e-300_wp) then
					relM=abs(Mraw-Mtarget) / abs(Mtarget)
				else
					relM=abs(Mraw-Mtarget)
				endif
				if (relN > 1.e-6_wp) then
					print *,'Chen-Lamb raw number integration error'
					print *,'mode/source bin = ',imode,isbin
					print *,'Nsrc,Nraw = ',Nsrc,Nraw
					print *,'relative error = ',relN
					error stop
				endif
				if (relM > 1.e-6_wp) then
					print *,'Chen-Lamb raw mass integration error'
					print *,'mode/source bin = ',imode,isbin
					print *,'Mtarget,Mraw = ',Mtarget,Mraw
					print *,'relative error = ',relM
					error stop
				endif
				corrN=Nsrc/Nraw				
				if (abs(Mtarget) > 1.e-300_wp) then
					if (abs(Mraw) <= 1.e-300_wp) then
						print *,'Chen-Lamb zero mapped mass'
						print *,'mode/source bin = ',imode,isbin
						error stop
					endif
					corrM=Mtarget/Mraw
				else
					corrM=1._wp
				endif
				! ------------------------------------------------------------
				! Find the highest-index destination bin that actually
				! receives a contribution from this source bin.  Its
				! number-fraction will be set as the *complement* of all
				! the others below, so that the fractions used to split
				! the auxiliary moments sum to exactly 1.0 in floating
				! point (a true partition of unity), rather than to
				! 1.0 +/- a few ULP.
				!
				! Without this, each frac=dN/Nsrc is rounded independently.
				! Individually these roundoffs are tiny (~1e-16), but they
				! do not cancel, and moments(isrc,:) can be many orders of
				! magnitude larger than the isrc,k component being
				! conserved (e.g. isrc's total mass moment is O(1) while a
				! trace chemical/auxiliary moment k is O(1e-8)). The
				! absolute roundoff picked up from partitioning the large
				! moments then shows up as a large *relative* error on the
				! small ones once everything is summed over all bins -
				! this is exactly the "moment = 9" failure mode reported.
				! ------------------------------------------------------------
				last_idbin=idlo
				do idbin=idlo,idhi
					if (dN_work(idbin) > 0._wp) last_idbin=idbin
				enddo

				Nmapped=0._wp
				Mmapped=0._wp
				frac_sum=0._wp
				do idbin=idlo,idhi
					if (dN_work(idbin) <= 0._wp) cycle
				
					! ----------------------------------------------------------
					! Enforce exact source-wise zeroth and first moment
					! conservation to floating-point precision.
					! ----------------------------------------------------------
					dN=dN_work(idbin)*corrN
					dM=dM_work(idbin)*corrM
					idest=(imode-1)*n_binst+idbin
					npart_tmp(idest)= npart_tmp(idest)+dN
					mass_tmp(idest)= mass_tmp(idest)+dM
					! ----------------------------------------------------------
					! Auxiliary moments follow the normalized particle-number
					! fraction.
					!
					! Because corrN=Nsrc/Nraw:
					!
					!       dN/Nsrc = dN_raw/Nraw
					!
					! so these fractions form a partition of unity for each
					! source category *mathematically*. Force that to also
					! be true numerically by closing the last contributing
					! bin's share with the exact remainder instead of another
					! independently-rounded division.
					! ----------------------------------------------------------
					if (idbin == last_idbin) then
						frac=1._wp-frac_sum
						if (frac < -100._wp*epsilon(1._wp)) then
							print *,'Chen-Lamb negative closing fraction'
							print *,'frac,frac_sum = ',frac,frac_sum
							error stop
						endif
						frac=max(frac,0._wp)
					else
						frac=dN/Nsrc
						frac_sum=frac_sum+frac
					endif
					moments_tmp(idest,:)= moments_tmp(idest,:) + moments(isrc,:)*frac
					Nmapped=Nmapped+dN
					Mmapped=Mmapped+dM
				enddo
				
				! ======================================================
				! Source-bin conservation checks
				! ======================================================
				tolN=1.e-10_wp*max(abs(Nsrc),1.e-300_wp)
				if (abs(Nmapped-Nsrc) > tolN) then
					print *,'Chen-Lamb source number error'
					print *,'mode/bin = ',imode,isbin
					print *,'original,mapped = ',Nsrc,Nmapped
					print *,'postgrowth limits = ',x1_new,x2_new
					print *,'positive support = ',xlo_new,xhi_new
					error stop
				endif
				tolM=1.e-10_wp*max(abs(Mtarget),tiny(1._wp))
				if (abs(Mmapped-Mtarget) > tolM) then
					print *,'Chen-Lamb source mass error'
					print *,'mode/bin = ',imode,isbin
					print *,'original,mapped = ',Mtarget,Mmapped
					print *,'postgrowth limits = ',x1_new,x2_new
					print *,'positive support = ',xlo_new,xhi_new
					error stop
				endif
			enddo
		enddo
		! ==============================================================
		! Replace phase state by the newly remapped state
		! ==============================================================
		do imode=1,n_mode
			do idbin=1,n_binst
				idest=(imode-1)*n_binst+idbin
				if (npart_tmp(idest) > tiny(1._wp)) then
					npart(idest)=npart_tmp(idest)
					mass_new(idest)= mass_tmp(idest)/npart(idest)
					moments(idest,:)= moments_tmp(idest,:)
					! Aerosol / chemical component mass per particle
					mbin(idest,1:n_comps)= &
						moments(idest,1:n_comps)/npart(idest)
	
					! Current water/ice mass per particle
					mbin(idest,n_comps+1)=mass_new(idest)
					! The destination mean must lie inside the destination
					! fixed mass interval.
					if ((mass_new(idest) < mbinedges(idbin,imode)).or. &
						(mass_new(idest) > mbinedges(idbin+1,imode))) then
						print *,'Chen-Lamb destination mean outside bin'
						print *,'mode/bin = ',imode,idbin
						print *,'lower,mean,upper = ', &
							mbinedges(idbin,imode), &
							mass_new(idest), mbinedges(idbin+1,imode)
						error stop
					endif
				else
					npart(idest)=0._wp
					moments(idest,:)=0._wp
					! Empty bins need a harmless representative mass.
					mass_new(idest)=0.5_wp*( &
						mbinedges(idbin,imode)+ &
						mbinedges(idbin+1,imode))
	
					mbin(idest,1:n_comps)=0._wp
					mbin(idest,n_comps+1)=mass_new(idest)
				endif
			enddo
		enddo
		! ==============================================================
		! Whole-phase conservation checks
		! ==============================================================
		Nafter=sum(npart)
		Mafter=sum(npart*mass_new)
		moments_after=sum(moments,dim=1)
	
		if (abs(Nafter-Nbefore) > &
			1.e-10_wp*max(Nbefore,1._wp)) then
			print *,'Chen-Lamb total number error = ', Nafter-Nbefore
			error stop
		endif
	
		if (abs(Mafter-Mbefore) > &
			1.e-10_wp*max(abs(Mbefore),tiny(1._wp))) then
			print *,'Chen-Lamb total mass error = ', Mafter-Mbefore
			error stop
		endif
		do k=1,n_moments
			! --------------------------------------------------------------
			! Auxiliary-moment conservation.
			!
			! Auxiliary moments are partitioned using normalized particle-
			! number fractions, with the final contributing fraction set as
			! the complement of all previous fractions.  They should
			! therefore conserve to normal floating-point accumulation error.
			! --------------------------------------------------------------
			tolM=max( &
				1.e-10_wp*max(abs(moments_before(k)),abs(moments_after(k))), &
				tiny(1._wp))
		
			if (abs(moments_after(k)-moments_before(k)) > tolM) then
				print *,'Chen-Lamb moment conservation error'
				print *,'moment = ',k
				print *,'before,after = ', moments_before(k),moments_after(k)
				print *,'absolute error = ', moments_after(k)-moments_before(k)
				print *,'relative error = ', &
					abs(moments_after(k)-moments_before(k)) / &
					max(abs(moments_before(k)),tiny(1._wp))
				error stop
			endif
		enddo
	end subroutine chen_lamb_growth_remap
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! chen_and_lamb_anc
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates the Chen-Lamb (1994) temperature/supersaturation-dependent 
	!> growth-ratio parameter and
	!>density of newly deposited ice.
	!>@param[in] t: temperature
	!>@param[in] qv: ambient water-vapour mixing ratio
	!>@param[in] qvsat: saturation water-vapour mixing ratio over ice
	!>@param[in] rhoa: air density
	!>@param[inout] gamma_t: Chen-Lamb growth-ratio parameter
	!>@param[inout] dep_density: density of newly deposited ice
    subroutine chen_and_lamb_anc(t,qv,qvsat,rhoa,gamma_t, dep_density)
        implicit none
        real(wp), intent(in) :: t,qv,qvsat,rhoa
        real(wp), intent(inout) :: gamma_t,dep_density

        real(wp) :: delta_rho,t1
        integer(i4b) :: i
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! calculate the inherent growth ratio - this is from a 17th order polynomial
        gamma_t=0._wp
        t1=min(max(t,243.15),273.15) ! range of fit
        do i=1,n_cl
            gamma_t=gamma_t+((t1-gam_mu_cl(1))/gam_mu_cl(2))**(n_cl-i)*gam_cl(i)
        enddo
        gamma_t=10._wp**gamma_t
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! equation 42 from Chen and Lamb (1994, JAS: The Theoretical Basis for 
        !   Parameterisation of Ice Crystal Habits)
        delta_rho=(qv-qvsat)*rhoa*1000._wp ! g/m^3
        dep_density=rhoice*exp(-3._wp*max(delta_rho-0.05_wp,0._wp)/gamma_t)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine chen_and_lamb_anc
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! chen_and_lamb_cap_fac
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates the ratio of spheroidal ice capacitance to the radius of an equal-volume sphere
	!>following Chen and Lamb (1994).
	!>@param[in] phi: spheroid aspect ratio c/a
	!>@return chen_and_lamb_cap_fac: capacitance correction factor relative to an equal-volume sphere
    function chen_and_lamb_cap_fac(phi)
        implicit none
        real(wp), intent(in) :: phi
        real(wp) :: chen_and_lamb_cap_fac
        real(wp) :: ecc, phi1

        ! convert between R and a - derived from equating volume of sphere to 
        ! volume of spheroid and taking the ratio of a / r
        phi1=max(phi,1.e-8_wp)
        
        ! factor to convert between a and capacitance
        if(phi1<0.99_wp) then
			! Oblate spheroid (plate)
			! Chen and Lamb (1994), Eq. 39:
			! C = a*ecc/asin(ecc)
			! R = a*phi^(1/3)
            ! see equation 39 of Chen and Lamb (1994)
            ecc=sqrt(1._wp-phi1**2)
            chen_and_lamb_cap_fac=ecc/asin(ecc)*phi1**(-onethird)
        elseif(phi1>1.01_wp) then
			! Prolate spheroid (column)
			! Chen and Lamb (1994), Eq. 40:
			! C = c*ecc/log[(1+ecc)*phi]
			! R = a*phi^(1/3), c=a*phi
            ! see equation 40 of Chen and Lamb (1994)
			ecc=sqrt(max(1._wp-phi1**(-2._wp),0._wp))
	
			chen_and_lamb_cap_fac = &
				phi1**(twothirds) * ecc / log((1._wp+ecc)*phi1)
        else
	        ! Spherical limit
            chen_and_lamb_cap_fac=1._wp
        endif
        
    end function chen_and_lamb_cap_fac
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
	! ============================================================================
	! mass_balance
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates total parcel water mass mixing ratio from vapour, liquid and optionally ice for
	!>mass-conservation diagnostics/correction.
	!>@param[in] neq: length of the warm-parcel state vector
	!>@param[in] neqice: length of the ice state vector
	!>@param[in] y: warm-parcel state vector
	!>@param[in] yice: ice mass state vector
	!>@param[in] npart: liquid-particle number concentrations
	!>@param[in] npartice: ice-particle number concentrations
	!>@param[inout] mass1: returned total water mass mixing ratio
	!>@param[in] n_bin_modew: number of liquid/warm bins
	!>@param[in] irh,ite,ipr: indices of relative humidity, temperature and pressure in y
	!>@param[in] ice_flag: switch indicating whether ice mass is included
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
       
       
       
       
	! ============================================================================
	! adjust_relative_humidity
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Adjusts vapour mixing ratio and relative humidity so that total water matches a target mass
	!>after a microphysical step.
	!>@param[in] mass1: target total-water mass before the step
	!>@param[in] mass2: calculated total-water mass after the step
	!>@param[inout] vapour_mass: vapour mixing ratio corrected for the mass imbalance
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[inout] rh: relative humidity corrected to the adjusted vapour mixing ratio
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


	! ============================================================================
	! ice_particle_properties_from_moments
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Recovers mean ice aspect ratio, monomer number, depositional density and rime mass from the
	!>conserved ice moments.
	!>@param[in] yice: ice mass per representative particle
	!>@param[in] npartice: ice-particle number concentration
	!>@param[in] moments: conserved liquid and ice moments
	!>@param[in] n_bin_modew: number of ice categories corresponding to the liquid grid
	!>@param[in] n_comps: number of aerosol composition components preceding the ice-property moments
	!>@param[out] phi: mean monomer aspect ratio
	!>@param[out] nump: mean monomer number per aggregate
	!>@param[out] rhoi: mean density of deposited/unrimed ice
	!>@param[out] rime: rime mass per representative particle
	subroutine ice_particle_properties_from_moments( &
		yice,npartice,moments,n_bin_modew,n_comps, &
		phi,nump,rhoi,rime)
		implicit none
		integer(i4b), intent(in) :: n_bin_modew,n_comps
		real(wp), dimension(:), intent(in) :: yice,npartice
		real(wp), dimension(:,:), intent(in) :: moments
		real(wp), dimension(n_bin_modew), intent(out) :: &
			phi,nump,rhoi,rime
		integer(i4b) :: i,ii
		real(wp), parameter :: small=1.e-60_wp
		real(wp) :: monomer_moment
		real(wp) :: volume_moment
		real(wp) :: rime_moment
		real(wp) :: unrimed_mass_moment
	
		do i=1,n_bin_modew
			ii=n_bin_modew+i
	
			if (npartice(i) > small) then
	
				! -------------------------------------------------------------
				! Aspect ratio
				!
				! n_comps+1 = phi-weighted monomer moment
				! n_comps+2 = monomer-number moment
				! -------------------------------------------------------------
				monomer_moment=moments(ii,n_comps+2)
				if (monomer_moment > small) then
					phi(i)=moments(ii,n_comps+1) / monomer_moment
					nump(i)=monomer_moment / npartice(i)
				else
					phi(i)=1._wp
					nump(i)=1._wp
				endif
				phi(i)=max(phi(i),1.e-6_wp)
				! Mean monomer count cannot physically be less than one
				nump(i)=max(nump(i),1._wp)
				! -------------------------------------------------------------
				! Rime mass
				!
				! n_comps+4 is the total rime-mass moment.
				! Convert to rime mass per representative particle.
				! -------------------------------------------------------------
				rime_moment=max(moments(ii,n_comps+4),0._wp)
				rime(i)=rime_moment/npartice(i)
				! Numerical protection
				rime(i)=min(rime(i),max(yice(i),0._wp))
				! -------------------------------------------------------------
				! Depositional-ice density
				!
				! n_comps+3 is the total deposited-ice volume moment.
				! rhoi = unrimed deposited mass / deposited volume
				! -------------------------------------------------------------
				volume_moment=moments(ii,n_comps+3)
				unrimed_mass_moment = &
					npartice(i)*yice(i)-rime_moment
	
				if ((volume_moment > small).and. &
					(unrimed_mass_moment > small)) then
					rhoi(i)=unrimed_mass_moment/volume_moment
				else
					rhoi(i)=rhoice
				endif
				if (isnan(rhoi(i)).or.(rhoi(i) <= 0._wp)) then
					rhoi(i)=rhoice
				endif
			else
				phi(i)=1._wp
				nump(i)=1._wp
				rhoi(i)=rhoice
				rime(i)=0._wp
			endif
		enddo
	end subroutine ice_particle_properties_from_moments
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! ice_vapour_growth_properties
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Calculates air/vapour quantities needed for ice deposition, updates Chen-Lamb deposition
	!>properties and reconstructs current ice-particle properties from moments.
	!>@param[inout] rhoa: air density
	!>@param[inout] qvsat: saturation vapour mixing ratio over ice
	!>@param[inout] qv: ambient vapour mixing ratio
	!>@param[in] neqice: length of the ice state vector
	!>@param[in] yice: ice mass state vector
	!>@param[in] n_bin_modew: number of ice categories
	!>@param[inout] npartice: ice-particle number concentration
	!>@param[inout] phi,nump,rhoi,rime: reconstructed mean ice properties
	!>@param[in] n_bin_mode: total liquid-plus-ice bin count
	!>@param[in] n_moms: number of conserved moments
	!>@param[in] n_comps: number of aerosol composition components
	!>@param[in] moments: conserved particle moments
	!>@param[inout] gamma_t: Chen-Lamb growth-ratio parameter
	!>@param[inout] dep_density: density assigned to newly deposited ice
	!>@param[in] t: temperature
	!>@param[in] p: pressure
	!>@param[in] rhi: relative humidity used to derive ambient vapour mixing ratio
    subroutine ice_vapour_growth_properties(rhoa,qvsat,qv,&
        neqice,yice,n_bin_modew,npartice,phi,nump,rhoi, rime, &
        n_bin_mode,n_moms,n_comps,moments, &
        gamma_t,dep_density,t,p,rhi)
    
        implicit none
        real(wp), intent(in) :: t,p,rhi
        real(wp), intent(inout) :: rhoa, qvsat,qv,gamma_t,dep_density
        integer(i4b), intent(in) :: neqice,n_bin_modew,n_bin_mode,n_moms,n_comps
        real(wp), dimension(neqice), intent(in) :: yice
        real(wp), dimension(n_bin_modew), intent(inout) :: npartice, phi, nump, rhoi, rime
        real(wp), dimension(n_bin_mode,n_moms), intent(in) :: moments
        
        
        rhoa = p / (ra*t)                             ! density of air
        qv = rhi*eps1*svp_liq(t) / (p-svp_liq(t))        ! qv_sat
        qvsat = (eps1*svp_ice(t) / (p-svp_ice(t)))

        call chen_and_lamb_anc(t, &
                qv,qvsat,rhoa,gamma_t, dep_density) ! set the deposition density
        
		call ice_particle_properties_from_moments( &
			yice(1:n_bin_modew), npartice, moments, n_bin_modew, n_comps, &
			phi,nump,rhoi,rime)
    
    end subroutine ice_vapour_growth_properties         
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
 
	! ============================================================================
	! update_volume_and_shape
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Updates deposited-ice volume and phi-weighted moments after vapour deposition/sublimation using
	!>the Chen-Lamb volume-ratio relation.
	!>@param[in] n_bin_modew: number of ice categories
	!>@param[in] n_bin_mode: total liquid-plus-ice bin count
	!>@param[in] n_moments: number of conserved moments
	!>@param[in] n_comps: number of aerosol composition components
	!>@param[inout] momtemp: temporary deposited-volume moment array
	!>@param[inout] moments: conserved moments to update
	!>@param[in] neqice: length of the ice mass state vector
	!>@param[in] yice: new ice mass state
	!>@param[in] yoldice: old ice mass state
	!>@param[in] gamma_t: Chen-Lamb growth-ratio parameter
	!>@param[in] dep_density: density of newly deposited ice
	!>@param[in] npartice: ice-particle number concentration
	subroutine update_volume_and_shape( n_bin_modew,n_bin_mode,n_moments,n_comps, &
		momtemp,moments,neqice,yice,yoldice, gamma_t,dep_density,npartice)
		implicit none
		integer(i4b), intent(in) :: &
			n_bin_modew,n_bin_mode,n_moments,n_comps,neqice
		real(wp), intent(in) :: gamma_t,dep_density
		real(wp), dimension(neqice), intent(in) :: yice,yoldice
		real(wp), dimension(n_bin_mode,n_moments), intent(inout) :: moments
		real(wp), dimension(n_bin_mode), intent(inout) :: momtemp
		real(wp), dimension(n_bin_modew), intent(in) :: npartice
	
		integer(i4b) :: i,ii
		real(wp) :: vold,vnew,vtrial,dv,phi_exponent,monomer_moment,tolV
		real(wp), parameter :: small=1.e-60_wp
	
		! ----------------------------------------------------------
		! Chen-Lamb aspect-ratio exponent
		!
		! phi_new/phi_old =
		!     (V_new/V_old)^((gamma-1)/(gamma+2))
		! ----------------------------------------------------------
		if (abs(gamma_t+2._wp) <= small) then
			error stop 'update_volume_and_shape: gamma_t + 2 = 0'
		endif
		if (dep_density <= 0._wp) then
			error stop 'update_volume_and_shape: dep_density <= 0'
		endif
		phi_exponent = (gamma_t-1._wp)/(gamma_t+2._wp)
		do i=1,n_bin_modew
			ii=n_bin_modew+i
			! ------------------------------------------------------
			! Empty ice category.
			!
			! No extensive volume or phi-weighted moment should
			! remain attached to a category containing no particles.
			! ------------------------------------------------------
			if (npartice(i) <= small) then
				momtemp(i)=0._wp
				moments(ii,n_comps+1)=0._wp
				moments(ii,n_comps+3)=0._wp
				cycle
			endif
			! ------------------------------------------------------
			! Existing deposited-ice volume moment.
			! ------------------------------------------------------
			vold=max(moments(ii,n_comps+3),0._wp)
			! ------------------------------------------------------
			! Change in deposited volume caused by vapour
			! deposition/sublimation during this timestep.
			!
			! yice and yoldice are mass per representative particle;
			! multiplying by npartice gives the mass mixing-ratio
			! increment for this category.
			! ------------------------------------------------------
			dv = (yice(i)-yoldice(i))* npartice(i)/dep_density
			vtrial=vold+dv
			! ------------------------------------------------------
			! A substantially negative volume indicates an
			! inconsistency rather than roundoff.
			! ------------------------------------------------------
			tolV=1000._wp*epsilon(1._wp)* max(abs(vold),abs(dv),1.e-300_wp)
			if (vtrial < -tolV) then
				print *,'Negative ice volume in update_volume_and_shape'
				print *,'bin = ',i
				print *,'npartice = ',npartice(i)
				print *,'yold,ynew = ',yoldice(i),yice(i)
				print *,'vold,dv,vnew = ',vold,dv,vtrial
				error stop
			endif
			! Tiny negative excursion is only floating-point noise.
			vnew=max(vtrial,0._wp)
			momtemp(i)=vnew
			! ------------------------------------------------------
			! Update the phi-weighted monomer moment.
			!
			! Only use the Chen-Lamb volume-ratio expression when
			! both old and new deposited volumes are well defined.
			! ------------------------------------------------------
			if ((vold > small).and.(vnew > small)) then
				moments(ii,n_comps+1) = moments(ii,n_comps+1) * &
					exp(phi_exponent*log(vnew/vold))
			else
				! --------------------------------------------------
				! No meaningful old/new volume ratio exists.
				!
				! Set phi=1.  Since
				!
				! phi = M_phi / M_nmon
				!
				! this requires M_phi=M_nmon, NOT M_phi=npartice.
				! --------------------------------------------------
				monomer_moment=max(moments(ii,n_comps+2), npartice(i))
				moments(ii,n_comps+1)=monomer_moment
			endif
			! ------------------------------------------------------
			! Save non-negative deposited-volume moment.
			! ------------------------------------------------------
			moments(ii,n_comps+3)=vnew
		enddo
	end subroutine update_volume_and_shape
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! reduce_rime
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Reduces rime mass in proportion to total particle mass lost during sublimation/evaporation while
	!>preventing negative rime mass.
	!>@param[in] ipart: number of ice categories
	!>@param[in] massnew: new total particle mass
	!>@param[in] massold: old total particle mass
	!>@param[inout] rimemass: rime mass per representative particle
    subroutine reduce_rime(ipart, massnew, massold, rimemass)
        implicit none
        integer(i4b), intent(in) :: ipart
        real(wp), intent(in), dimension(ipart) :: massold, massnew
        real(wp), intent(inout), dimension(ipart) :: rimemass
        
		where (massold > tiny(1._wp))		
			rimemass = max(rimemass,0._wp) * &
				max(min(massnew,massold),0._wp) / massold
		elsewhere
			rimemass = 0._wp
		end where
    end subroutine reduce_rime
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
                
	! ============================================================================
	! bin_microphysics
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Advances one BMM microphysics timestep, including condensational/depositional growth,
	!>nucleation, bin remapping, entrainment and water-mass correction.
	!>@param[in] func1: warm-parcel ODE tendency callback
	!>@param[in] func2: ice-parcel ODE tendency callback
	!>@param[in] func3: moving-grid/non-fixed ice-nucleation callback
	!>@param[in] func4: fixed-grid/non-collisional ice-formation callback
    subroutine bin_microphysics(func1,func2,func3,func4)
    use numerics_type
    use numerics, only : zeroin, dvode, fmin
    use sce, only : qsmall
    implicit none
    real(wp) :: mass1, mass2, deltam, vapour_mass, liquid_mass, x1,x2 , cpm, &
        var, dummy, gamma_t, dep_density, rhoa, qv, qvsat, wv, &
    	ql, qtot_m,qtot, eps2=1.e-4_wp,hmin=0.0_wp,htry=1.e-1_wp
    real(wp), dimension(parcel1%n_bin_modew) :: stk, vd, impaction_time, loss_rate
    integer(i4b) :: iloc, i
    
	! ============================================================================
	! func1
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Interface-compatible wrapper for the warm-parcel ODE tendency routine used by bin_microphysics.
	!>@param[inout] neq: number of ODE equations
	!>@param[inout] tt: ODE integration time
	!>@param[inout] y: solution vector
	!>@param[inout] ydot: derivative vector
	!>@param[inout] rpar: real ODE workspace
	!>@param[inout] ipar: integer ODE workspace
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
	! ============================================================================
	! func2
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Interface-compatible wrapper for the ice-parcel ODE tendency routine used by bin_microphysics.
	!>@param[inout] neq: number of ODE equations
	!>@param[inout] tt: ODE integration time
	!>@param[inout] y: solution vector
	!>@param[inout] ydot: derivative vector
	!>@param[inout] rpar: real ODE workspace
	!>@param[inout] ipar: integer ODE workspace
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
	! ============================================================================
	! func3
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Interface-compatible ice-nucleation callback used by bin_microphysics for the non-fixed/moving
	!>particle representation.
	!>@param[inout] npart: liquid-particle number concentration
	!>@param[inout] npartice: ice-particle number concentration
	!>@param[in] mwat: liquid-water mass per representative particle
	!>@param[in] mbin2: aerosol-component masses in liquid particles
	!>@param[inout] mbin2_ice: aerosol-component and ice masses in ice particles
	!>@param[in] rhobin,nubin,kappabin,molwbin: aerosol-component thermodynamic properties
	!>@param[inout] moments: conserved liquid/ice moments
	!>@param[inout] t: parcel temperature
	!>@param[in] p: pressure
	!>@param[in] sz: number of aerosol composition components
	!>@param[in] sz2: number of liquid/warm bins
	!>@param[in] sz3: number of conserved moments
	!>@param[inout] yice: ice mass per representative particle
	!>@param[inout] rh: relative humidity
	!>@param[in] dt: timestep
    interface
        subroutine func3(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin,moments, &
                         t,p,sz,sz2,sz3,yice,rh,dt) 
            use numerics_type
            implicit none
            real(wp), intent(inout) :: t
            real(wp), intent(in) :: p,dt
            real(wp), dimension(sz2), intent(inout) :: npart,npartice
            real(wp), dimension(sz2), intent(in) :: mwat
            real(wp), dimension(sz2,sz), intent(in) :: mbin2, &
                                                  rhobin,nubin,kappabin,molwbin
            real(wp), dimension(2*sz2,sz3), intent(inout) :: moments
            integer(i4b), intent(in) :: sz,sz2, sz3
            real(wp), dimension(sz2,sz+1), intent(inout) :: mbin2_ice
            real(wp), intent(inout), dimension(sz2) :: yice
            real(wp), intent(inout) :: rh
        end subroutine func3
    end interface
	! ============================================================================
	! func4
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Interface-compatible fixed-grid non-collisional ice-formation callback used by bin_microphysics.
	!>@param[inout] npart: liquid-particle number concentration
	!>@param[inout] npartice: ice-particle number concentration
	!>@param[in] mwat: liquid-water mass per representative particle
	!>@param[in] mbin2: aerosol-component masses in liquid particles
	!>@param[inout] mbin2_ice: aerosol-component and ice masses in ice particles
	!>@param[in] rhobin,nubin,kappabin,molwbin: aerosol-component thermodynamic properties
	!>@param[inout] moments: conserved liquid/ice moments
	!>@param[in] medges: fixed mass-bin edges
	!>@param[inout] t: parcel temperature
	!>@param[in] p: pressure
	!>@param[in] nbins1: number of fixed bins per mode
	!>@param[in] ncomps: number of aerosol composition components
	!>@param[in] nbinw: number of liquid/warm bins
	!>@param[in] nmoms: number of conserved moments
	!>@param[in] nmodes: number of external modes
	!>@param[inout] yice: ice mass per representative particle
	!>@param[inout] rh: relative humidity
	!>@param[in] dt: timestep
	!>@param[in] sce_flag: SCE switch
	!>@param[in] mode1_flag: switch for mode-1 secondary-ice treatment
	!>@param[in] ice_nucleation_mech: logical switches for Koop, INAS, DeMott and Daily/DCMEX
    interface
        subroutine func4(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin,moments,medges, &
                         t,p,nbins1,ncomps,nbinw,nmoms,nmodes,yice,rh,dt,sce_flag, &
                         mode1_flag, ice_nucleation_mech_in) 
            use numerics_type
            import :: N_INUC_MECH
            implicit none
            real(wp), intent(inout) :: t
            real(wp), intent(in) :: p,dt
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
            logical, dimension(N_INUC_MECH), intent(in) :: ice_nucleation_mech_in
            real(wp), intent(inout) :: rh
        end subroutine func4
    end interface
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mass balance                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	parcel1%yold=parcel1%y ! store old
    if(adiabatic_prof.and.(.not.chamber_forcing_active)) then
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

	! =====================================================================
	! Numerical bin treatment of liquid condensation / evaporation.
	!
	! yold contains the pre-growth liquid-water masses.
	! y contains the accepted post-growth masses from DVODE.
	!
	! This is deliberately BEFORE chamber loss and entrainment.
	! =====================================================================
	call apply_growth_bin_scheme( parcel1%npart, &
		parcel1%yold(1:parcel1%n_bin_modew), &
		parcel1%y(1:parcel1%n_bin_modew), &
		parcel1%moments(1:parcel1%n_bin_modew,:), parcel1%mbin)
	
	if (parcel1%ice_flag.eq.1) then	
		! For liquid particles both of these mass moments are
		! the current liquid-water mass carried by the category.
		parcel1%moments(1:parcel1%n_bin_modew,parcel1%n_comps+4) = &
			parcel1%npart * parcel1%y(1:parcel1%n_bin_modew)
	
		parcel1%moments(1:parcel1%n_bin_modew,parcel1%n_comps+5) = &
			parcel1%npart * parcel1%y(1:parcel1%n_bin_modew)
	endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Routines for adjusting for entrainment          					   !                                           
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call entrainment(thresh_to_start_hom_mix)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     stk = (6._wp*parcel1%y(1:parcel1%n_bin_modew)/(pi*rhow))**(twothirds) * &
!     	235._wp/(18._wp*1.8e-5_wp*0.1)
!     stk = min(stk,1.0_wp)
!     vd = (stk*235._wp)/(1.0+stk)
!     impaction_time = 0.1_wp/vd
!     loss_rate=stk/impaction_time
!     	
! 	parcel1%npart=parcel1%npart*exp(-loss_rate*parcel1%dt)


    if(.not.adiabatic_prof) then
        ! Add aerosol carried by the newly entrained environmental-air
        ! fraction.  Existing parcel particles have already been diluted by
        ! dilute_send inside entrainment(); the environmental contribution is
        ! therefore (1-dilute_send) times the fixed t=0 environmental PSD.
        if (entrain_aerosol) then
            where ((parcel1%npart + (1._wp-dilute_send)*parcel1%npart_ent) > 0._wp)
                parcel1%y(1:parcel1%n_bin_modew)= &
                    parcel1%y(1:parcel1%n_bin_modew)*parcel1%npart + &
                    (1._wp-dilute_send)* &
                    parcel1%mbin_ent(1:parcel1%n_bin_modew,parcel1%n_comps+1) * &
                    parcel1%npart_ent
                parcel1%y(1:parcel1%n_bin_modew)= &
                    parcel1%y(1:parcel1%n_bin_modew)/ &
                    (parcel1%npart + (1._wp-dilute_send)*parcel1%npart_ent)
            end where
            parcel1%moments=parcel1%moments + &
                (1._wp-dilute_send)*parcel1%moments_ent
            parcel1%npart=parcel1%npart + &
                (1._wp-dilute_send)*parcel1%npart_ent
        endif

        ! Add aerosol residuals produced by the discrete inhomogeneous event.
        ! This is independent of whether the next timestep has already switched
        ! to homogeneous mixing.
        if (l_inhom_event) then
            where ((parcel1%npart + parcel1%npart_temp) .gt. 0._wp)
                parcel1%y(1:parcel1%n_bin_modew)= &
                    parcel1%y(1:parcel1%n_bin_modew)*parcel1%npart + &
                    parcel1%mbin_temp(1:parcel1%n_bin_modew,parcel1%n_comps+1) * &
                    parcel1%npart_temp
                parcel1%y(1:parcel1%n_bin_modew)= &
                    parcel1%y(1:parcel1%n_bin_modew)/ &
                    max(parcel1%npart + parcel1%npart_temp,1.e-30_wp)
            end where
            parcel1%moments=parcel1%moments + parcel1%moments_temp
            parcel1%npart=parcel1%npart + parcel1%npart_temp
        endif

        ! Keep the per-particle component masses consistent even when the
        ! full-moving scheme does not call moving_centre below.  Entrainment
        ! changes number and extensive component masses independently.
        do i=1,parcel1%n_bin_modew
            if (parcel1%npart(i).gt.tiny(1._wp)) then
                parcel1%mbin(i,1:parcel1%n_comps)= &
                    parcel1%moments(i,1:parcel1%n_comps)/parcel1%npart(i)
                parcel1%mbin(i,parcel1%n_comps+1)=parcel1%y(i)
            endif
        enddo

        ! Remap only when the selected numerical treatment requires it.
        if (parcel1%bin_scheme_flag.ne.BIN_FULL_MOVING) then
            call moving_centre(parcel1%n_bin_mode,parcel1%n_bin_modew,&
                    parcel1%n_bins1,parcel1%n_modes, parcel1%n_comps, &
                    parcel1%imoms+parcel1%n_comps, parcel1%npart, &
                    parcel1%y(1:parcel1%n_bin_modew), &
                    parcel1%moments(1:parcel1%n_bin_modew,:), &
                    parcel1%mbin,parcel1%mbinedges)
        endif
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
                parcel1%moments(1:2*parcel1%n_bin_modew, &
                	1:parcel1%n_comps+parcel1%imoms), &
                parcel1%mbinedges(:,:), &
                parcel1%y(parcel1%ite), &
                parcel1%y(parcel1%ipr),&
                parcel1%n_bins1, &
                n_comps,parcel1%n_bin_modew,parcel1%imoms+n_comps,parcel1%n_modes, &
                parcel1%yice(1:parcel1%n_bin_modew), &
                parcel1%y(parcel1%irh), parcel1%dt,sce_flag, &
                mode1_flag, ice_nucleation_mech) 
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
            parcel1%phi,parcel1%nump,parcel1%rhoi, parcel1%rime, &
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

		! =====================================================================
		! Numerical bin treatment of ice deposition / sublimation.
		!
		! yoldice contains the pre-growth ice mass.
		! yice contains the accepted post-growth ice mass.
		!
		! Ice shape, deposited volume and rime moments have already been
		! updated above, so they can now be transported conservatively with
		! the remapped ice population.
		! =====================================================================
		call apply_growth_bin_scheme( parcel1%npartice, &
			parcel1%yoldice(1:parcel1%n_bin_modew), &
			parcel1%yice(1:parcel1%n_bin_modew), &
			parcel1%moments(parcel1%n_bin_modew+1:parcel1%n_bin_mode,:), &
			parcel1%mbinice)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    endif



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mass balance                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(adiabatic_prof.and.(.not.chamber_forcing_active)) then
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
    
    
    
	! ============================================================================
	! inhomogeneous_liquid_reservoir
	! ============================================================================
	!>@brief
	!>Returns the liquid-water mixing ratio carried by ACTIVATED warm-bin
	!>particles before the finite entrainment dilution is applied.  Interstitial
	!>haze/aerosol water is excluded: extreme-inhomogeneous cloud mixing removes
	!>whole cloud droplets, whereas unactivated aerosol undergoes only ordinary
	!>parcel dilution.  The same activation criterion is used by ndrop.
    subroutine inhomogeneous_liquid_reservoir(ql_reservoir)
        implicit none
        real(wp), intent(out) :: ql_reservoir
        integer(i4b) :: i,n

        n=parcel1%n_bin_modew
        ql_reservoir=0._wp

        do i=1,n
            if (parcel1%npart(i).le.tiny(1._wp)) cycle
            ! Extreme-inhomogeneous evaporation acts on cloud droplets, not
            ! on the interstitial/haze aerosol population.  Use the same
            ! Koehler/FHH activation criterion as the ndrop diagnostic.
            if (.not.particle_is_activated(i,parcel1%y(i), &
                parcel1%y(parcel1%ite))) cycle
            ql_reservoir=ql_reservoir + parcel1%npart(i)*parcel1%y(i)
        enddo
    end subroutine inhomogeneous_liquid_reservoir


	! ============================================================================
	! diagnose_inhomogeneous_fraction_state
	! ============================================================================
	!>@brief
	!>Diagnoses the thermodynamic/finite-entrainment state for prescribed
	!>fractions of the DILUTED liquid and ice reservoirs removed completely.
	!>The routine has no side effects and is therefore safe to call repeatedly
	!>from the one-dimensional root solves.
	!>
	!>fliq=1 means all selected liquid remaining after ordinary dilution is
	!>evaporated.  fice=1 means all diluted ice is sublimated.  Sensible mixing,
	!>latent cooling and the finite P&K bubble/plume geometry are evaluated
	!>using prescribed D=exp(-integral(mu dz)) from the shared entrainment path.
    subroutine diagnose_inhomogeneous_fraction_state(fliq,fice,ql0,qi0, &
        dilution,tnew,qv_new,radius_new,w_new,rhw_new,rhi_new)
        implicit none
        real(wp), intent(in) :: fliq,fice,ql0,qi0,dilution
        real(wp), intent(out) :: tnew,qv_new,radius_new,w_new,rhw_new,rhi_new

        integer(i4b) :: ihaze
        real(wp) :: told,p,qv_old,svp0,rm_old,rho_old,rho_new, &
            cpm,ql_rem,qi_rem,denom,wold,r_old,dz_event,qsw,qsi, &
            ql_haze,ql_haze_old,qevap_liq,qevap_ice,rh_guess

        told=parcel1%y(parcel1%ite)
        p=parcel1%y(parcel1%ipr)
        wold=parcel1%y(parcel1%iw)
        r_old=parcel1%y(parcel1%ira)
        dz_event=max(parcel1%y(parcel1%iz)-parcel1%zlast,0._wp)

        svp0=svp_liq(told)
        qv_old=parcel1%y(parcel1%irh)*eps1*svp0/(p-svp0)
        rm_old=ra+rv*qv_old
        rho_old=p/(rm_old*told)

        ! The finite dilution factor is prescribed by the shared entrainment
        ! integral, D=exp(-integral(mu dz)).  It is therefore independent of
        ! the amount of liquid/ice subsequently removed by the chosen
        ! inhomogeneous-mixing closure.
        !
        ! Released liquid-droplet residuals retain a finite equilibrium haze-
        ! water mass.  Solve that small contribution self-consistently because
        ! the equilibrium residual water depends on the final T and RH, while
        ! those thermodynamic variables depend on the NET evaporated water.
        ql_haze=0._wp
        do ihaze=1,20
            ql_haze_old=ql_haze

            qevap_liq=max(dilution*fliq*ql0-ql_haze,0._wp)
            ! Ice residual haze water is not yet represented separately; retain
            ! the existing complete-sublimation treatment for the ice fraction.
            qevap_ice=dilution*fice*qi0

            qv_new=dilution*qv_old+(1._wp-dilution)*wvenv_send + &
                qevap_liq+qevap_ice
            qv_new=max(qv_new,0._wp)

            ql_rem=dilution*(1._wp-fliq)*ql0+ql_haze
            qi_rem=dilution*(1._wp-fice)*qi0
            cpm=max(cp+qv_new*cpv+ql_rem*cpw+qi_rem*cpi,0.5_wp*cp)

            ! Environmental sensible mixing plus latent cooling from the NET
            ! condensate mass transferred to vapour in this event.
            tnew=told+(1._wp-dilution)*(tenv_send-told) - &
                (lv*qevap_liq+ls*qevap_ice)/cpm
            if (tnew.le.100._wp) error stop &
                'Unphysical temperature in inhomogeneous mixing state'

            qsw=eps1*svp_liq(tnew)/(p-svp_liq(tnew))
            rh_guess=qv_new/max(qsw,tiny(1._wp))
            call released_liquid_haze_water(dilution*fliq,tnew,rh_guess, &
                ql_haze)

            if (abs(ql_haze-ql_haze_old).le. &
                1.e-10_wp*max(ql0,1.e-20_wp)) exit
        enddo

        rho_new=p/((ra+rv*qv_new)*tnew)

        ! Retain the existing finite P&K geometry/momentum update.  Crucially,
        ! this geometry no longer rediagnoses D; the same prescribed dilution
        ! is used by both homogeneous and inhomogeneous mixing.
        if (bubble_flag) then
            radius_new=r_old+onethird*dz_event*ent_rate - &
                onethird*r_old/rho_old*(rho_new-rho_old)
            radius_new=max(radius_new,1.e-12_wp)
            w_new=wold*(1._wp+(1._wp-radius_new**3*rho_new/ &
                max(r_old**3*rho_old,tiny(1._wp)))/gam_fac_ent2)
        else
            radius_new=r_old+0.5_wp*dz_event* &
                (ent_rate*(1._wp+gam_fac_ent2)/gam_fac_ent2) - &
                0.5_wp*r_old/rho_old*(rho_new-rho_old)
            radius_new=max(radius_new,1.e-12_wp)
            denom=gam_fac_ent2+radius_new**2*rho_new/ &
                max(r_old**2*rho_old,tiny(1._wp))
            w_new=wold*(1._wp+gam_fac_ent2)/denom
        endif

        qsw=eps1*svp_liq(tnew)/(p-svp_liq(tnew))
        rhw_new=qv_new/max(qsw,tiny(1._wp))
        if (tnew.lt.ttr .and. parcel1%ice_flag.eq.1) then
            qsi=eps1*svp_ice(tnew)/(p-svp_ice(tnew))
            rhi_new=qv_new/max(qsi,tiny(1._wp))
        else
            ! Above freezing there is no distinct ice-saturation branch.
            rhi_new=huge(1._wp)
        endif
    end subroutine diagnose_inhomogeneous_fraction_state


	! ============================================================================
	! solve_inhomogeneous_liquid_fraction
	! ============================================================================
	!>@brief
	!>Bisection solve for the fraction of the diluted warm-bin liquid reservoir
	!>that must evaporate to recover the PRE-EVENT parcel RH with respect to
	!>water.  If all
	!>available liquid is insufficient the returned fraction is one and the final
	!>RH is left below the target according to the available water/energy budget.
    subroutine solve_inhomogeneous_liquid_fraction(ql0,qi0,rh_target,fliq, &
        dilution,tnew,qv_new,radius_new,w_new,rhw_new,rhi_new)
        implicit none
        real(wp), intent(in) :: ql0,qi0,rh_target,dilution
        real(wp), intent(out) :: fliq,tnew,qv_new,radius_new,w_new, &
            rhw_new,rhi_new
        integer(i4b) :: iter
        real(wp) :: lo,hi,mid,rlo,dum_t,dum_qv,dum_r,dum_w, &
            dum_rhw,dum_rhi

        call diagnose_inhomogeneous_fraction_state(0._wp,0._wp,ql0,qi0, &
            dilution,tnew,qv_new,radius_new,w_new,rhw_new,rhi_new)
        if (ql0.le.tiny(1._wp) .or. rhw_new.ge.rh_target) then
            fliq=0._wp
            return
        endif

        call diagnose_inhomogeneous_fraction_state(1._wp,0._wp,ql0,qi0, &
            dilution,dum_t,dum_qv,dum_r,dum_w,dum_rhw,dum_rhi)
        if (dum_rhw.lt.rh_target) then
            fliq=1._wp
            tnew=dum_t; qv_new=dum_qv
            radius_new=dum_r; w_new=dum_w
            rhw_new=dum_rhw; rhi_new=dum_rhi
            return
        endif

        lo=0._wp
        hi=1._wp
        rlo=rhw_new-rh_target
        do iter=1,60
            mid=0.5_wp*(lo+hi)
            call diagnose_inhomogeneous_fraction_state(mid,0._wp,ql0,qi0, &
                dilution,dum_t,dum_qv,dum_r,dum_w,dum_rhw,dum_rhi)
            if (abs(dum_rhw-rh_target).le.1.e-8_wp) exit
            if ((dum_rhw-rh_target)*rlo.le.0._wp) then
                hi=mid
            else
                lo=mid
                rlo=dum_rhw-rh_target
            endif
        enddo
        fliq=mid
        tnew=dum_t; qv_new=dum_qv
        radius_new=dum_r; w_new=dum_w
        rhw_new=dum_rhw; rhi_new=dum_rhi
    end subroutine solve_inhomogeneous_liquid_fraction


	! ============================================================================
	! solve_inhomogeneous_ice_fraction
	! ============================================================================
	!>@brief
	!>After ALL warm-bin liquid has evaporated, solves for the fraction of diluted
	!>ice that must sublimate to reach ice saturation (RHi=1).  If all ice is
	!>insufficient, fice=1 and the final RH follows from the total budget.
    subroutine solve_inhomogeneous_ice_fraction(ql0,qi0,fice,dilution,tnew, &
        qv_new,radius_new,w_new,rhw_new,rhi_new)
        implicit none
        real(wp), intent(in) :: ql0,qi0,dilution
        real(wp), intent(out) :: fice,tnew,qv_new,radius_new,w_new, &
            rhw_new,rhi_new
        integer(i4b) :: iter
        real(wp) :: lo,hi,mid,rlo,dum_t,dum_qv,dum_r,dum_w, &
            dum_rhw,dum_rhi

        call diagnose_inhomogeneous_fraction_state(1._wp,0._wp,ql0,qi0, &
            dilution,tnew,qv_new,radius_new,w_new,rhw_new,rhi_new)
        if (qi0.le.tiny(1._wp) .or. rhi_new.ge.1._wp) then
            fice=0._wp
            return
        endif

        call diagnose_inhomogeneous_fraction_state(1._wp,1._wp,ql0,qi0, &
            dilution,dum_t,dum_qv,dum_r,dum_w,dum_rhw,dum_rhi)
        if (dum_rhi.lt.1._wp) then
            fice=1._wp
            tnew=dum_t; qv_new=dum_qv
            radius_new=dum_r; w_new=dum_w
            rhw_new=dum_rhw; rhi_new=dum_rhi
            return
        endif

        lo=0._wp
        hi=1._wp
        rlo=rhi_new-1._wp
        do iter=1,60
            mid=0.5_wp*(lo+hi)
            call diagnose_inhomogeneous_fraction_state(1._wp,mid,ql0,qi0, &
                dilution,dum_t,dum_qv,dum_r,dum_w,dum_rhw,dum_rhi)
            if (abs(dum_rhi-1._wp).le.1.e-8_wp) exit
            if ((dum_rhi-1._wp)*rlo.le.0._wp) then
                hi=mid
            else
                lo=mid
                rlo=dum_rhi-1._wp
            endif
        enddo
        fice=mid
        tnew=dum_t; qv_new=dum_qv
        radius_new=dum_r; w_new=dum_w
        rhw_new=dum_rhw; rhi_new=dum_rhi
    end subroutine solve_inhomogeneous_ice_fraction



	! ============================================================================
	! equilibrium_residual_water_mass
	! ============================================================================
	!>Return the stable (pre-activation) Koehler/FHH equilibrium water mass for
	!>a released aerosol residual at the supplied temperature and liquid-water RH.
	!>The dry component masses are supplied explicitly because a released
	!>hydrometeor can carry a different dry mass from the current warm-bin mean.
	!>Using this equilibrium haze-water mass prevents all evaporated residuals
	!>from being forced into the first wet-mass bin.
    subroutine equilibrium_residual_water_mass(ibin,mcomp,t_resid,rh_resid,mwat_eq)
        use numerics, only : zeroin, fmin
        implicit none
        integer(i4b), intent(in) :: ibin
        real(wp), dimension(:), intent(in) :: mcomp
        real(wp), intent(in) :: t_resid,rh_resid
        real(wp), intent(out) :: mwat_eq

        integer(i4b) :: n_sel_save
        real(wp) :: rh_act_save,mult_save,t_save,nwcrit,nweq,rh_use
        real(wp), dimension(size(mcomp)) :: mbin_save

        ! A residual with no nonvolatile material has no aerosol identity to
        ! preserve.  This should not normally occur for released droplets.
        if (sum(mcomp).le.tiny(1._wp)) then
            mwat_eq=0._wp
            return
        endif

        n_sel_save=n_sel
        rh_act_save=rh_act
        mult_save=mult
        t_save=parcel1%t
        mbin_save=parcel1%mbin(ibin,1:size(mcomp))

        n_sel=ibin
        parcel1%t=t_resid
        parcel1%mbin(ibin,1:size(mcomp))=mcomp

        ! Find the Koehler/FHH maximum, then solve on the stable haze branch
        ! exactly as in the initial aerosol equilibration.  Cap at 0.999 so a
        ! residual released into saturated/supersaturated air is represented as
        ! a near-saturated haze particle and can reactivate naturally afterwards.
        rh_act=0._wp
        mult=-1._wp
        select case(kappa_flag)
        case(0)
            nwcrit=fmin(1.e-50_wp,1.e1_wp,koehler02,1.e-30_wp)
            rh_use=min(max(rh_resid,1.e-8_wp),0.999_wp)
            rh_act=rh_use
            mult=1._wp
            nweq=zeroin(1.e-30_wp,nwcrit,koehler02,1.e-30_wp)
        case(1)
            nwcrit=fmin(1.e-50_wp,1.e1_wp,kkoehler02,1.e-30_wp)
            rh_use=min(max(rh_resid,1.e-8_wp),0.999_wp)
            rh_act=rh_use
            mult=1._wp
            nweq=zeroin(1.e-30_wp,nwcrit,kkoehler02,1.e-30_wp)
        case default
            parcel1%mbin(ibin,1:size(mcomp))=mbin_save
            parcel1%t=t_save
            mult=mult_save
            rh_act=rh_act_save
            n_sel=n_sel_save
            error stop 'Unknown kappa_flag in equilibrium_residual_water_mass'
        end select

        mwat_eq=max(nweq*molw_water,0._wp)

        parcel1%mbin(ibin,1:size(mcomp))=mbin_save
        parcel1%t=t_save
        mult=mult_save
        rh_act=rh_act_save
        n_sel=n_sel_save
    end subroutine equilibrium_residual_water_mass

	! ============================================================================
	! released_liquid_haze_water
	! ============================================================================
	!>@brief
	!>Returns the liquid-water mixing ratio retained on aerosol residuals from
	!>the completely evaporated ACTIVATED warm-bin population.  The supplied
	!>factor is D*fliq, so the returned haze water is on the same pre-event
	!>population basis as the inhomogeneous liquid reservoir.  If residual
	!>aerosol release is disabled, all removed droplet water is transferred to
	!>vapour and this contribution is zero.
    subroutine released_liquid_haze_water(factor,t_resid,rh_resid,ql_haze)
        implicit none
        real(wp), intent(in) :: factor,t_resid,rh_resid
        real(wp), intent(out) :: ql_haze

        integer(i4b) :: i,j,n
        real(wp) :: residual_water
        real(wp), dimension(parcel1%n_comps) :: mcomp

        ql_haze=0._wp
        if (.not.release_aerosol) return
        if (factor.le.tiny(1._wp)) return

        n=parcel1%n_bin_modew
        do i=1,n
            if (parcel1%npart(i).le.tiny(1._wp)) cycle
            if (.not.particle_is_activated(i,parcel1%y(i), &
                parcel1%y(parcel1%ite))) cycle

            do j=1,parcel1%n_comps
                mcomp(j)=parcel1%moments(i,j)/parcel1%npart(i)
            enddo

            call equilibrium_residual_water_mass(i,mcomp,t_resid,rh_resid, &
                residual_water)
            ql_haze=ql_haze+factor*parcel1%npart(i)*residual_water
        enddo
    end subroutine released_liquid_haze_water

	! ============================================================================
	! prepare_released_hydrometeor_aerosol
	! ============================================================================
	!>@brief
	!>Builds one warm-bin aerosol-residual population from the fractions of
	!>warm particles whose liquid water evaporates completely and ice particles
	!>that sublimate completely.  Nonvolatile component masses and cumulative Nin
	!>are retained;
	!>ice-only/current-ice state, including n_demott, is reset.  If
	!>release_aerosol=.false. the residual population is intentionally discarded.
    subroutine prepare_released_hydrometeor_aerosol(liq_factor,ice_factor, &
        t_resid,rh_resid,force_release)
        implicit none
        real(wp), intent(in) :: liq_factor,ice_factor,t_resid,rh_resid
        logical, intent(in), optional :: force_release
        integer(i4b) :: i,j,ii,n
        real(wp) :: residual_water
        logical :: do_release

        n=parcel1%n_bin_modew
        parcel1%npart_temp=0._wp
        parcel1%moments_temp=0._wp
        parcel1%mbin_temp=0._wp
        parcel1%mbinedges_temp=parcel1%mbinedges

        do_release=release_aerosol
        if (present(force_release)) do_release=force_release
        if (.not.do_release) return

        ! A fraction of each ACTIVATED old warm-bin population belongs to the
        ! part of the cloud whose droplets evaporate completely.  Interstitial
        ! aerosol is not part of this residual source.  Returning the activated
        ! residual population conserves aerosol number and nonvolatile mass.
        if (liq_factor.gt.tiny(1._wp)) then
            do i=1,n
                if (parcel1%npart(i).le.tiny(1._wp)) cycle
                ! Only activated cloud droplets can be completely evaporated
                ! by the extreme-inhomogeneous closure.  Interstitial aerosol
                ! is merely diluted by D and must not be removed by fliq.
                if (.not.particle_is_activated(i,parcel1%y(i), &
                    parcel1%y(parcel1%ite))) cycle
                parcel1%npart_temp(i)=parcel1%npart_temp(i)+ &
                    liq_factor*parcel1%npart(i)
                parcel1%moments_temp(i,1:parcel1%n_comps)= &
                    parcel1%moments_temp(i,1:parcel1%n_comps)+ &
                    liq_factor*parcel1%moments(i,1:parcel1%n_comps)
                if (parcel1%ice_flag.eq.1 .and. parcel1%n_inp_classes.gt.0) then
                    parcel1%moments_temp(i,parcel1%iinp_start: &
                        parcel1%iinp_start+parcel1%n_inp_classes-1)= &
                        parcel1%moments_temp(i,parcel1%iinp_start: &
                        parcel1%iinp_start+parcel1%n_inp_classes-1)+ &
                        liq_factor*parcel1%moments(i,parcel1%iinp_start: &
                        parcel1%iinp_start+parcel1%n_inp_classes-1)
                endif
            enddo
        endif

        ! Ice residuals are returned to the warm aerosol grid.  Aggregated
        ! aerosol remains coagulated because component masses are transferred
        ! ensemble-wise rather than splitting an aggregate into old monomers.
        if (parcel1%ice_flag.eq.1 .and. ice_factor.gt.tiny(1._wp)) then
            do i=1,n
                if (parcel1%npartice(i).le.tiny(1._wp)) cycle
                ii=i+n
                parcel1%npart_temp(i)=parcel1%npart_temp(i)+ &
                    ice_factor*parcel1%npartice(i)
                parcel1%moments_temp(i,1:parcel1%n_comps)= &
                    parcel1%moments_temp(i,1:parcel1%n_comps)+ &
                    ice_factor*parcel1%moments(ii,1:parcel1%n_comps)
                if (parcel1%n_inp_classes.gt.0) then
                    parcel1%moments_temp(i,parcel1%iinp_start: &
                        parcel1%iinp_start+parcel1%n_inp_classes-1)= &
                        parcel1%moments_temp(i,parcel1%iinp_start: &
                        parcel1%iinp_start+parcel1%n_inp_classes-1)+ &
                        ice_factor*parcel1%moments(ii,parcel1%iinp_start: &
                        parcel1%iinp_start+parcel1%n_inp_classes-1)
                endif
            enddo
        endif

        do i=1,n
            if (parcel1%npart_temp(i).le.tiny(1._wp)) cycle
            do j=1,parcel1%n_comps
                parcel1%mbin_temp(i,j)= &
                    parcel1%moments_temp(i,j)/parcel1%npart_temp(i)
            enddo
            call equilibrium_residual_water_mass(i, &
                parcel1%mbin_temp(i,1:parcel1%n_comps),t_resid,rh_resid, &
                residual_water)
            parcel1%mbin_temp(i,parcel1%n_comps+1)=residual_water

            if (parcel1%ice_flag.eq.1) then
                ! Warm residual convention: phi and nmon default to one per
                ! residual particle; ice volume/rime state and n_demott vanish.
                parcel1%moments_temp(i,parcel1%n_comps+1)=parcel1%npart_temp(i)
                parcel1%moments_temp(i,parcel1%n_comps+2)=parcel1%npart_temp(i)
                parcel1%moments_temp(i,parcel1%n_comps+3)=0._wp
                parcel1%moments_temp(i,parcel1%n_comps+4)= &
                    parcel1%npart_temp(i)*residual_water
                parcel1%moments_temp(i,parcel1%n_comps+5)= &
                    parcel1%npart_temp(i)*residual_water
                parcel1%moments_temp(i,parcel1%idemott)=0._wp
            endif
        enddo

        if ((sum(parcel1%npart_temp).gt.tiny(1._wp)).and. &
            (parcel1%bin_scheme_flag.ne.BIN_FULL_MOVING)) then
            call moving_centre(parcel1%n_bin_mode,parcel1%n_bin_modew, &
                parcel1%n_bins1,parcel1%n_modes,parcel1%n_comps, &
                parcel1%imoms+parcel1%n_comps,parcel1%npart_temp, &
                parcel1%mbin_temp(:,parcel1%n_comps+1), &
                parcel1%moments_temp(1:n,:),parcel1%mbin_temp, &
                parcel1%mbinedges_temp)
        endif
    end subroutine prepare_released_hydrometeor_aerosol


	! ============================================================================
	! apply_inhomogeneous_event
	! ============================================================================
	!>@brief
	!>Applies one discrete extreme-inhomogeneous event using the agreed closure:
	!>all existing particles first dilute by D.  If mixed RHw is below water
	!>saturation but at/above ice saturation, warm-bin liquid alone evaporates
	!>towards the PRE-EVENT RH.  If the mixed state is below ice saturation,
	!>liquid is still tried first; only if all warm-bin liquid is exhausted and
	!>the state remains below ice saturation is ice sublimated to RHi=1.  If a
	!>phase budget is insufficient, the resulting subsaturated RH is retained.
	!>The finite dilution D is supplied by the shared accumulated entrainment
	!>integral rather than rediagnosed from the finite geometry.
    subroutine apply_inhomogeneous_event(dilution)
        implicit none
        real(wp), intent(in) :: dilution
        integer(i4b) :: i,n
        real(wp) :: ql0,qi0,rh0,fliq,fice,tnew,qv_new, &
            radius_new,w_new,rhw_new,rhi_new,liq_res_factor,ice_res_factor, &
            factor,svp_new

        n=parcel1%n_bin_modew
        rh0=parcel1%y(parcel1%irh)
        call inhomogeneous_liquid_reservoir(ql0)
        qi0=0._wp
        if (parcel1%ice_flag.eq.1) qi0=sum(parcel1%npartice*parcel1%yice)

        ! First diagnose pure air mixing/dilution with no instantaneous phase
        ! loss.  This mixed RH decides which condensate reservoirs are allowed
        ! to participate in the extreme-inhomogeneous adjustment.
        call diagnose_inhomogeneous_fraction_state(0._wp,0._wp,ql0,qi0, &
            dilution,tnew,qv_new,radius_new,w_new,rhw_new,rhi_new)
        fliq=0._wp
        fice=0._wp

        if (rhw_new.lt.1._wp) then
            if (rhi_new.ge.1._wp) then
                ! Water-subsaturated but ice-saturated/supersaturated: liquid
                ! alone can evaporate.  Try to recover the pre-event parcel RH;
                ! if the complete warm-bin liquid budget is insufficient, keep
                ! the lower RH set by the available water and energy.
                call solve_inhomogeneous_liquid_fraction(ql0,qi0,rh0,fliq, &
                    dilution,tnew,qv_new,radius_new,w_new,rhw_new,rhi_new)
            else
                ! Below ice saturation: still use the full warm-bin liquid
                ! reservoir first.  If liquid alone can restore the pre-event
                ! RH, no ice is touched.
                call solve_inhomogeneous_liquid_fraction(ql0,qi0,rh0,fliq, &
                    dilution,tnew,qv_new,radius_new,w_new,rhw_new,rhi_new)

                if (fliq.ge.1._wp-1.e-10_wp .and. rhi_new.lt.1._wp .and. &
                    parcel1%ice_flag.eq.1 .and. qi0.gt.tiny(1._wp)) then
                    ! All warm-bin liquid has gone and the parcel remains below
                    ! ice saturation.  Sublimate ice only until RHi=1, or all
                    ! ice is exhausted if that target cannot be reached.
                    call solve_inhomogeneous_ice_fraction(ql0,qi0,fice, &
                        dilution,tnew,qv_new,radius_new,w_new,rhw_new,rhi_new)
                endif
            endif
        endif

        dilute_send=dilution
        liq_res_factor=dilution*fliq
        ice_res_factor=dilution*fice

        ! Build residual aerosol BEFORE changing the prognostic populations so
        ! component masses/Nin are taken from the complete pre-event particles.
        call prepare_released_hydrometeor_aerosol(liq_res_factor, &
            ice_res_factor,tnew,rhw_new)

		! Construct the final extreme-inhomogeneous population directly from
		! the pre-event state.  D accounts for ordinary parcel dilution and
		! (1-fliq) retains the fraction of activated droplets that survive
		! complete evaporation.  These factors form one coupled event rather
		! than two sequential population updates.
        do i=1,n
            ! All warm particles undergo ordinary parcel dilution.  The
            ! additional (1-fliq) loss applies ONLY to activated droplets.
            ! This preserves the small interstitial/haze aerosol population.
            factor=dilution
            if (parcel1%npart(i).gt.tiny(1._wp)) then
                if (particle_is_activated(i,parcel1%y(i), &
                    parcel1%y(parcel1%ite))) then
                    factor=dilution*(1._wp-fliq)
                endif
            endif
            parcel1%npart(i)=parcel1%npart(i)*factor
            parcel1%moments(i,:)=parcel1%moments(i,:)*factor
        enddo

        ! Ice always dilutes; when the mixed state remains below ice saturation
        ! after exhausting liquid, the same all-or-nothing convention is used
        ! for the diagnosed sublimated fraction fice.
        if (parcel1%ice_flag.eq.1) then
            factor=dilution*(1._wp-fice)
            parcel1%npartice=parcel1%npartice*factor
            parcel1%moments(n+1:2*n,:)=parcel1%moments(n+1:2*n,:)*factor
        endif

        parcel1%y(parcel1%ite)=tnew
        parcel1%y(parcel1%ira)=radius_new
        parcel1%y(parcel1%iw)=w_new

        svp_new=svp_liq(tnew)
        parcel1%y(parcel1%irh)=qv_new*parcel1%y(parcel1%ipr)/ &
            max((eps1+qv_new)*svp_new,tiny(1._wp))
        parcel1%y(parcel1%irh)=max(parcel1%y(parcel1%irh),0._wp)

        l_inhom_event=.true.
    end subroutine apply_inhomogeneous_event


	! ============================================================================
	! entrainment
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Coordinates environmental entrainment.  Before tth, positive
	!>entrain_period selects discrete extreme-inhomogeneous liquid-mixing events;
	!>at and after tth the model switches permanently to continuous homogeneous
	!>P&K mixing.  Environmental aerosol and evaporated-drop residual aerosol are
	!>controlled independently by entrain_aerosol and release_aerosol.
	!>@param[in] tth: time [s] at which homogeneous mixing starts
    subroutine entrainment(tth)
        use numerics_type
        use numerics, only : zeroin, fmin
        implicit none
        real(wp), intent(in) :: tth

        real(wp) :: pe,qve,te,dummy,flux_new,flux_old,rm_new,rm_old, &
            svp1,wv,rhenv,d_dummy, mu_old, mu_new, dz_ent
        integer(i4b) :: iloc,i
        logical :: inhom_step,event_due,crossing_threshold

        if (adiabatic_prof) return

        dilute_send=1._wp
        l_inhom_event=.false.
        parcel1%npart_temp=0._wp
        parcel1%moments_temp=0._wp
        parcel1%mbin_temp=0._wp

        ! Start every environmental-aerosol source from the immutable t=0
        ! distribution.  Working arrays may be hydrated/remapped below without
        ! allowing the environmental source spectrum to drift with time.
        parcel1%npart_ent=parcel1%npart_ent0
        parcel1%mbin_ent=parcel1%mbin_ent0
        parcel1%moments_ent=parcel1%moments_ent0
        parcel1%mbinedges_ent=parcel1%mbinedges_ent0

        ! Environmental state at the midpoint parcel pressure.  Parcel and
        ! environment are assumed to share pressure at a given height.
        pe=0.5_wp*(parcel1%yold(parcel1%ipr)+parcel1%y(parcel1%ipr))
        iloc=find_pos(parcel1%p_sound(1:n_levels_s),pe)
        iloc=min(n_levels_s-1,iloc)
        iloc=max(1,iloc)
        call poly_int(parcel1%p_sound(iloc:iloc+1), &
            parcel1%q_sound(1,iloc:iloc+1), &
            max(pe,parcel1%p_sound(n_levels_s)),qve,dummy)
        wvenv_send=qve
        call poly_int(parcel1%p_sound(iloc:iloc+1), &
            parcel1%t_sound(iloc:iloc+1), &
            max(pe,parcel1%p_sound(n_levels_s)),te,dummy)
        tenv_send=te
        rhenv=qve/(eps1*svp_liq(te)/(pe-svp_liq(te)))

        ! l_inhom records which pathway was active DURING the just-completed
        ! warm ODE step.  In inhomogeneous mode mu=0 throughout every ODE step;
        ! entrainment is accumulated and applied only at discrete events.
        inhom_step=l_inhom

        if (inhom_step) then
            entrain_count=entrain_count+1

            ! Accumulate exactly the same finite-dilution measure used by the
            ! homogeneous pathway.  Inhomogeneous mixing stores this integral
            ! until an event, whereas homogeneous mixing applies it each step.
            mu_old=ent_rate/max(parcel1%yold(parcel1%ira),tiny(1._wp))
            mu_new=ent_rate/max(parcel1%y(parcel1%ira),tiny(1._wp))
            dz_ent=max(parcel1%y(parcel1%iz)- &
                parcel1%yold(parcel1%iz),0._wp)
            entrain_integral=entrain_integral + &
                0.5_wp*(mu_old+mu_new)*dz_ent

            ! When the threshold is crossed, apply one final discrete event so
            ! entrainment accumulated since zlast is not lost, then use
            ! homogeneous P&K mixing from the next timestep onward.
            crossing_threshold=parcel1%tt.ge.tth
            event_due=(entrain_count.ge.max(entrain_period,1)) .or. &
                crossing_threshold

            if (event_due) then
                dilute_send=exp(-max(entrain_integral,0._wp))
                call apply_inhomogeneous_event(dilute_send)
                entrain_integral=0._wp
                entrain_count=0
                parcel1%zlast=parcel1%y(parcel1%iz)
            endif

            if (crossing_threshold) then
                l_inhom=.false.
                entrain_count=0
                entrain_integral=0._wp
            else
                l_inhom=.true.
            endif

        else
            ! Continuous homogeneous P&K thermodynamics were already integrated
            ! inside fparcelwarm.  Commit the corresponding concentration
            ! dilution to the prognostic particle populations now that the ODE
            ! step has been accepted.
            svp1=svp_liq(parcel1%yold(parcel1%ite))
            wv=eps1*parcel1%yold(parcel1%irh)*svp1/ &
                (parcel1%yold(parcel1%ipr)-svp1)
            rm_old=ra+wv*rv
            svp1=svp_liq(parcel1%y(parcel1%ite))
            wv=eps1*parcel1%y(parcel1%irh)*svp1/ &
                (parcel1%y(parcel1%ipr)-svp1)
            rm_new=ra+wv*rv

            if (bubble_flag) then
                flux_old=fourthirds*pi*parcel1%yold(parcel1%ira)**3 * &
                    parcel1%yold(parcel1%ipr)/ &
                    (parcel1%yold(parcel1%ite)*rm_old)
                flux_new=fourthirds*pi*parcel1%y(parcel1%ira)**3 * &
                    parcel1%y(parcel1%ipr)/(parcel1%y(parcel1%ite)*rm_new)
            else
                flux_old=pi*parcel1%yold(parcel1%ira)**2 * &
                    parcel1%yold(parcel1%ipr)/ &
                    (parcel1%yold(parcel1%ite)*rm_old) * &
                    parcel1%yold(parcel1%iw)
                flux_new=pi*parcel1%y(parcel1%ira)**2 * &
                    parcel1%y(parcel1%ipr)/(parcel1%y(parcel1%ite)*rm_new) * &
                    parcel1%y(parcel1%iw)
            endif

			mu_old=ent_rate/max(parcel1%yold(parcel1%ira),tiny(1._wp))
			mu_new=ent_rate/max(parcel1%y(parcel1%ira),tiny(1._wp))
			dz_ent=max(parcel1%y(parcel1%iz)- &
					   parcel1%yold(parcel1%iz),0._wp)
            entrain_integral=0.5_wp*(mu_old+mu_new)*dz_ent
			dilute_send=exp(-max(entrain_integral,0._wp))

            parcel1%npart=parcel1%npart*dilute_send
            parcel1%moments=parcel1%moments*dilute_send
            if (parcel1%ice_flag.eq.1) &
                parcel1%npartice=parcel1%npartice*dilute_send

            ! Once homogeneous mixing has begun it remains the selected pathway.
            l_inhom=.false.
            entrain_count=0
            entrain_integral=0._wp
        endif

        ! entrain_aerosol controls ONLY whether environmental aerosol accompanies
        ! the entrained environmental-air fraction.  It never changes D or the
        ! thermodynamic entrainment itself.
        if (entrain_aerosol .and. (1._wp-dilute_send).gt.tiny(1._wp)) then
            select case(kappa_flag)
            case(0)
                do i=1,parcel1%n_bin_modew
                    if (parcel1%npart_ent(i).le.0._wp) cycle
                    n_sel=i
                    rh_act=0._wp
                    mult=-1._wp
                    dummy=fmin(1.e-50_wp,1.e1_wp,koehler02_ent,1.e-30_wp)
                    mult=1._wp
                    rh_act=koehler02_ent(dummy)
                    rh_act=min(rhenv,rh_act)
                    d_dummy=zeroin(1.e-30_wp,dummy,koehler02_ent,1.e-30_wp)* &
                        molw_water
                    parcel1%mbin_ent(i,n_comps+1)=d_dummy
                enddo
            case(1)
                do i=1,parcel1%n_bin_modew
                    if (parcel1%npart_ent(i).le.0._wp) cycle
                    n_sel=i
                    rh_act=0._wp
                    mult=-1._wp
                    dummy=fmin(1.e-50_wp,1.e1_wp,kkoehler02_ent,1.e-30_wp)
                    mult=1._wp
                    rh_act=kkoehler02_ent(dummy)
                    rh_act=min(rhenv,rh_act)
                    d_dummy=zeroin(1.e-30_wp,dummy,kkoehler02_ent,1.e-30_wp)* &
                        molw_water
                    parcel1%mbin_ent(i,n_comps+1)=d_dummy
                enddo
            case default
                error stop 'Unknown kappa_flag in entrainment'
            end select

            if (parcel1%ice_flag.eq.1) then
                parcel1%moments_ent(1:parcel1%n_bin_modew,n_comps+1)= &
                    parcel1%npart_ent
                parcel1%moments_ent(1:parcel1%n_bin_modew,n_comps+2)= &
                    parcel1%npart_ent
                parcel1%moments_ent(1:parcel1%n_bin_modew,n_comps+4)= &
                    parcel1%npart_ent*parcel1%mbin_ent(:,n_comps+1)
                parcel1%moments_ent(1:parcel1%n_bin_modew,n_comps+5)= &
                    parcel1%npart_ent*parcel1%mbin_ent(:,n_comps+1)
            endif

            if (parcel1%bin_scheme_flag.ne.BIN_FULL_MOVING) then
                call moving_centre(parcel1%n_bin_mode,parcel1%n_bin_modew, &
                    parcel1%n_bins1,parcel1%n_modes,parcel1%n_comps, &
                    parcel1%imoms+parcel1%n_comps,parcel1%npart_ent, &
                    parcel1%mbin_ent(:,parcel1%n_comps+1), &
                    parcel1%moments_ent(1:parcel1%n_bin_modew,:), &
                    parcel1%mbin_ent,parcel1%mbinedges_ent)
            endif
        endif
    end subroutine entrainment
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	


    ! ============================================================================
    ! chamber_released_haze_water
    ! ============================================================================
    !>@brief
    !>Diagnoses equilibrium haze water carried by aerosol residuals released
    !>when prescribed fractions of activated liquid droplets and ice particles
    !>are completely evaporated/sublimated by chamber boundary-layer mixing.
    !>Unlike atmospheric entrainment, chamber BL recirculation always returns
    !>the nonvolatile aerosol residual to the airborne chamber population.
    subroutine chamber_released_haze_water(fliq,fice,t_resid,rh_resid, &
        qhaze_liq,qhaze_ice)
        implicit none
        real(wp), intent(in) :: fliq,fice,t_resid,rh_resid
        real(wp), intent(out) :: qhaze_liq,qhaze_ice
        integer(i4b) :: i,j,n,ii
        real(wp) :: residual_water
        real(wp), dimension(parcel1%n_comps) :: mcomp

        n=parcel1%n_bin_modew
        qhaze_liq=0._wp
        qhaze_ice=0._wp

        if (fliq.gt.tiny(1._wp)) then
            do i=1,n
                if (parcel1%npart(i).le.tiny(1._wp)) cycle
                if (.not.particle_is_activated(i,parcel1%y(i), &
                    parcel1%y(parcel1%ite))) cycle
                do j=1,parcel1%n_comps
                    mcomp(j)=parcel1%moments(i,j)/parcel1%npart(i)
                enddo
                call equilibrium_residual_water_mass(i,mcomp,t_resid,rh_resid, &
                    residual_water)
                qhaze_liq=qhaze_liq+fliq*parcel1%npart(i)*residual_water
            enddo
        endif

        if (parcel1%ice_flag.eq.1 .and. fice.gt.tiny(1._wp)) then
            do i=1,n
                if (parcel1%npartice(i).le.tiny(1._wp)) cycle
                ii=n+i
                do j=1,parcel1%n_comps
                    mcomp(j)=parcel1%moments(ii,j)/parcel1%npartice(i)
                enddo
                call equilibrium_residual_water_mass(i,mcomp,t_resid,rh_resid, &
                    residual_water)
                qhaze_ice=qhaze_ice+fice*parcel1%npartice(i)*residual_water
            enddo
        endif
    end subroutine chamber_released_haze_water


    ! ============================================================================
    ! diagnose_chamber_bl_inhomogeneous_state
    ! ============================================================================
    !>@brief
    !>Side-effect-free thermodynamic diagnosis for extreme-inhomogeneous chamber
    !>BL mixing.  fliq/fice are fractions of the CURRENT activated-liquid and
    !>ice populations that completely disappear as hydrometeors.  External BL
    !>water exchange has already produced qv_air/t_air.  Hydrometeor evaporation
    !>only redistributes water internally, so it cannot change total water.
    subroutine diagnose_chamber_bl_inhomogeneous_state(fliq,fice,ql_act,ql_total, &
        qi0,qv_air,t_air,p,tnew,qv_new,rhw_new,rhi_new,qhaze_liq,qhaze_ice)
        implicit none
        real(wp), intent(in) :: fliq,fice,ql_act,ql_total,qi0,qv_air,t_air,p
        real(wp), intent(out) :: tnew,qv_new,rhw_new,rhi_new,qhaze_liq,qhaze_ice
        integer(i4b) :: iter
        real(wp) :: qh_l_old,qh_i_old,qevap_liq,qevap_ice,ql_rem,qi_rem,cpm, &
            qsw,qsi,rh_guess

        qhaze_liq=0._wp
        qhaze_ice=0._wp
        tnew=t_air

        do iter=1,30
            qh_l_old=qhaze_liq
            qh_i_old=qhaze_ice

            qevap_liq=max(fliq*ql_act-qhaze_liq,0._wp)
            qevap_ice=max(fice*qi0-qhaze_ice,0._wp)
            qv_new=max(qv_air+qevap_liq+qevap_ice,0._wp)

            ql_rem=max(ql_total-fliq*ql_act,0._wp)+qhaze_liq+qhaze_ice
            qi_rem=max((1._wp-fice)*qi0,0._wp)
            cpm=max(cp+qv_new*cpv+ql_rem*cpw+qi_rem*cpi,0.5_wp*cp)

            if (chamber_force_temperature) then
                tnew=t_air
            else
                ! Liquid residual water never vaporised; ice residual water is
                ! treated as sublimation followed by condensation onto aerosol.
                tnew=t_air-(lv*qevap_liq+ls*fice*qi0-lv*qhaze_ice)/cpm
            endif
            if (tnew.le.100._wp) error stop &
                'Unphysical temperature in chamber BL inhomogeneous mixing'

            qsw=eps1*svp_liq(tnew)/max(p-svp_liq(tnew),tiny(1._wp))
            rh_guess=qv_new/max(qsw,tiny(1._wp))
            call chamber_released_haze_water(fliq,fice,tnew,rh_guess, &
                qhaze_liq,qhaze_ice)
            qhaze_liq=min(max(qhaze_liq,0._wp),max(fliq*ql_act,0._wp))
            qhaze_ice=min(max(qhaze_ice,0._wp),max(fice*qi0,0._wp))

            if (abs(qhaze_liq-qh_l_old)+abs(qhaze_ice-qh_i_old).le. &
                1.e-10_wp*max(ql_total+qi0,1.e-20_wp)) exit
        enddo

        ! Re-evaluate with the converged residual-water estimate.
        qevap_liq=max(fliq*ql_act-qhaze_liq,0._wp)
        qevap_ice=max(fice*qi0-qhaze_ice,0._wp)
        qv_new=max(qv_air+qevap_liq+qevap_ice,0._wp)
        ql_rem=max(ql_total-fliq*ql_act,0._wp)+qhaze_liq+qhaze_ice
        qi_rem=max((1._wp-fice)*qi0,0._wp)
        cpm=max(cp+qv_new*cpv+ql_rem*cpw+qi_rem*cpi,0.5_wp*cp)
        if (chamber_force_temperature) then
            tnew=t_air
        else
            tnew=t_air-(lv*qevap_liq+ls*fice*qi0-lv*qhaze_ice)/cpm
        endif

        qsw=eps1*svp_liq(tnew)/max(p-svp_liq(tnew),tiny(1._wp))
        rhw_new=qv_new/max(qsw,tiny(1._wp))
        if (tnew.lt.ttr .and. parcel1%ice_flag.eq.1) then
            qsi=eps1*svp_ice(tnew)/max(p-svp_ice(tnew),tiny(1._wp))
            rhi_new=qv_new/max(qsi,tiny(1._wp))
        else
            rhi_new=huge(1._wp)
        endif
    end subroutine diagnose_chamber_bl_inhomogeneous_state


    ! ============================================================================
    ! solve_chamber_bl_liquid_fraction
    ! ============================================================================
    subroutine solve_chamber_bl_liquid_fraction(ql_act,ql_total,qi0,qv_air,t_air,p, &
        rh_target,fmax,fliq,tnew,qv_new,rhw_new,rhi_new)
        implicit none
        real(wp), intent(in) :: ql_act,ql_total,qi0,qv_air,t_air,p,rh_target,fmax
        real(wp), intent(out) :: fliq,tnew,qv_new,rhw_new,rhi_new
        integer(i4b) :: iter
        real(wp) :: lo,hi,mid,rlo,dum_t,dum_qv,dum_rhw,dum_rhi,dum_qlh,dum_qih, &
            qlh,qih

        call diagnose_chamber_bl_inhomogeneous_state(0._wp,0._wp,ql_act,ql_total, &
            qi0,qv_air,t_air,p,tnew,qv_new,rhw_new,rhi_new,qlh,qih)
        if (ql_act.le.tiny(1._wp) .or. fmax.le.tiny(1._wp) .or. &
            rhw_new.ge.rh_target) then
            fliq=0._wp
            return
        endif

        call diagnose_chamber_bl_inhomogeneous_state(fmax,0._wp,ql_act,ql_total, &
            qi0,qv_air,t_air,p,dum_t,dum_qv,dum_rhw,dum_rhi,dum_qlh,dum_qih)
        if (dum_rhw.lt.rh_target) then
            fliq=fmax
            tnew=dum_t; qv_new=dum_qv; rhw_new=dum_rhw; rhi_new=dum_rhi
            return
        endif

        lo=0._wp
        hi=fmax
        rlo=rhw_new-rh_target
        mid=0._wp
        do iter=1,60
            mid=0.5_wp*(lo+hi)
            call diagnose_chamber_bl_inhomogeneous_state(mid,0._wp,ql_act,ql_total, &
                qi0,qv_air,t_air,p,dum_t,dum_qv,dum_rhw,dum_rhi,dum_qlh,dum_qih)
            if (abs(dum_rhw-rh_target).le.1.e-8_wp) exit
            if ((dum_rhw-rh_target)*rlo.le.0._wp) then
                hi=mid
            else
                lo=mid
                rlo=dum_rhw-rh_target
            endif
        enddo
        fliq=mid
        tnew=dum_t; qv_new=dum_qv; rhw_new=dum_rhw; rhi_new=dum_rhi
    end subroutine solve_chamber_bl_liquid_fraction


    ! ============================================================================
    ! solve_chamber_bl_ice_fraction
    ! ============================================================================
    subroutine solve_chamber_bl_ice_fraction(fliq,ql_act,ql_total,qi0,qv_air,t_air,p, &
        fmax,fice,tnew,qv_new,rhw_new,rhi_new)
        implicit none
        real(wp), intent(in) :: fliq,ql_act,ql_total,qi0,qv_air,t_air,p,fmax
        real(wp), intent(out) :: fice,tnew,qv_new,rhw_new,rhi_new
        integer(i4b) :: iter
        real(wp) :: lo,hi,mid,rlo,dum_t,dum_qv,dum_rhw,dum_rhi,dum_qlh,dum_qih, &
            qlh,qih

        call diagnose_chamber_bl_inhomogeneous_state(fliq,0._wp,ql_act,ql_total, &
            qi0,qv_air,t_air,p,tnew,qv_new,rhw_new,rhi_new,qlh,qih)
        if (qi0.le.tiny(1._wp) .or. fmax.le.tiny(1._wp) .or. rhi_new.ge.1._wp) then
            fice=0._wp
            return
        endif

        call diagnose_chamber_bl_inhomogeneous_state(fliq,fmax,ql_act,ql_total, &
            qi0,qv_air,t_air,p,dum_t,dum_qv,dum_rhw,dum_rhi,dum_qlh,dum_qih)
        if (dum_rhi.lt.1._wp) then
            fice=fmax
            tnew=dum_t; qv_new=dum_qv; rhw_new=dum_rhw; rhi_new=dum_rhi
            return
        endif

        lo=0._wp
        hi=fmax
        rlo=rhi_new-1._wp
        mid=0._wp
        do iter=1,60
            mid=0.5_wp*(lo+hi)
            call diagnose_chamber_bl_inhomogeneous_state(fliq,mid,ql_act,ql_total, &
                qi0,qv_air,t_air,p,dum_t,dum_qv,dum_rhw,dum_rhi,dum_qlh,dum_qih)
            if (abs(dum_rhi-1._wp).le.1.e-8_wp) exit
            if ((dum_rhi-1._wp)*rlo.le.0._wp) then
                hi=mid
            else
                lo=mid
                rlo=dum_rhi-1._wp
            endif
        enddo
        fice=mid
        tnew=dum_t; qv_new=dum_qv; rhw_new=dum_rhw; rhi_new=dum_rhi
    end subroutine solve_chamber_bl_ice_fraction


    ! ============================================================================
    ! merge_released_aerosol_into_warm
    ! ============================================================================
    !>@brief
    !>Merges the temporary residual-aerosol population into the prognostic warm
    !>population while conserving number, nonvolatile component moments, INP
    !>moments and residual haze water.
    subroutine merge_released_aerosol_into_warm()
        implicit none
        integer(i4b) :: i,j,n
        real(wp) :: nnew

        n=parcel1%n_bin_modew
        if (sum(parcel1%npart_temp).le.tiny(1._wp)) return

        do i=1,n
            nnew=parcel1%npart(i)+parcel1%npart_temp(i)
            if (nnew.gt.tiny(1._wp)) then
                parcel1%y(i)=(parcel1%y(i)*parcel1%npart(i)+ &
                    parcel1%mbin_temp(i,parcel1%n_comps+1)*parcel1%npart_temp(i))/nnew
            endif
        enddo
        parcel1%moments(1:n,:)=parcel1%moments(1:n,:)+parcel1%moments_temp(1:n,:)
        parcel1%npart=parcel1%npart+parcel1%npart_temp

        do i=1,n
            if (parcel1%npart(i).le.tiny(1._wp)) cycle
            do j=1,parcel1%n_comps
                parcel1%mbin(i,j)=parcel1%moments(i,j)/parcel1%npart(i)
            enddo
            parcel1%mbin(i,parcel1%n_comps+1)=parcel1%y(i)
        enddo

        if (parcel1%bin_scheme_flag.ne.BIN_FULL_MOVING) then
            call moving_centre(parcel1%n_bin_mode,parcel1%n_bin_modew, &
                parcel1%n_bins1,parcel1%n_modes,parcel1%n_comps, &
                parcel1%imoms+parcel1%n_comps,parcel1%npart, &
                parcel1%y(1:n),parcel1%moments(1:n,:),parcel1%mbin, &
                parcel1%mbinedges)
        endif

        parcel1%npart_temp=0._wp
        parcel1%moments_temp=0._wp
        parcel1%mbin_temp=0._wp
    end subroutine merge_released_aerosol_into_warm


    ! ============================================================================
    ! apply_chamber_bl_exchange
    ! ============================================================================
    !>@brief
    !>Applies one chamber boundary-layer recirculation operator.  Modes 1 and 2
    !>use the SAME external total-water exchange
    !>
    !>  dq_BL = (qv_BL-qv) [1-exp(-dt/tau_BL)] .
    !>
    !>Mode 1 leaves the resulting humidity deficit to the ordinary diffusional
    !>microphysics (homogeneous response).  Mode 2 instead permits at most the
    !>same recirculated fraction of activated droplets/ice to disappear
    !>completely (extreme inhomogeneous response).  Internal phase changes are
    !>water-conserving, and a final budget closure enforces identical total-water
    !>change in the two modes to roundoff.
    subroutine apply_chamber_bl_exchange()
        use numerics, only : find_pos, poly_int
        implicit none
        integer(i4b) :: i,n,iloc
        real(wp) :: p,t0,rh0,svp0,qv0,qv_bl,t_bl,qvs_bl,fmix,dq_wall,qv_air,t_air, &
            ql_total,ql_act,qi0,qtot_target,rh_target,fliq,fice,tnew,qv_new, &
            rhw_new,rhi_new,qh_l,qh_i,qcond_new,qv_budget,svp_new,factor,var,dummy, &
            tquery

        parcel1%qchamber_bl_step=0._wp
        if (chamber_bl_mix.eq.0) return
        if (chamber_bl_tau.le.0._wp) error stop 'chamber_bl_tau must be > 0'

        n=parcel1%n_bin_modew
        p=parcel1%y(parcel1%ipr)
        t0=parcel1%y(parcel1%ite)
        rh0=parcel1%y(parcel1%irh)
        svp0=svp_liq(t0)
        qv0=eps1*rh0*svp0/max(p-svp0,tiny(1._wp))

        ql_total=sum(parcel1%npart*parcel1%y(1:n))
        call inhomogeneous_liquid_reservoir(ql_act)
        qi0=0._wp
        if (parcel1%ice_flag.eq.1) qi0=sum(parcel1%npartice*parcel1%yice(1:n))

        ! Fraction of chamber air that completes a BL excursion during dt.
        fmix=1._wp-exp(-parcel1%dt/chamber_bl_tau)
        fmix=max(0._wp,min(fmix,1._wp))

        ! Air enters the BL carrying the bulk vapour mixing ratio.  The BL
        ! temperature is either a fixed offset from the evolving gas temperature
        ! or the measured wall-temperature time series.  There is deliberately
        ! no prescribed BL RH reservoir.
        select case(chamber_bl_temp_mode)
        case(0)
            t_bl=t0+chamber_bl_temp_offset
        case(1)
            if (n_levels_c.lt.2) error stop 'Wall-temperature BL mode requires chamber data'
            tquery=min(max(parcel1%tt,time_chamber(1)),time_chamber(n_levels_c))
            iloc=find_pos(time_chamber(1:n_levels_c),tquery)
            iloc=max(1,min(iloc,n_levels_c-1))
            call poly_int(time_chamber(iloc:iloc+1),wall_temp_chamber(iloc:iloc+1), &
                tquery,var,dummy)
            t_bl=var
        case default
            error stop 'Invalid chamber_bl_temp_mode'
        end select

        ! Conserve qv during the BL temperature excursion unless the air would
        ! exceed liquid saturation at T_bl.  In that case excess vapour is
        ! assumed to condense on the wall and is removed from airborne total
        ! water.  The wall is not an infinite source of vapour, so dq_wall <= 0.
        qvs_bl=eps1*svp_liq(t_bl)/max(p-svp_liq(t_bl),tiny(1._wp))
        qv_bl=min(qv0,qvs_bl)
        dq_wall=fmix*(qv_bl-qv0)
        dq_wall=max(dq_wall,-qv0)
        qv_air=max(qv0+dq_wall,0._wp)

        ! If gas temperature is observationally forced, do not add a second
        ! bulk sensible-heat tendency.  Otherwise mix the BL temperature back
        ! into the well-mixed chamber in proportion to fmix.
        if (chamber_force_temperature) then
            t_air=t0
        else
            t_air=t0+fmix*(t_bl-t0)
        endif

        qtot_target=qv0+ql_total+qi0+dq_wall
        parcel1%qchamber_bl_step=-dq_wall
        parcel1%qchamber_bl=parcel1%qchamber_bl-dq_wall

        fliq=0._wp
        fice=0._wp
        tnew=t_air
        qv_new=qv_air

        if (chamber_bl_mix.eq.2) then
            ! Extreme-inhomogeneous response: no more than the recirculated
            ! fraction fmix of hydrometeors can completely evaporate/sublimate.
            ! Residual aerosol is always returned to the airborne population.
            rh_target=min(max(rh0,0._wp),1._wp)
            call diagnose_chamber_bl_inhomogeneous_state(0._wp,0._wp, &
                ql_act,ql_total,qi0,qv_air,t_air,p,tnew,qv_new,rhw_new,rhi_new, &
                qh_l,qh_i)

            if (rhw_new.lt.rh_target) then
                call solve_chamber_bl_liquid_fraction(ql_act,ql_total,qi0, &
                    qv_air,t_air,p,rh_target,fmix,fliq,tnew,qv_new,rhw_new,rhi_new)

                if (parcel1%ice_flag.eq.1 .and. tnew.lt.ttr .and. &
                    qi0.gt.tiny(1._wp) .and. rhi_new.lt.1._wp .and. &
                    (ql_act.le.tiny(1._wp) .or. fliq.ge.fmix-1.e-10_wp)) then
                    call solve_chamber_bl_ice_fraction(fliq,ql_act,ql_total,qi0, &
                        qv_air,t_air,p,fmix,fice,tnew,qv_new,rhw_new,rhi_new)
                endif
            endif

            call prepare_released_hydrometeor_aerosol(fliq,fice,tnew,rhw_new,.true.)

            do i=1,n
                factor=1._wp
                if (parcel1%npart(i).gt.tiny(1._wp)) then
                    if (particle_is_activated(i,parcel1%y(i),t0)) factor=1._wp-fliq
                endif
                parcel1%npart(i)=parcel1%npart(i)*factor
                parcel1%moments(i,:)=parcel1%moments(i,:)*factor
            enddo
            if (parcel1%ice_flag.eq.1) then
                factor=1._wp-fice
                parcel1%npartice=parcel1%npartice*factor
                parcel1%moments(n+1:2*n,:)=parcel1%moments(n+1:2*n,:)*factor
            endif

            call merge_released_aerosol_into_warm()
        elseif (chamber_bl_mix.ne.1) then
            error stop 'Invalid chamber_bl_mix in apply_chamber_bl_exchange'
        endif

        ! Internal evaporation/condensation only redistributes airborne water.
        ! Close the budget exactly to the external wall sink dq_wall, ensuring
        ! homogeneous and inhomogeneous modes have the same total-water change.
        qcond_new=sum(parcel1%npart*parcel1%y(1:n))
        if (parcel1%ice_flag.eq.1) qcond_new=qcond_new+ &
            sum(parcel1%npartice*parcel1%yice(1:n))
        qv_budget=qtot_target-qcond_new
        if (qv_budget.lt.-1.e-12_wp) error stop &
            'Chamber BL phase adjustment exceeded available total water'
        qv_new=max(qv_budget,0._wp)

        parcel1%y(parcel1%ite)=tnew
        svp_new=svp_liq(tnew)
        parcel1%y(parcel1%irh)=qv_new*max(p-svp_new,tiny(1._wp))/ &
            max(eps1*svp_new,tiny(1._wp))
        parcel1%y(parcel1%irh)=max(parcel1%y(parcel1%irh),0._wp)
        parcel1%qtot=qtot_target

        if (parcel1%ice_flag.eq.1) then
            parcel1%yice(parcel1%itei)=parcel1%y(parcel1%ite)
            parcel1%yice(parcel1%irhi)=parcel1%y(parcel1%irh)
        endif
    end subroutine apply_chamber_bl_exchange


    ! ============================================================================
    ! chamber_fan_loss_rate
    ! ============================================================================
    !>@brief
    !>Sigmoid fan-blade collection rate.  The half-maximum diameter is defined
    !>at chamber_fan_rpm_ref and shifts as RPM**(-1/2), consistent with
    !>Stokes-number scaling Stk ~ D^2 U for a fixed fan geometry.
    real(wp) function chamber_fan_loss_rate(d) result(rate)
        implicit none
        real(wp), intent(in) :: d
        real(wp) :: d50,ratio

        rate=0._wp
        if (chamber_fan_loss.eq.0 .or. d.le.0._wp) return
        d50=chamber_fan_loss_d50_ref*sqrt(chamber_fan_rpm_ref/chamber_fan_rpm)
        ratio=(d50/max(d,tiny(1._wp)))**chamber_fan_loss_exp
        rate=chamber_fan_loss_kmax/(1._wp+ratio)
    end function chamber_fan_loss_rate


    ! ============================================================================
    ! chamber_wall_loss_rate
    ! ============================================================================
    !>@brief
    !>Non-gravitational deposition loss to cylindrical chamber side walls and
    !>ceiling using the smooth-surface three-layer model of Lai & Nazaroff
    !>(2000).  The floor is deliberately excluded because gravitational floor
    !>loss is already represented by apply_particle_fallout.  The model uses
    !>Brownian and turbulent diffusion and the effect of settling away from the
    !>downward-facing ceiling.  For very large r+ outside the intended smooth-
    !>wall aerosol range, r+ is capped at the edge of the viscous-layer fit.
    real(wp) function chamber_wall_loss_rate(d,vsettle,t,p) result(rate)
        implicit none
        real(wp), intent(in) :: d,vsettle,t,p
        real(wp), parameter :: kb=1.380649e-23_wp
        real(wp) :: mu,rhoa,nu,lambda_air,cc,db,sc,scm13,rplus,a,b,integ, &
                    vdvert,vdceil,xarg,radius,volume,aside,aceil,sqrt3

        rate=0._wp
        if (chamber_wall_loss.eq.0 .or. d.le.0._wp) return
        if (chamber_diameter.le.0._wp .or. chamber_height.le.0._wp) return

        mu=viscosity_air(t)
        rhoa=p/(ra*t)
        nu=mu/rhoa
        lambda_air=6.6e-8_wp*(101325._wp/p)*(t/293.15_wp)
        cc=1._wp+2._wp*lambda_air/d*(1.257_wp+0.4_wp* &
            exp(-1.1_wp*d/(2._wp*lambda_air)))
        db=kb*t*cc/(3._wp*pi*mu*d)
        if (db.le.0._wp .or. nu.le.0._wp) return

        sc=max(nu/db,1._wp)
        scm13=sc**(-1._wp/3._wp)
        rplus=0.5_wp*d*chamber_wall_ustar/nu
        rplus=max(0._wp,min(rplus,4.299_wp))
        sqrt3=sqrt(3._wp)

        a=0.5_wp*log((10.92_wp*scm13+4.3_wp)**3 / &
            (1._wp/sc+0.0609_wp)) + sqrt3*atan((8.6_wp-10.92_wp*scm13)/ &
            (sqrt3*10.92_wp*scm13))
        b=0.5_wp*log((10.92_wp*scm13+rplus)**3 / &
            (1._wp/sc+7.669e-4_wp*rplus**3)) + sqrt3*atan((2._wp*rplus- &
            10.92_wp*scm13)/(sqrt3*10.92_wp*scm13))
        integ=max(3.64_wp*sc**(2._wp/3._wp)*(a-b)+39._wp,tiny(1._wp))

        vdvert=chamber_wall_ustar/integ
        if (vsettle.le.tiny(1._wp)) then
            vdceil=vdvert
        else
            xarg=max(vsettle,0._wp)*integ/chamber_wall_ustar
            if (xarg.lt.1.e-6_wp) then
                vdceil=vdvert
            elseif (xarg.gt.50._wp) then
                vdceil=0._wp
            else
                vdceil=max(vsettle,0._wp)/(exp(xarg)-1._wp)
            endif
        endif

        radius=0.5_wp*chamber_diameter
        volume=pi*radius**2*chamber_height
        aside=2._wp*pi*radius*chamber_height
        aceil=pi*radius**2
        rate=(aside*vdvert+aceil*vdceil)/volume
        rate=max(rate,0._wp)
    end function chamber_wall_loss_rate


    ! ============================================================================
    ! apply_chamber_particle_losses
    ! ============================================================================
    !>@brief
    !>Applies independent fan-blade impaction and non-gravitational wall
    !>deposition to every airborne warm/aerosol and ice category.  Fan loss uses
    !>a saturating sigmoid in current particle diameter; wall loss uses the
    !>Lai-Nazaroff smooth-wall deposition velocity with chamber geometry and
    !>friction velocity.  Gravitational floor settling remains in
    !>apply_particle_fallout.  All extensive moments use the same exact
    !>survival fraction as number.
    subroutine apply_chamber_particle_losses()
        implicit none
        integer(i4b) :: i,n
        real(wp) :: kfan,kwall,ktot,survival,nold,nremoved,qremoved, &
            fanfrac,wallfrac,t,p

        if (chamber_fan_loss.eq.0 .and. chamber_wall_loss.eq.0) return
        n=parcel1%n_bin_modew
        t=parcel1%y(parcel1%ite)
        p=parcel1%y(parcel1%ipr)

        do i=1,n
            if (parcel1%npart(i).le.tiny(1._wp)) cycle
            kfan=chamber_fan_loss_rate(parcel1%dw(i))
            kwall=chamber_wall_loss_rate(parcel1%dw(i),parcel1%vel(i),t,p)
            ktot=kfan+kwall
            if (ktot.le.0._wp) cycle

            survival=exp(-ktot*parcel1%dt)
            nold=parcel1%npart(i)
            nremoved=nold*(1._wp-survival)
            qremoved=parcel1%y(i)*nremoved
            fanfrac=kfan/ktot
            wallfrac=kwall/ktot

            parcel1%nfan_liq=parcel1%nfan_liq+fanfrac*nremoved
            parcel1%qfan_liq=parcel1%qfan_liq+fanfrac*qremoved
            parcel1%nwall_liq=parcel1%nwall_liq+wallfrac*nremoved
            parcel1%qwall_liq=parcel1%qwall_liq+wallfrac*qremoved

            parcel1%npart(i)=nold*survival
            parcel1%moments(i,:)=parcel1%moments(i,:)*survival
        enddo

        if (parcel1%ice_flag.eq.1) then
            do i=1,n
                if (parcel1%npartice(i).le.tiny(1._wp)) cycle
                kfan=chamber_fan_loss_rate(parcel1%dwice(i))
                kwall=chamber_wall_loss_rate(parcel1%dwice(i),parcel1%vel(n+i),t,p)
                ktot=kfan+kwall
                if (ktot.le.0._wp) cycle

                survival=exp(-ktot*parcel1%dt)
                nold=parcel1%npartice(i)
                nremoved=nold*(1._wp-survival)
                qremoved=parcel1%yice(i)*nremoved
                fanfrac=kfan/ktot
                wallfrac=kwall/ktot

                parcel1%nfan_ice=parcel1%nfan_ice+fanfrac*nremoved
                parcel1%qfan_ice=parcel1%qfan_ice+fanfrac*qremoved
                parcel1%nwall_ice=parcel1%nwall_ice+wallfrac*nremoved
                parcel1%qwall_ice=parcel1%qwall_ice+wallfrac*qremoved

                parcel1%npartice(i)=nold*survival
                parcel1%moments(n+i,:)=parcel1%moments(n+i,:)*survival
            enddo
        endif
    end subroutine apply_chamber_particle_losses


	! ============================================================================
	! update_terminal_velocities
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Recomputes current wet diameters, liquid terminal velocities and, when ice is enabled, ice
	!>properties, projected areas and terminal velocities.
	subroutine update_terminal_velocities()
		implicit none
		integer(i4b) :: n
	
		n=parcel1%n_bin_modew
		! Current wet diameter and particle density
		call wetdiam( &
			parcel1%y(1:n), &
			parcel1%mbin(:,1:parcel1%n_comps), &
			parcel1%rhobin(:,1:parcel1%n_comps), &
			n, parcel1%dw, parcel1%rhoat)
		! liquid
		call terminal01( &
			parcel1%vel(1:n), &
			parcel1%dw, &
			parcel1%rhoat, &
			parcel1%y(parcel1%ite), &
			parcel1%y(parcel1%ipr), &
			parcel1%nre(1:n), &
			parcel1%cd(1:n), n)
		! ice
		if (parcel1%ice_flag.eq.1) then

			call ice_particle_properties_from_moments( &
				parcel1%yice(1:n), parcel1%npartice, parcel1%moments, &
				n, parcel1%n_comps, &
				parcel1%phi, parcel1%nump, parcel1%rhoi, parcel1%rime)		
		
			call terminal02( &
				parcel1%vel(n+1:2*n), &
				parcel1%yice(1:n), &
				parcel1%y(parcel1%ite), &
				parcel1%y(parcel1%ipr), &
				parcel1%phi, &
				parcel1%rhoi, &
				parcel1%nump, &
				parcel1%rime, &
				parcel1%nre(n+1:2*n), n, parcel1%dwice,parcel1%areaice)
		endif
	end subroutine update_terminal_velocities
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! ============================================================================
	! apply_particle_fallout
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Applies finite-residence-depth fallout to liquid and ice populations using exponential survival
	!>and removes all extensive moments by the same fraction.
	!> Particle fallout from finite-depth parcel
	!>
	!> Treat parcel as a well-mixed volume of vertical depth residence_depth.
	!> Particles are removed with residence time:
	!>
	!>      tau = residence_depth / terminal_velocity
	!>
	!> Hence:
	!>
	!>      N(t+dt) = N(t) exp(-Vt dt / residence_depth)
	!>
	!> Per-particle masses and properties are unchanged.  All extensive moments
	!> must be reduced by the same factor as particle number.
	subroutine apply_particle_fallout()
		implicit none
		integer(i4b) :: i,n
		real(wp) :: survival, qfall,nfall,fall_depth
		if (.not.fallout_flag) return
		
        ! For chamber runs with specified geometry, use the physical chamber
        ! height as V/A_floor.  Generic parcel/LES runs retain residence_depth.
        if (n_levels_c.gt.0 .and. chamber_height.gt.0._wp) then
            if (chamber_height.le.0._wp) error stop 'chamber_height must be > 0'
        elseif (residence_depth.le.0._wp) then
            error stop 'residence_depth must be > 0'
        endif
		n=parcel1%n_bin_modew
        fall_depth=residence_depth
        if (n_levels_c.gt.0 .and. chamber_height.gt.0._wp) fall_depth=chamber_height
		parcel1%qfall_step_liq=0._wp
		parcel1%qfall_step_ice=0._wp
		!------------------------------------------------------------------
		! Liquid particles
		!------------------------------------------------------------------
		do i=1,n
			survival=exp( &
				-max(parcel1%vel(i),0._wp) * &
				 parcel1%dt/fall_depth)	
			! Number removed during this timestep
			nfall=parcel1%npart(i)*(1._wp-survival)
			! Liquid-water mass removed
			qfall=parcel1%y(i)*nfall
			! Timestep diagnostics
			parcel1%qfall_step_liq = &
				parcel1%qfall_step_liq + qfall
			! Cumulative diagnostics
			parcel1%qfall_liq = &
				parcel1%qfall_liq + qfall
			parcel1%nfall_liq = &
				parcel1%nfall_liq + nfall
			! Remaining particles
			parcel1%npart(i)= &
				parcel1%npart(i)*survival
			! All extensive moments lose the same fraction
			parcel1%moments(i,:)= &
				parcel1%moments(i,:)*survival
		enddo
		!------------------------------------------------------------------
		! Ice particles
		!------------------------------------------------------------------
		if (parcel1%ice_flag.eq.1) then
			do i=1,n
				survival=exp( &
					-max(parcel1%vel(n+i),0._wp) * &
					 parcel1%dt/fall_depth)
				nfall=parcel1%npartice(i)*(1._wp-survival)
				qfall=parcel1%yice(i)*nfall
				parcel1%qfall_step_ice = &
					parcel1%qfall_step_ice + qfall
				parcel1%qfall_ice = &
					parcel1%qfall_ice + qfall
				parcel1%nfall_ice = &
					parcel1%nfall_ice + nfall
				parcel1%npartice(i)= &
					parcel1%npartice(i)*survival
				parcel1%moments(n+i,:)= &
					parcel1%moments(n+i,:)*survival
			enddo
		endif
	end subroutine apply_particle_fallout
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
	! ============================================================================
	! check
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Checks a NetCDF status code and terminates with the corresponding NetCDF error message on failure.
	!>@param[in] status: NetCDF status code returned by an nf90 call
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
    
    
	! ============================================================================
	! output
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Writes the current BMM parcel, particle, moment, optical and fallout diagnostics to the NetCDF
	!>output file.
	!>output 1 time-step of model
	!>@param[inout] new_file: true when the NetCDF output file and variables must be created; reset
	!>after initial creation
	!>@param[in] outputfile: path/name of the NetCDF output file
    subroutine output(new_file,outputfile)

    use numerics_type
    use netcdf
	use adt_scattering_mod, only : beta_ext_adt_bin
	use numerics, only : fmin
	
    implicit none
    logical, intent(inout) :: new_file
    character (len=*),intent(in) :: outputfile
    real(wp) :: phi, sd2, sd3, deff, precip
    real(wp) :: m0_liq,m1_liq,m2_liq,m3_liq,dmean_liq,dvol_liq, &
        dsigma_liq,rel_disp_liq
    real(wp) :: m0_ice,m1_ice,m2_ice,dmean_ice,dsigma_ice,rel_disp_ice
    real(wp) :: svp1,qv,rm,rhod,beta_ext, beta_abs,beta_ext_ice, beta_abs_ice
	real(wp) :: phi_mean,nmon_mean,rhoi_mean, nact, test, &
		mcrit
	real(wp) :: denom
	real(wp) :: fallrate_liq,fallrate_ice
	integer(i4b) :: i
    
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
        
        if(.not.adiabatic_prof) then
			! define variable: radius of parcel
			call check( nf90_def_var(io1%ncid, "rad_par", NF90_DOUBLE, &
						(/io1%x_dimid/), io1%varid) )
			! get id to a_dimid
			call check( nf90_inq_varid(io1%ncid, "rad_par", io1%a_dimid) )
			! units
			call check( nf90_put_att(io1%ncid, io1%a_dimid, &
					   "units", "m") )
		endif                   




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
                   "units", "#/kg") )
                   
        ! Liquid-drop bulk size diagnostics. These are calculated from the
        ! same activated population used for ndrop and from the instantaneous
        ! wet volume-equivalent diameter.
        call check( nf90_def_var(io1%ncid, "deff", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        call check( nf90_put_att(io1%ncid,io1%varid,"units","m") )
        call check( nf90_put_att(io1%ncid,io1%varid,"long_name", &
                   "liquid effective diameter M3/M2 of activated drops") )

        call check( nf90_def_var(io1%ncid, "dmean_liq", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        call check( nf90_put_att(io1%ncid,io1%varid,"units","m") )
        call check( nf90_put_att(io1%ncid,io1%varid,"long_name", &
                   "number-weighted mean wet diameter of activated liquid drops") )

        call check( nf90_def_var(io1%ncid, "dvol_liq", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        call check( nf90_put_att(io1%ncid,io1%varid,"units","m") )
        call check( nf90_put_att(io1%ncid,io1%varid,"long_name", &
                   "volume-mean diameter (M3/M0)^(1/3) of activated liquid drops") )

        call check( nf90_def_var(io1%ncid, "rel_disp_liq", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        call check( nf90_put_att(io1%ncid,io1%varid,"units","1") )
        call check( nf90_put_att(io1%ncid,io1%varid,"long_name", &
                   "relative dispersion sigma_D/mean_D of activated liquid drops") )

        ! Native warm-particle PSD coordinates. nwat retains all warm particles;
        ! nliq contains only activated liquid drops. dwet is the instantaneous
        ! wet diameter associated with each native model bin. These are moving
        ! wet diameters, not fixed wet-size intervals, so no diagnostic wet-bin
        ! width is written; post-processing should rebin (number,dwet) onto a
        ! fixed logarithmic or instrument grid before calculating dN/dlogD.
        call check( nf90_def_var(io1%ncid, "dwet", NF90_DOUBLE, &
                    (/io1%bin2_dimid,io1%mode_dimid,io1%x_dimid/), io1%varid) )
        call check( nf90_put_att(io1%ncid,io1%varid,"units","m") )
        call check( nf90_put_att(io1%ncid,io1%varid,"long_name", &
                   "instantaneous wet volume-equivalent diameter of all warm particles") )
        call check( nf90_put_att(io1%ncid,io1%varid,"comment", &
                   "Use with nwat for the complete warm PSD; includes aerosol, haze and activated drops") )

        call check( nf90_def_var(io1%ncid, "nliq", NF90_DOUBLE, &
                    (/io1%bin2_dimid,io1%mode_dimid,io1%x_dimid/), io1%varid) )
        call check( nf90_put_att(io1%ncid,io1%varid,"units","kg-1") )
        call check( nf90_put_att(io1%ncid,io1%varid,"long_name", &
                   "activated liquid-drop number mixing ratio in each native bin") )

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
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "long_name", "number mixing ratio of all warm particles") )
        call check( nf90_put_att(io1%ncid, io1%a_dimid, "comment", &
                   "Includes unactivated aerosol, haze particles and activated "// &
                   "liquid drops; pair with dwet for the complete warm PSD") )

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

        ! define variable: precip
        call check( nf90_def_var(io1%ncid, "precip", NF90_DOUBLE, &
                    (/io1%x_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "precip", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "mm hr-1") )
                   
                   
  		if(fallout_flag) then
			! Cumulative liquid-water fallout
			call check(nf90_def_var(io1%ncid,"qfall_liq",NF90_DOUBLE, &
				(/io1%x_dimid/),io1%varid))
			call check(nf90_put_att(io1%ncid,io1%varid, &
				"units","kg kg-1"))
			call check(nf90_put_att(io1%ncid,io1%varid, &
				"long_name","cumulative liquid water removed by fallout"))
			! Cumulative liquid particle-number fallout
			call check(nf90_def_var(io1%ncid,"nfall_liq",NF90_DOUBLE, &
				(/io1%x_dimid/),io1%varid))
			call check(nf90_put_att(io1%ncid,io1%varid, &
				"units","kg-1"))
			call check(nf90_put_att(io1%ncid,io1%varid, &
				"long_name","cumulative liquid particle number removed by fallout"))
			! Actual timestep-mean liquid fallout rate
			call check(nf90_def_var(io1%ncid,"fallrate_liq",NF90_DOUBLE, &
				(/io1%x_dimid/),io1%varid))
			call check(nf90_put_att(io1%ncid,io1%varid, &
				"units","mm h-1"))
			call check(nf90_put_att(io1%ncid,io1%varid, &
				"long_name","liquid fallout rate from residence-time sink"))
        endif

        if (chamber_bl_mix.gt.0) then
            call check(nf90_def_var(io1%ncid,"qchamber_bl",NF90_DOUBLE, &
                (/io1%x_dimid/),io1%varid))
            call check(nf90_put_att(io1%ncid,io1%varid,"units","kg kg-1"))
            call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                "cumulative airborne water transferred to wall by BL saturation cap"))
            call check(nf90_def_var(io1%ncid,"qchamber_bl_step",NF90_DOUBLE, &
                (/io1%x_dimid/),io1%varid))
            call check(nf90_put_att(io1%ncid,io1%varid,"units","kg kg-1"))
            call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                "airborne water transferred to wall by BL saturation cap in previous timestep"))
        endif

        if (chamber_fan_loss.gt.0) then
            call check(nf90_def_var(io1%ncid,"qfan_liq",NF90_DOUBLE, &
                (/io1%x_dimid/),io1%varid))
            call check(nf90_put_att(io1%ncid,io1%varid,"units","kg kg-1"))
            call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                "cumulative liquid water removed by chamber fan"))
            call check(nf90_def_var(io1%ncid,"nfan_liq",NF90_DOUBLE, &
                (/io1%x_dimid/),io1%varid))
            call check(nf90_put_att(io1%ncid,io1%varid,"units","kg-1"))
            call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                "cumulative warm particle number removed by chamber fan"))
        endif

        if (chamber_wall_loss.gt.0) then
            call check(nf90_def_var(io1%ncid,"qwall_liq",NF90_DOUBLE, &
                (/io1%x_dimid/),io1%varid))
            call check(nf90_put_att(io1%ncid,io1%varid,"units","kg kg-1"))
            call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                "cumulative liquid water removed by non-gravitational chamber wall deposition"))
            call check(nf90_def_var(io1%ncid,"nwall_liq",NF90_DOUBLE, &
                (/io1%x_dimid/),io1%varid))
            call check(nf90_put_att(io1%ncid,io1%varid,"units","kg-1"))
            call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                "cumulative warm particle number removed by non-gravitational chamber wall deposition"))
        endif

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
                       "units", "#/kg") )  
                        
            ! Native ice size coordinate. dwice is maximum dimension (Dmax),
            ! not a volume-equivalent diameter. nicem contains the associated
            ! ice number mixing ratio on the same native bins.
            call check( nf90_def_var(io1%ncid, "dmaxice", NF90_DOUBLE, &
                        (/io1%bin2_dimid,io1%mode_dimid,io1%x_dimid/), io1%varid) )
            call check( nf90_put_att(io1%ncid,io1%varid,"units","m") )
            call check( nf90_put_att(io1%ncid,io1%varid,"long_name", &
                       "instantaneous maximum dimension Dmax of ice particles") )

            call check( nf90_def_var(io1%ncid, "dmean_ice", NF90_DOUBLE, &
                        (/io1%x_dimid/), io1%varid) )
            call check( nf90_put_att(io1%ncid,io1%varid,"units","m") )
            call check( nf90_put_att(io1%ncid,io1%varid,"long_name", &
                       "number-weighted mean ice maximum dimension") )

            call check( nf90_def_var(io1%ncid, "rel_disp_ice", NF90_DOUBLE, &
                        (/io1%x_dimid/), io1%varid) )
            call check( nf90_put_att(io1%ncid,io1%varid,"units","1") )
            call check( nf90_put_att(io1%ncid,io1%varid,"long_name", &
                       "relative dispersion sigma_Dmax/mean_Dmax of ice particles") )

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
                       
            if(fallout_flag) then
				call check(nf90_def_var(io1%ncid,"qfall_ice",NF90_DOUBLE, &
					(/io1%x_dimid/),io1%varid))
				call check(nf90_put_att(io1%ncid,io1%varid, &
					"units","kg kg-1"))
				call check(nf90_put_att(io1%ncid,io1%varid, &
					"long_name","cumulative ice water removed by fallout"))
				call check(nf90_def_var(io1%ncid,"nfall_ice",NF90_DOUBLE, &
					(/io1%x_dimid/),io1%varid))
				call check(nf90_put_att(io1%ncid,io1%varid, &
					"units","kg-1"))
				call check(nf90_put_att(io1%ncid,io1%varid, &
					"long_name","cumulative ice particle number removed by fallout"))
				call check(nf90_def_var(io1%ncid,"fallrate_ice",NF90_DOUBLE, &
					(/io1%x_dimid/),io1%varid))
				call check(nf90_put_att(io1%ncid,io1%varid, &
					"units","mm h-1"))
				call check(nf90_put_att(io1%ncid,io1%varid, &
					"long_name","ice fallout rate from residence-time sink"))
            endif

            if (chamber_fan_loss.gt.0) then
                call check(nf90_def_var(io1%ncid,"qfan_ice",NF90_DOUBLE, &
                    (/io1%x_dimid/),io1%varid))
                call check(nf90_put_att(io1%ncid,io1%varid,"units","kg kg-1"))
                call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                    "cumulative ice water removed by chamber fan"))
                call check(nf90_def_var(io1%ncid,"nfan_ice",NF90_DOUBLE, &
                    (/io1%x_dimid/),io1%varid))
                call check(nf90_put_att(io1%ncid,io1%varid,"units","kg-1"))
                call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                    "cumulative ice particle number removed by chamber fan"))
            endif

            if (chamber_wall_loss.gt.0) then
                call check(nf90_def_var(io1%ncid,"qwall_ice",NF90_DOUBLE, &
                    (/io1%x_dimid/),io1%varid))
                call check(nf90_put_att(io1%ncid,io1%varid,"units","kg kg-1"))
                call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                    "cumulative ice water removed by non-gravitational chamber wall deposition"))
                call check(nf90_def_var(io1%ncid,"nwall_ice",NF90_DOUBLE, &
                    (/io1%x_dimid/),io1%varid))
                call check(nf90_put_att(io1%ncid,io1%varid,"units","kg-1"))
                call check(nf90_put_att(io1%ncid,io1%varid,"long_name", &
                    "cumulative ice particle number removed by non-gravitational chamber wall deposition"))
            endif
        endif
        
        call check( nf90_enddef(io1%ncid) )

        ! write variable: mbinedges
        call check( nf90_inq_varid(io1%ncid, "mbinedges", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, parcel1%mbinedges))

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

	if(.not.adiabatic_prof) then
		! write variable: ql
		call check( nf90_inq_varid(io1%ncid, "rad_par", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, &
			parcel1%y(parcel1%ira), start = (/io1%icur/)))	
	endif
	
	svp1=svp_liq(parcel1%y(parcel1%ite))
	qv=eps1*parcel1%y(parcel1%irh)*svp1 / &
		(parcel1%y(parcel1%ipr)-svp1)
	rm=ra+qv*rv
	! Dry-air density, consistent with number/mass mixing ratios
	rhod=parcel1%y(parcel1%ipr) / &
		(rm*parcel1%y(parcel1%ite))
		
	if (use_adt_optics) then
		call beta_ext_adt_bin(parcel1%dw, parcel1%npart, rhod, &
							  optics_wavelength, .false., beta_ext, beta_abs)
	else
		! Original geometric-optics approximation: Qext = 2
		beta_ext = 2._wp*pi*rhod*sum(onequarter*parcel1%dw**2 * parcel1%npart)	
		beta_abs = 0._wp
	endif
	
	if(fallout_flag) then
		fallrate_liq=0._wp
		fallrate_ice=0._wp	
		fallrate_liq = &
			parcel1%qfall_step_liq * &
			rhod*merge(chamber_height,residence_depth, &
                n_levels_c.gt.0 .and. chamber_height.gt.0._wp)/rhow / &
			parcel1%dt * &
			1000._wp*3600._wp
		if(ice_flag.eq.1) then
			fallrate_ice = &
				parcel1%qfall_step_ice * &
				rhod*merge(chamber_height,residence_depth, &
                n_levels_c.gt.0 .and. chamber_height.gt.0._wp)/rhow / &
				parcel1%dt * &
				1000._wp*3600._wp
		endif
	endif


 
    if(ice_flag .eq. 0) then
		! write variable: beta_ext
		call check( nf90_inq_varid(io1%ncid, "beta_ext", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, beta_ext, &
					start = (/io1%icur/)))
	endif
	
    ! Refresh diagnostic particle diameters from the accepted state. This is
    ! diagnostic-only and does not change particle number or mass moments.
    call update_terminal_velocities()

    ! Diagnose activated drops once and use exactly the same population for
    ! ndrop, bulk size moments, relative dispersion and the native nliq PSD.
    parcel1%ndrop=0._wp
    do i=1,parcel1%n_bin_modew
        if (parcel1%npart(i) <= 0._wp) cycle
        if (particle_is_activated(i,parcel1%y(i), &
            parcel1%y(parcel1%ite))) parcel1%ndrop(i)=parcel1%npart(i)
    enddo
    nact=sum(parcel1%ndrop)

    call check( nf90_inq_varid(io1%ncid, "ndrop", io1%varid ) )
    call check( nf90_put_var(io1%ncid,io1%varid,nact,start=(/io1%icur/)))

    m0_liq=sum(parcel1%ndrop)
    m1_liq=sum(parcel1%ndrop*parcel1%dw)
    m2_liq=sum(parcel1%ndrop*parcel1%dw**2)
    m3_liq=sum(parcel1%ndrop*parcel1%dw**3)
    if (m0_liq > tiny(1._wp)) then
        dmean_liq=m1_liq/m0_liq
        dvol_liq=(m3_liq/m0_liq)**onethird
        dsigma_liq=sqrt(max(m2_liq/m0_liq-dmean_liq**2,0._wp))
        if (dmean_liq > tiny(1._wp)) then
            rel_disp_liq=dsigma_liq/dmean_liq
        else
            rel_disp_liq=0._wp
        endif
    else
        dmean_liq=0._wp
        dvol_liq=0._wp
        rel_disp_liq=0._wp
    endif
    if (m2_liq > tiny(1._wp)) then
        deff=m3_liq/m2_liq
    else
        deff=0._wp
    endif

    call check( nf90_inq_varid(io1%ncid, "deff", io1%varid ) )
    call check( nf90_put_var(io1%ncid,io1%varid,deff,start=(/io1%icur/)))
    call check( nf90_inq_varid(io1%ncid, "dmean_liq", io1%varid ) )
    call check( nf90_put_var(io1%ncid,io1%varid,dmean_liq,start=(/io1%icur/)))
    call check( nf90_inq_varid(io1%ncid, "dvol_liq", io1%varid ) )
    call check( nf90_put_var(io1%ncid,io1%varid,dvol_liq,start=(/io1%icur/)))
    call check( nf90_inq_varid(io1%ncid, "rel_disp_liq", io1%varid ) )
    call check( nf90_put_var(io1%ncid,io1%varid,rel_disp_liq,start=(/io1%icur/)))

    call check( nf90_inq_varid(io1%ncid, "dwet", io1%varid ) )
    call check( nf90_put_var(io1%ncid,io1%varid, &
        reshape(parcel1%dw,(/parcel1%n_bins1,parcel1%n_modes/)), &
        start=(/1,1,io1%icur/)))
    call check( nf90_inq_varid(io1%ncid, "nliq", io1%varid ) )
    call check( nf90_put_var(io1%ncid,io1%varid, &
        reshape(parcel1%ndrop,(/parcel1%n_bins1,parcel1%n_modes/)), &
        start=(/1,1,io1%icur/)))

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


	! liquid precipitation		
	precip = rhod*sum(parcel1%vel(1:parcel1%n_bin_modew)* &
		parcel1%y(1:parcel1%n_bin_modew)* &
		parcel1%npart(1:parcel1%n_bin_modew)/rhow*1000._wp*3600._wp)

	if(ice_flag .eq. 0) then
		! write variable: precip
		call check( nf90_inq_varid(io1%ncid, "precip", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, precip, &
					start = (/io1%icur/)))
	endif

	if(fallout_flag) then
		call check(nf90_inq_varid(io1%ncid,"qfall_liq",io1%varid))
		call check(nf90_put_var(io1%ncid,io1%varid, &
			parcel1%qfall_liq,start=(/io1%icur/)))

		call check(nf90_inq_varid(io1%ncid,"nfall_liq",io1%varid))
		call check(nf90_put_var(io1%ncid,io1%varid, &
			parcel1%nfall_liq,start=(/io1%icur/)))
			
		call check(nf90_inq_varid(io1%ncid,"fallrate_liq",io1%varid))
		call check(nf90_put_var(io1%ncid,io1%varid, &
			fallrate_liq,start=(/io1%icur/)))
	endif

    if (chamber_bl_mix.gt.0) then
        call check(nf90_inq_varid(io1%ncid,"qchamber_bl",io1%varid))
        call check(nf90_put_var(io1%ncid,io1%varid, &
            parcel1%qchamber_bl,start=(/io1%icur/)))
        call check(nf90_inq_varid(io1%ncid,"qchamber_bl_step",io1%varid))
        call check(nf90_put_var(io1%ncid,io1%varid, &
            parcel1%qchamber_bl_step,start=(/io1%icur/)))
    endif

    if (chamber_fan_loss.gt.0) then
        call check(nf90_inq_varid(io1%ncid,"qfan_liq",io1%varid))
        call check(nf90_put_var(io1%ncid,io1%varid, &
            parcel1%qfan_liq,start=(/io1%icur/)))
        call check(nf90_inq_varid(io1%ncid,"nfan_liq",io1%varid))
        call check(nf90_put_var(io1%ncid,io1%varid, &
            parcel1%nfan_liq,start=(/io1%icur/)))
    endif

    if (chamber_wall_loss.gt.0) then
        call check(nf90_inq_varid(io1%ncid,"qwall_liq",io1%varid))
        call check(nf90_put_var(io1%ncid,io1%varid, &
            parcel1%qwall_liq,start=(/io1%icur/)))
        call check(nf90_inq_varid(io1%ncid,"nwall_liq",io1%varid))
        call check(nf90_put_var(io1%ncid,io1%varid, &
            parcel1%nwall_liq,start=(/io1%icur/)))
    endif

    if(ice_flag .eq. 1) then
    	if (use_adt_optics) then
			! Ice: use BMM's calculated projected area
			call beta_ext_adt_bin(parcel1%dwice, parcel1%npartice, rhod, &
								  optics_wavelength, .true., &
								  beta_ext_ice, beta_abs_ice, &
								  parcel1%areaice)
    		beta_ext = beta_ext+beta_ext_ice
    		beta_abs = beta_abs+beta_abs_ice
    	else
    		! Original geometric-optics approximation: Qext = 2
			beta_ext = beta_ext + 2._wp*rhod*sum(parcel1%areaice*parcel1%npartice)
			beta_abs = 0._wp
    	endif
    	
		! write variable: beta_ext
		call check( nf90_inq_varid(io1%ncid, "beta_ext", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, beta_ext, &
					start = (/io1%icur/)))
    	
    	
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

        call check( nf90_inq_varid(io1%ncid, "dmaxice", io1%varid ) )
        call check( nf90_put_var(io1%ncid,io1%varid, &
            reshape(parcel1%dwice,(/parcel1%n_bins1,parcel1%n_modes/)), &
            start=(/1,1,io1%icur/)))

        m0_ice=sum(parcel1%npartice)
        m1_ice=sum(parcel1%npartice*parcel1%dwice)
        m2_ice=sum(parcel1%npartice*parcel1%dwice**2)
        if (m0_ice > tiny(1._wp)) then
            dmean_ice=m1_ice/m0_ice
            dsigma_ice=sqrt(max(m2_ice/m0_ice-dmean_ice**2,0._wp))
            if (dmean_ice > tiny(1._wp)) then
                rel_disp_ice=dsigma_ice/dmean_ice
            else
                rel_disp_ice=0._wp
            endif
        else
            dmean_ice=0._wp
            rel_disp_ice=0._wp
        endif
        call check( nf90_inq_varid(io1%ncid, "dmean_ice", io1%varid ) )
        call check( nf90_put_var(io1%ncid,io1%varid,dmean_ice,start=(/io1%icur/)))
        call check( nf90_inq_varid(io1%ncid, "rel_disp_ice", io1%varid ) )
        call check( nf90_put_var(io1%ncid,io1%varid,rel_disp_ice,start=(/io1%icur/)))
    
        call check( nf90_inq_varid(io1%ncid, "mice", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            reshape(parcel1%yice(1:parcel1%n_bin_modew),(/n_bins,n_mode/)), &
             start = (/1,1,io1%icur/)))

        ! write variable: phi
		denom=sum(parcel1%moments( &
			parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
			parcel1%n_comps+2))			
		if(denom > tiny(1._wp)) then
			phi_mean=sum(parcel1%moments( &
				parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
				parcel1%n_comps+1))/denom
		else
			phi_mean=0._wp
		endif
        call check( nf90_inq_varid(io1%ncid, "phi", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            phi_mean, start = (/io1%icur/)))
    
        ! write variable: nmon
		denom=sum(parcel1%npartice)		
		if(denom > tiny(1._wp)) then
			nmon_mean=sum(parcel1%moments( &
				parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
				parcel1%n_comps+2))/denom
		else
			nmon_mean=0._wp
		endif
        call check( nf90_inq_varid(io1%ncid, "nmon", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            nmon_mean, start = (/io1%icur/)))
    
        ! write variable: rhoi
		denom=sum(parcel1%moments( &
			parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
			parcel1%n_comps+3))		
		if(denom > tiny(1._wp)) then
			rhoi_mean = &
				(sum(parcel1%yice*parcel1%npartice) - &
				 sum(parcel1%moments( &
					 parcel1%n_bin_modew+1:parcel1%n_bin_mode, &
					 parcel1%n_comps+4))) / denom
		else
			rhoi_mean=0._wp
		endif
        call check( nf90_inq_varid(io1%ncid, "rhoi", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            rhoi_mean, start = (/io1%icur/)))
            
        call check( nf90_inq_varid(io1%ncid, "nicem", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            reshape(parcel1%npartice(1:parcel1%n_bin_modew), &
            (/parcel1%n_bins1,parcel1%n_modes/)), start = (/1,1,io1%icur/)))

        call check( nf90_inq_varid(io1%ncid, "maeri", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            reshape(parcel1%mbinice(:,&
                1:parcel1%n_comps), &
            (/parcel1%n_bins1,parcel1%n_modes,parcel1%n_comps/)), start = (/1,1,1,io1%icur/)))

            
		! ice precipitation
		precip = precip + rhod*sum( &
			parcel1%vel(parcel1%n_bin_modew+1:parcel1%n_bin_mode) * &
			parcel1%yice(1:parcel1%n_bin_modew) * &
			parcel1%npartice(1:parcel1%n_bin_modew) / &
			rhow*1000._wp*3600._wp)
    	
		! write variable: precip
		call check( nf90_inq_varid(io1%ncid, "precip", io1%varid ) )
		call check( nf90_put_var(io1%ncid, io1%varid, precip, &
					start = (/io1%icur/)))

		if(fallout_flag) then					
			call check(nf90_inq_varid(io1%ncid,"qfall_ice",io1%varid))
			call check(nf90_put_var(io1%ncid,io1%varid, &
				parcel1%qfall_ice,start=(/io1%icur/)))
			
			call check(nf90_inq_varid(io1%ncid,"nfall_ice",io1%varid))
			call check(nf90_put_var(io1%ncid,io1%varid, &
				parcel1%nfall_ice,start=(/io1%icur/)))
			
			call check(nf90_inq_varid(io1%ncid,"fallrate_ice",io1%varid))
			call check(nf90_put_var(io1%ncid,io1%varid, &
				fallrate_ice,start=(/io1%icur/)))
		endif

        if (chamber_fan_loss.gt.0) then
            call check(nf90_inq_varid(io1%ncid,"qfan_ice",io1%varid))
            call check(nf90_put_var(io1%ncid,io1%varid, &
                parcel1%qfan_ice,start=(/io1%icur/)))
            call check(nf90_inq_varid(io1%ncid,"nfan_ice",io1%varid))
            call check(nf90_put_var(io1%ncid,io1%varid, &
                parcel1%nfan_ice,start=(/io1%icur/)))
        endif

        if (chamber_wall_loss.gt.0) then
            call check(nf90_inq_varid(io1%ncid,"qwall_ice",io1%varid))
            call check(nf90_put_var(io1%ncid,io1%varid, &
                parcel1%qwall_ice,start=(/io1%icur/)))
            call check(nf90_inq_varid(io1%ncid,"nwall_ice",io1%varid))
            call check(nf90_put_var(io1%ncid,io1%varid, &
                parcel1%nwall_ice,start=(/io1%icur/)))
        endif
            
    endif
    

    call check( nf90_close(io1%ncid) )


    io1%icur=io1%icur+1
    end subroutine output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	! ============================================================================
	! map_to_sce
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Copies the current BMM liquid and optional ice populations into the combined arrays used by the
	!>SCE routines.
	!>@param[in] ice_flag: switch indicating whether ice populations are included
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

	! ============================================================================
	! update_collision_kernel
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Recomputes the full symmetric collision kernel from current liquid/ice dimensions, projected
	!>areas, masses, terminal velocities and air properties.
	subroutine update_collision_kernel()	
		use sce, only : collision_air_properties, &
						collision_kernel_pair, &
						long_collection_efficiency_pair

		implicit none
		integer(i4b) :: i,j
		integer(i4b) :: n,nall
		real(wp) :: eff,kernel
		real(wp) :: va,prefac,lambda_air
		real(wp), dimension(parcel1%n_bin_mode) :: dcoll
		real(wp), dimension(parcel1%n_bin_mode) :: acoll
		real(wp), dimension(parcel1%n_bin_mode) :: dlong
		real(wp), dimension(parcel1%n_bin_mode) :: mtot
	
		n=parcel1%n_bin_modew
		nall=parcel1%n_bin_mode
	
		! ----------------------------------------------------------------------
		! Current collision dimensions
		! ----------------------------------------------------------------------
		! Liquid: current wet volume-equivalent diameter
		dcoll(1:n)=parcel1%dw
		acoll(1:n)=pi*onequarter*parcel1%dw**2

		! For Long efficiency, use actual liquid-drop diameter
		dlong(1:n)=parcel1%dw
		if (parcel1%ice_flag.eq.1) then
			! Ice physical characteristic size
			dcoll(n+1:nall)=parcel1%dwice		
			! Ice projected area
			acoll(n+1:nall)=parcel1%areaice
			
			! Melted-equivalent diameter for applying the temporary
			! liquid-drop Long collection efficiency to ice.
			!
			! Use ICE MASS ONLY, not total mass including aerosol.
			dlong(n+1:nall) = (6._wp*max( &
					parcel1%mbinall(n+1:nall,parcel1%n_comps+1), &
					0._wp) / (pi*rhow))**onethird
			
		endif
		! ----------------------------------------------------------------------
		! Total particle mass
		!
		! mbinall columns:
		!   1:n_comps      aerosol components
		!   n_comps+1      water or ice mass
		! ----------------------------------------------------------------------
		mtot=sum(parcel1%mbinall(:,1:parcel1%n_comps+1),dim=2)
		! ----------------------------------------------------------------------
		! Air properties are identical for every collision pair.
		! Calculate them only ONCE.
		! ----------------------------------------------------------------------
		call collision_air_properties( &
			parcel1%y(parcel1%ite), &
			parcel1%y(parcel1%ipr), &
			va,prefac,lambda_air)
	
		! ----------------------------------------------------------------------
		! Collision kernel.
		!
		! Only calculate one triangle because K(i,j) = K(j,i).
		! ----------------------------------------------------------------------
		do j=1,nall
			do i=1,j
	
				! -------------------------------------------------------------
				! Gravitational collection efficiency
				!
				! FOR NOW:
				! retain the old Long treatment for every pair while checking
				! that the code refactor itself works.
				!
				! The phase-specific LL / LI / II efficiencies should be
				! introduced here afterwards.
				! -------------------------------------------------------------
				if ((i <= n).and.(j <= n)) then
					! Liquid-liquid
					eff=long_collection_efficiency_pair( &
						dlong(i),dlong(j))
				else
					! TEMPORARY mixed/ice behaviour during refactor validation.
					!
					! Do not regard this as the final ice collection efficiency.
					eff=long_collection_efficiency_pair( &
						dlong(i),dlong(j))
				endif
				! -------------------------------------------------------------
				! Complete pair kernel:
				!
				! gravitational
				! turbulent inertial
				! turbulent shear
				! Brownian
				! Brownian enhancement
				! -------------------------------------------------------------
				call collision_kernel_pair( &
					parcel1%y(parcel1%ite), &
					dcoll(i), &
					dcoll(j), &
					acoll(i), &
					acoll(j), &
					mtot(i), &
					mtot(j), &
					parcel1%vel(i), &
					parcel1%vel(j), &
					parcel1%nre(i), &
					parcel1%nre(j), &
					eff, &
					va,prefac,lambda_air, &
					kernel)
	
				parcel1%ecoll(i,j)=kernel
				! Kernel is symmetric
				parcel1%ecoll(j,i)=kernel
	
			enddo
		enddo
	end subroutine update_collision_kernel
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! ============================================================================
	! map_to_bmm
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Copies the combined post-SCE liquid and optional ice populations back into the native BMM arrays
	!>and solution vectors.
	!>@param[in] ice_flag: switch indicating whether ice populations are included
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


	! ============================================================================
	! adjust_t_rh
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Applies latent heat of freezing/melting to temperature and recomputes relative humidity at fixed
	!>vapour mixing ratio.
	!>@param[in] totmass: net water mass undergoing the phase change represented by the latent-heating
	!>adjustment
	!>@param[inout] t: parcel temperature
	!>@param[inout] rh: parcel relative humidity
	!>@param[in] p: parcel pressure
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



	! ============================================================================
	! bmm_driver
	! ============================================================================
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>Main BMM timestep driver coordinating output, diffusional microphysics, SCE collection/secondary
	!>ice, terminal velocities, fallout and entrainment.
	!>@param[in] sce_flag: switch selecting stochastic collection
	!>@param[in] hm_flag: switch for Hallett-Mossop secondary-ice production
	!>@param[in] break_flag: selector for collisional ice-breakup treatment
	!>@param[in] mode1_flag: switch for mode-1 secondary-ice production
	!>@param[in] mode2_flag: switch for mode-2 secondary-ice production
    subroutine bmm_driver(sce_flag,hm_flag,break_flag,mode1_flag, mode2_flag)
    use numerics_type
    use sce, only : sce_microphysics, sce_sip_microphysics, qsmall
    implicit none
    integer(i4b), intent(in) :: sce_flag
    logical, intent(in) :: hm_flag, mode1_flag, mode2_flag
    integer(i4b), intent(in) :: break_flag
    integer(i4b) :: i, j, nt
    real(wp) :: rhoa
    
	! calculate terminal velocities (needed because of the output call on 1st time-step)
	call update_terminal_velocities()

    nt=ceiling(runtime / real(dt,kind=wp))
    do i=1,nt
        ! Output state left by previous timestep
        call output(io1%new_file,outputfile)
        
        if ((updraft_type==2).and.(parcel1%TT>t_thresh)) parcel1%y(parcel1%iw)=0._wp


		! Condensation, deposition, nucleation, etc.
		! These calculate temporary terminal velocities internally
		! when ventilation is required.
        call bin_microphysics(fparcelwarm, fparcelcold, & 
            icenucleation, noncollisional_iceformation)

        ! Optional chamber boundary-layer recirculation.  Both BL modes use
        ! the same external total-water exchange; only the phase response differs.
        if (chamber_bl_mix.gt.0) then
            call apply_chamber_bl_exchange()
        endif

		! Particle masses/shapes have now changed.  BIN_FULL_MOVING is
		! deliberately NOT projected to the fixed SCE grid here; the SCE
		! moving-pivot gain treatment updates its representative masses
		! directly.
		call update_terminal_velocities()

        if(sce_flag.gt.0) then
        	! calculate the density of air
        	rhoa=parcel1%y(parcel1%ipr)/(parcel1%y(parcel1%ite)*ra)
            
            ! Map the BMM variables to the SCE variables
            call map_to_sce(ice_flag)
			! Build the collision kernel from the CURRENT BMM state
			call update_collision_kernel()              
              
            ! One timestep of the selected SCE model.  Both the basic SCE
            ! (sce_flag=1) and the secondary-ice SCE (sce_flag=2) use the
            ! fully moving gain treatment when bin_scheme_flag=BIN_FULL_MOVING.
            if(sce_flag.eq.1) then
                call sce_microphysics(parcel1%n_bins1,parcel1%n_bin_mode,&
                	parcel1%n_bin_modew, parcel1%n_comps+parcel1%imoms,&
                                parcel1%npartall,parcel1%moments,parcel1%momenttype, &
                                parcel1%ecoll,parcel1%indexc, &
                                parcel1%mbinall(:,parcel1%n_comps+1),parcel1%dt, &
                                parcel1%y(parcel1%ite),rhoa,parcel1%totaddto, &
                                parcel1%bin_scheme_flag.eq.BIN_FULL_MOVING)
            elseif(sce_flag.eq.2) then
                call sce_sip_microphysics(parcel1%n_bins1,parcel1%n_bin_mode,&
                            parcel1%n_bin_modew,parcel1%n_comps, parcel1%n_comps+&
                            parcel1%imoms,&
                            parcel1%npartall,parcel1%moments,parcel1%momenttype, &
                            parcel1%ecoll,parcel1%indexc, &
                            parcel1%mbinall(:,parcel1%n_comps+1),parcel1%vel,parcel1%dt, &
                            parcel1%y(parcel1%ite),rhoa, parcel1%totaddto, &
                            mass_fragment1, mass_fragment2, mass_fragment3, &
                            hm_flag,break_flag,mode1_flag, mode2_flag, &
                            parcel1%bin_scheme_flag.eq.BIN_FULL_MOVING )
            endif

			! latent heat of fusion
			if(ice_flag.eq.1) then
				if(.not.chamber_force_temperature) then
					call adjust_t_rh(parcel1%totaddto,parcel1%y(parcel1%ite), &
							parcel1%y(parcel1%irh), parcel1%y(parcel1%ipr))
				endif
				parcel1%yice(parcel1%irhi)=parcel1%y(parcel1%irh)
				parcel1%yice(parcel1%itei)=parcel1%y(parcel1%ite)
			endif
                                        
            ! redefine the mass of each component of aerosol
            do j=1,parcel1%n_comps
                where (parcel1%npartall(:).gt.qsmall)
                    parcel1%mbinall(:,j)=parcel1%moments(:,j)/parcel1%npartall(:)
                end where
            enddo                                                         
            ! map SCE to BMM
            call map_to_bmm(ice_flag)

			! SCE has changed masses/numbers/composition, so update
			! velocities again for the new accepted state.
			call update_terminal_velocities()
            
        endif    

        !--------------------------------------------------------------
        ! Optional chamber-specific physical particle losses.  These are
        ! distinct from gravitational sedimentation/fallout below.
        !--------------------------------------------------------------
        if (chamber_fan_loss.gt.0 .or. chamber_wall_loss.gt.0) then
            call apply_chamber_particle_losses()
        endif

        !--------------------------------------------------------------
        ! Sedimentation / finite parcel residence time
        !
        ! If SCE is off, velocities are those calculated after
        ! bin_microphysics.
        !
        ! If SCE is on, velocities have just been recalculated after SCE.
        !--------------------------------------------------------------
        if(fallout_flag) then
            call apply_particle_fallout()
        endif


        ! --------------------------------------------------------------
        ! Prescribed cloud-top stopping criterion.  z_ctop is a generic
        ! parcel run-control limit and is independent of the removed
        ! Sanchez cloud-top entrainment pathway.  Check it only after
        ! completing the accepted timestep so condensation/deposition,
        ! SCE and fallout remain operator-consistent.  The final output
        ! call below writes this last accepted state.
        ! --------------------------------------------------------------
        if ((.not.chamber_forcing_active) .and. (z_ctop.gt.0._wp)) then
            if (parcel1%y(parcel1%iz).ge.z_ctop) then
                parcel1%break_flag=.true.
            endif
        endif

        ! break-out if flag has been set
        if(parcel1%break_flag) exit
    enddo
    ! output to file
    call output(io1%new_file,outputfile)
    
    
    end subroutine bmm_driver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end module bmm	

