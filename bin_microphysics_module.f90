	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>code to allocate arrays, and call activation 
	module bmm
    use nrtype
    use nr, only : locate, polint
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the bin microphysics model

    implicit none
        ! constants for the bin microphysics model

        real(sp), parameter :: r_gas=8.314_sp, molw_a=29.e-3_sp,molw_water=18.e-3_sp, &
                                cp=1005.0_sp, cpv=1870._sp, cpw=4.27e3_sp, cpi=2104.6_sp, &
                                grav=9.81_sp, &
        						lv=2.5e6_sp, ls=2.837e6_sp, lf=ls-lv, ttr=273.15_sp, &
        						joules_in_an_erg=1.0e-7_sp,joules_in_a_cal=4.187e0_sp, &
        						rhow=1000._sp, ra=r_gas/molw_a,rv=r_gas/molw_water , &
        						eps1=ra/rv, rhoice=910._sp
        						

        type parcel
            ! variables for bin model
            integer(i4b) :: n_bins1,n_modes,n_comps, n_bin_mode, n_bin_mode1, &
                            n_sound, ice_flag
            real(sp) :: dt
            real(sp), dimension(:,:), allocatable :: q_sound
            real(sp), dimension(:), allocatable :: t_sound, z_sound, rh_sound, &
                                                    p_sound, theta_q_sound
            real(sp) :: z,p,t,w,rh, qinit, t_cbase, q_cbase, p_cbase, z_cbase, &
                        t_ctop, q_ctop, p_ctop, z_ctop, theta_q_cbase, theta_q_ctop, &
                        x_ent, theta_q
                        
                        
            ! liquid water
            real(sp), dimension(:), allocatable :: d, maer, npart, rho_core, &
                            rh_eq, rhoat, dw, da_dt, ndrop
            real(sp), dimension(:,:), allocatable :: mbin, rhobin, &
                                        nubin,molwbin,kappabin ! all bins x all comps                                
            ! variables for ODE:                    
            integer(i4b) :: neq, itol, ipr, ite, iz, iw, irh, &
                            itask, istate, iopt, mf, lrw, liw
            integer(i4b), dimension(:), allocatable :: iwork, ipar
            real(sp) :: rtol, tt, tout
            real(sp), dimension(:), allocatable :: y, yold, atol, rwork, rpar
            
            
            ! ice water
            real(sp), dimension(:), allocatable :: dice, maerice, npartice, rho_coreice, &
                            rh_eqice, rhoatice, dwice, da_dtice, nice, &
                            phi, rhoi, nump, rime
            real(sp), dimension(:,:), allocatable :: mbinice, rhobinice, &
                                        nubinice,molwbinice,kappabinice ! all bins x all comps                                
            ! variables for ODE:                    
            integer(i4b) :: neqice, itolice, ipri, itei, izi, iwi, irhi, &
                            itaskice, istateice, ioptice, mfice, lrwice, liwice
            integer(i4b), dimension(:), allocatable :: iworkice, iparice
            real(sp) :: rtolice, ttice, toutice
            real(sp), dimension(:), allocatable :: yice, yoldice, atolice, rworkice, rparice
            
            
            logical :: break_flag=.false.
            
        end type parcel

        type sounding
            ! variables for grid
            integer(i4b) :: n_levels
            real(sp), dimension(:,:), allocatable :: q
            real(sp), dimension(:), allocatable :: theta, p, z, rh
        end type sounding


        type io
            ! variables for io
            integer(i4b) :: ncid, varid, x_dimid, bin_dimid, mode_dimid, y_dimid, z_dimid, &
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
        real(sp) :: ent_rate, dmina,dmaxa
        real(sp) :: zinit,tpert,winit,tinit,pinit,rhinit,z_ctop, alpha_therm, alpha_cond, &
                    alpha_therm_ice, alpha_dep
        integer(i4b) :: microphysics_flag=0, kappa_flag,updraft_type, vent_flag, &
                        ice_flag=0, bin_scheme_flag=1
        logical :: use_prof_for_tprh
        real(sp) :: dz,dt, runtime
        ! sounding spec
        real(sp) :: psurf, tsurf
        integer(i4b), parameter :: nlevels_r=1000
        integer(i4b), parameter :: nq=3
        integer(i4b) :: n_levels_s, idum, n_sel
        real(sp) :: mult, rh_act
        real(sp), allocatable, dimension(:,:) :: q_read !nq x nlevels_r
        real(sp), allocatable, dimension(:) :: theta_read,rh_read,  z_read
        ! aerosol setup
        integer(i4b) :: n_intern, n_mode,n_sv,sv_flag,n_bins,n_comps
        ! aerosol_spec
        real(sp), allocatable, dimension(:,:) :: n_aer1,d_aer1,sig_aer1, mass_frac_aer1
        real(sp), allocatable, dimension(:) ::  molw_core1,density_core1,nu_core1, &
                                        kappa_core1
        real(sp), allocatable, dimension(:) :: org_content1, molw_org1, kappa_org1, &
                                    density_org1, delta_h_vap1,nu_org1,log_c_star1
        
        

        ! variables for model
        real(sp) :: theta_surf,theta_init, &
            theta_q_sat,t1old, p111, w_cb, n_dummy, d_dummy, x2old=1.0_sp
        logical :: set_theta_q_cb_flag=.true.


        character (len=200) :: outputfile='output'

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
		use nrtype
		implicit none
		integer(i4b), intent(in) :: n_intern, n_mode, n_sv, n_bins,n_comps, nq, &
		                            n_levels_s
		real(sp), dimension(:), allocatable, intent(inout) :: theta_read,rh_read,z_read, &
		                        org_content1,molw_org1,kappa_org1, &
		                        density_org1,delta_h_vap1,nu_org1,log_c_star1
		real(sp), dimension(:,:), allocatable, intent(inout) :: q_read, &
		                        n_aer1,d_aer1,sig_aer1,mass_frac_aer1
		real(sp), dimension(:), allocatable, intent(inout) :: molw_core1,density_core1, &
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
        namelist /run_vars/ outputfile, runtime, dt, &
                    zinit,tpert,use_prof_for_tprh,winit,tinit,pinit,rhinit, &
                    microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                    kappa_flag, updraft_type,adiabatic_prof, vert_ent, &
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
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine read_in_bmm_namelist







	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! initialise arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>interpolates the sounding to the grid
    subroutine initialise_bmm_arrays()
    use nrtype
    use nr, only : locate, polint, rkqs, odeint, zbrent, brent

    implicit none
    real(sp) :: num, ntot, number_per_bin, test, var1, &
                eps2, z1, z2, htry, hmin, var, dummy
    real(sp), dimension(1) :: p1, z11
    real(sp) :: p11, p22, rm, cpm
    integer(i4b) :: i,j,k, AllocateStatus, iloc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set variables and allocate arrays in parcel                                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parcel1%n_sound=n_levels_s   
    parcel1%n_bins1=n_bins    
    parcel1%n_modes=n_mode
    parcel1%n_comps=n_comps
    parcel1%n_bin_mode=n_bins*n_mode
    parcel1%n_bin_mode1=(n_bins+1)*n_mode
    parcel1%z=zinit
    parcel1%p=pinit
    parcel1%t=tinit
    parcel1%w=winit
    parcel1%rh=rhinit
    parcel1%dt=dt
    allocate( parcel1%d(1:parcel1%n_bin_mode1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%maer(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%npart(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%mbin(1:parcel1%n_bin_mode,1:n_comps+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%rho_core(1:parcel1%n_modes), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

    allocate( parcel1%rhobin(1:parcel1%n_bin_mode,1:n_comps), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%nubin(1:parcel1%n_bin_mode,1:n_comps), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%molwbin(1:parcel1%n_bin_mode,1:n_comps), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%kappabin(1:parcel1%n_bin_mode,1:n_comps), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

    allocate( parcel1%rh_eq(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%rhoat(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%dw(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%da_dt(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%ndrop(1:parcel1%n_bin_mode), STAT = AllocateStatus)
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
        parcel1%rho_core(i) = 1._sp/var1
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
        number_per_bin=ntot/real(n_bins,sp)
        parcel1%npart(1+(k-1)*n_bins:(k)*n_bins)=number_per_bin
        parcel1%d(1+(k-1)*(n_bins+1))=dmina
        do i=1,n_bins
            d_dummy=parcel1%d(i+(k-1)*(n_bins+1))
            n_dummy=number_per_bin*(1._sp-1.e-15_sp)
            parcel1%d(i+1+(k-1)*(n_bins+1))= zbrent(find_upper_diameter, &
                        d_dummy*0.9_sp,dmaxa*2._sp,1.e-30_sp)
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
                pi/6._sp*(0.5_sp*(parcel1%d(i+1)+parcel1%d(i)))**3 * &
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
    parcel1%t_sound(1)=theta_read(1)*(psurf/1.e5_sp)**(ra/cp) 
    parcel1%p_sound(1)=psurf
    
    eps2=1.e-5_sp
    htry=10._sp
    hmin=1.e-2_sp
    p1=psurf

    ! integrate hydrostatic equation
    do i=2,n_levels_s
        call odeint(p1,z_read(i-1),z_read(i),eps2,htry,hmin,hydrostatic1b,rkqs)
        parcel1%t_sound(i)=theta_read(i)*(p1(1)/1.e5_sp)**(ra/cp)  
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
        iloc=locate(parcel1%z_sound(1:n_levels_s),parcel1%z)
        iloc=min(n_levels_s-1,iloc)
        iloc=max(1,iloc)
        ! linear interp t
        call polint(parcel1%z_sound(iloc:iloc+1), parcel1%t_sound(iloc:iloc+1), &
                    min(parcel1%z,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%t=var +tpert
        ! linear interp rh
        call polint(parcel1%z_sound(iloc:iloc+1), parcel1%p_sound(iloc:iloc+1), &
                    min(parcel1%z,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%p=var 
        ! linear interp rh
        call polint(parcel1%z_sound(iloc:iloc+1), parcel1%rh_sound(iloc:iloc+1), &
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
		theta_init=parcel1%t*(1.e5_sp/parcel1%p)**(ra/cp)
		! calculate p required so that qs(t,p) = q init
		parcel1%p_cbase=zbrent(cloud_base,parcel1%p, 100._sp, 1.e-30_sp)
		rm=ra+parcel1%qinit*rv
		cpm=cp+parcel1%qinit*cpv
		parcel1%t_cbase=theta_init*(parcel1%p_cbase/1.e5_sp)**(rm/cpm)
		parcel1%theta_q_cbase= &
		    calc_theta_q3(parcel1%t_cbase,parcel1%p_cbase,parcel1%q_cbase)
		print *,'Cloud-base t, p: ',parcel1%t_cbase,parcel1%p_cbase
		
		
		
		! now, find the height of cloud-base along dry adiabat:
		theta_surf=tsurf*(1.e5_sp/psurf)**(ra/cp)
		p11=parcel1%p
		z11(1)=parcel1%z
		p22=parcel1%p_cbase
		htry=p22-p11
		eps2=1.e-5_sp
		call odeint(z11,p11,p22,eps2,htry,hmin,hydrostatic1,rkqs)
		parcel1%z_cbase=z11(1)

        ! cloud-top properties from sounding:
        parcel1%z_ctop=z_ctop
        ! interpolate to find t
        iloc=locate(parcel1%z_sound(1:n_levels_s),parcel1%z_ctop)
        iloc=min(n_levels_s-1,iloc)
        iloc=max(1,iloc)
        ! linear interp t
        call polint(parcel1%z_sound(iloc:iloc+1), parcel1%t_sound(iloc:iloc+1), &
                    min(parcel1%z_ctop,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%t_ctop=var
        ! linear interp p
        call polint(parcel1%z_sound(iloc:iloc+1), parcel1%p_sound(iloc:iloc+1), &
                    min(parcel1%z_ctop,parcel1%z_sound(n_levels_s)), var,dummy)        
        parcel1%p_ctop=var
        ! linear interp q
        call polint(parcel1%z_sound(iloc:iloc+1), parcel1%rh_sound(iloc:iloc+1), &
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
            do i=1,parcel1%n_bin_mode
                n_sel=i
                rh_act=0._sp !min(parcel1%rh,0.999_sp)
                mult=-1._sp
                ! has to be less than the peak moles of water at activation
                d_dummy=brent(1.e-50_sp,1.e-12_sp,1.e1_sp, koehler02,1.e-30_sp,test)
                rh_act=min(parcel1%rh,0.999_sp)
                mult=1._sp
                d_dummy=zbrent(koehler02,1.e-30_sp, test, 1.e-30_sp)*molw_water 
                parcel1%mbin(i,n_comps+1)= d_dummy
            enddo
!             call koehler01(parcel1%t,parcel1%mbin(:,n_comps+1),&
!                 parcel1%mbin,parcel1%rhobin,&
!                 parcel1%nubin,parcel1%molwbin,parcel1%n_bin_mode,&
!                 parcel1%rh_eq,parcel1%rhoat,parcel1%dw) 
!             print *,parcel1%rh_eq
        case(1)
            do i=1,parcel1%n_bin_mode
                n_sel=i
                rh_act=0._sp !min(parcel1%rh,0.999_sp)
                mult=-1._sp
                ! has to be less than the peak moles of water at activation
                d_dummy=brent(1.e-50_sp,1.e-12_sp,1.e1_sp, kkoehler02,1.e-30_sp,test)
                rh_act=min(parcel1%rh,0.999_sp)
                mult=1._sp
                d_dummy=zbrent(kkoehler02,1.e-30_sp, test, 1.e-30_sp)*molw_water 
                parcel1%mbin(i,n_comps+1)= d_dummy
            enddo
!             call kkoehler01(parcel1%t,parcel1%mbin(:,n_comps+1),&
!                 parcel1%mbin,parcel1%rhobin,&
!                 parcel1%kappabin,parcel1%molwbin,parcel1%n_bin_mode,&
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
    parcel1%neq=parcel1%n_bin_mode+5 ! p,t,rh,z,w
    parcel1%tt=0._sp
    parcel1%tout=parcel1%tt+parcel1%dt
    parcel1%itol=2
    parcel1%rtol=1.e-4_sp
    allocate( parcel1%y(parcel1%neq), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"
    allocate( parcel1%yold(parcel1%neq), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"
    allocate( parcel1%atol(parcel1%neq), stat = allocatestatus)
    if (allocatestatus /= 0) stop "*** not enough memory ***"
    
    parcel1%atol(1:parcel1%n_bin_mode)=1.e-25_sp
    
    parcel1%ipr=parcel1%n_bin_mode+1 ! pressure
    parcel1%ite=parcel1%n_bin_mode+2 ! temperarture
    parcel1%irh=parcel1%n_bin_mode+3 ! rh
    parcel1%iz =parcel1%n_bin_mode+4 ! altitude
    parcel1%iw =parcel1%n_bin_mode+5 ! vertical wind
    
    parcel1%atol(parcel1%ipr)=10._sp
    parcel1%atol(parcel1%ite)=1.e-4_sp
    parcel1%atol(parcel1%irh)=1.e-8_sp
    parcel1%atol(parcel1%iz) =2.e-2_sp
    parcel1%atol(parcel1%iw) =2.e-2_sp
    
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
    parcel1%iwork(6) = 100 ! max steps
    parcel1%iwork(7) = 10 ! max message printed per problem
    parcel1%iwork(5) = 5 ! order
    parcel1%rwork(5) = 1.e-3_sp ! initial time-step
    parcel1%rwork(6) = dt ! max time-step
    parcel1%rwork(7) = 1.e-9_sp ! min time-step allowed
    parcel1%rwork(14) = 2._sp ! tolerance scale factor
    
    ! put water in solution vector and set p, t, rh, z, w
    parcel1%y(1:parcel1%n_bin_mode)=parcel1%mbin(:,n_comps+1)
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
        allocate( parcel1%maerice(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%npartice(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%mbinice(1:parcel1%n_bin_mode,1:n_comps+1), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%rho_coreice(1:parcel1%n_modes), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

        allocate( parcel1%rhobinice(1:parcel1%n_bin_mode,1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%nubinice(1:parcel1%n_bin_mode,1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%molwbinice(1:parcel1%n_bin_mode,1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%kappabinice(1:parcel1%n_bin_mode,1:n_comps), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

        allocate( parcel1%rh_eqice(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%rhoatice(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%dwice(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%da_dtice(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%nice(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        
        allocate( parcel1%phi(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%rhoi(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        allocate( parcel1%nump(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"     
        allocate( parcel1%rime(1:parcel1%n_bin_mode), STAT = AllocateStatus)
        if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        
        parcel1%phi=1._sp
        parcel1%rhoi=rhoice
        parcel1%nump=1._sp
        parcel1%rime=0._sp
        
                
        parcel1%rho_coreice(:) = parcel1%rho_core(:)
        
        parcel1%npartice=0._sp
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
        parcel1%neqice=parcel1%n_bin_mode+4 ! p,t,rh,w
        parcel1%ttice=0._sp
        parcel1%toutice=parcel1%tout
        parcel1%itolice=2
        parcel1%rtolice=1.e-3_sp
        allocate( parcel1%yice(parcel1%neqice), stat = allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory ***"
        allocate( parcel1%yoldice(parcel1%neqice), stat = allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory ***"
        allocate( parcel1%atolice(parcel1%neqice), stat = allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory ***"
    
        parcel1%atolice(1:parcel1%n_bin_mode)=1.e-25_sp
    
        parcel1%ipri=parcel1%n_bin_mode+1 ! pressure
        parcel1%itei=parcel1%n_bin_mode+2 ! temperature
        parcel1%irhi=parcel1%n_bin_mode+3 ! rh
        parcel1%iwi =parcel1%n_bin_mode+4 ! vertical wind
    
        parcel1%atolice(parcel1%ipri)=10._sp
        parcel1%atolice(parcel1%itei)=1.e-4_sp
        parcel1%atolice(parcel1%irhi)=1.e-8_sp
        parcel1%atolice(parcel1%iwi) =2.e-2_sp
    
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
        parcel1%iworkice(6) = 100 ! max steps
        parcel1%iworkice(7) = 10 ! max message printed per problem
        parcel1%iworkice(5) = 5 ! order
        parcel1%rworkice(5) = 1.e-3_sp ! initial time-step
        parcel1%rworkice(6) = dt ! max time-step
        parcel1%rworkice(7) = 1.e-9_sp ! min time-step allowed
        parcel1%rworkice(14) = 2._sp ! tolerance scale factor
    
        ! put water in solution vector and set p, t, rh, z, w
        parcel1%yice(1:parcel1%n_bin_mode)=parcel1%mbinice(:,n_comps+1)
        parcel1%yice(parcel1%ipri)=parcel1%p
        parcel1%yice(parcel1%itei)=parcel1%t
        parcel1%yice(parcel1%irhi)=parcel1%rh
        parcel1%yice(parcel1%iwi) =parcel1%w
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end subroutine initialise_bmm_arrays
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
    real(sp), intent(in) :: dmin,dmax
    real(sp), intent(in), dimension(n_intern) :: n_aer1,d_aer1,sig_aer1
    integer(i4b), intent(in) :: n_intern
    real(sp), intent(inout) :: num
    
    integer(i4b) :: i
       
    num=0._sp                                 
    do i=1,n_intern
        num=num+n_aer1(i)*(0.5_sp*erfc(-log(dmax/d_aer1(i))/sqrt(2._sp)/sig_aer1(i) ) - &
            0.5_sp*erfc(-log(dmin/d_aer1(i))/sqrt(2._sp)/sig_aer1(i) ))
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
        use nrtype
        implicit none
        real(sp), intent(in) :: x
        real(sp) :: find_upper_diameter, num
        
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
	use nrtype
	implicit none
	real(sp), intent(in) :: p
	real(sp) :: t,qs, cloud_base, rm, cpm
	
	rm=ra+rv*parcel1%qinit
	cpm=cp+cpv**parcel1%qinit
	
	t=theta_init*(p/1.e5_sp)**(rm/cpm)
	qs=eps1*svp_liq(t)/(p-svp_liq(t))
	
	cloud_base=qs-parcel1%q_cbase
	
	end function cloud_base
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	subroutine hydrostatic1(p,z,dzdp)
	use nrtype
	implicit none
	real(sp), intent(in) :: p
	real(sp), dimension(:), intent(in) :: z
	real(sp), dimension(:), intent(out) :: dzdp
	real(sp) :: t
	
	t=theta_surf*(p/1.e5_sp)**(ra/cp)
	dzdp(1)=-(ra*t) / (grav*p)
	
	end subroutine hydrostatic1

	subroutine hydrostatic1a(z,p,dpdz)
	use nrtype
	implicit none
	real(sp), intent(in) :: z
	real(sp), dimension(:), intent(in) :: p
	real(sp), dimension(:), intent(out) :: dpdz
	real(sp) :: t
	
	t=theta_surf*(p(1)/1.e5_sp)**(ra/cp)
	dpdz(1)=-(grav*p(1)) / (ra*t) 
	
	end subroutine hydrostatic1a

	subroutine hydrostatic1b(z,p,dpdz)
	use nrtype
	implicit none
	real(sp), intent(in) :: z
	real(sp), dimension(:), intent(in) :: p
	real(sp), dimension(:), intent(out) :: dpdz
	real(sp) :: t, var, dummy, theta
	integer(i4b) :: iloc
	
	! interpolate to find theta
    iloc=locate(z_read(1:n_levels_s),z)
    iloc=min(n_levels_s-1,iloc)
    iloc=max(1,iloc)
    ! linear interp theta
    call polint(z_read(iloc:iloc+1), theta_read(iloc:iloc+1), &
                min(z,z_read(n_levels_s)), var,dummy)
    theta=var     
                
	t=theta*(p(1)/1.e5_sp)**(ra/cp)
	dpdz(1)=-(grav*p(1)) / (ra*t) 
	
	end subroutine hydrostatic1b

	subroutine hydrostatic2(p,z,dzdp)
	use nrtype
	use nr, only : zbrent
	implicit none
	real(sp), intent(in) :: p
	real(sp), dimension(:), intent(in) :: z
	real(sp), dimension(:), intent(out) :: dzdp
	real(sp) :: t
	
	p111=p
	t=theta_surf*(p111/1.e5_sp)**(ra/cp)
	t=zbrent(calc_theta_q,t,t1old*1.01_sp,1.e-5_sp)
!	print *,'hi',t,calc_theta_q(t)
	! find the temperature by iteration
	dzdp(1)=-(ra*t) / (grav*p)
	
	end subroutine hydrostatic2

	subroutine hydrostatic2a(z,p,dpdz)
	use nrtype
	use nr, only : zbrent
	implicit none
	real(sp), intent(in) :: z
	real(sp), dimension(:), intent(in) :: p
	real(sp), dimension(:), intent(out) :: dpdz
	real(sp) :: t
	
	p111=p(1)
	t=theta_surf*(p111/1.e5_sp)**(ra/cp)
	t=zbrent(calc_theta_q,t,t1old*1.01_sp,1.e-5_sp)
!	print *,'hi',t,calc_theta_q(t)
	! find the temperature by iteration
	dpdz(1)=-(grav*p(1))/(ra*t)
	
	end subroutine hydrostatic2a
	

	function calc_theta_q(t111)
	use nrtype
	implicit none
	real(sp), intent(in) :: t111
	real(sp) :: calc_theta_q
	real(sp) :: qs,rm,cpm
	qs=eps1*svp_liq(t111)/(p111-svp_liq(t111))
	rm=ra+rv*qs
	cpm=cp+cpv*qs
	calc_theta_q=t111*(1.e5_sp/p111)**(rm/cpm)*exp(lv*qs/cpm/t111)-theta_q_sat

	end function calc_theta_q     

	function calc_theta_q2(p)
	use nrtype
	implicit none
	real(sp), intent(in) :: p
	real(sp) :: calc_theta_q2
	real(sp) :: ws
	ws=eps1*svp_liq(t1old)/(p-svp_liq(t1old))
	calc_theta_q2=t1old*(1e5_sp/p)**(ra/cp)*exp(lv*ws/cp/t1old)-theta_q_sat

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
	use nrtype
	implicit none
	real(sp), intent(in) :: t,p,q
	real(sp) :: calc_theta_q3
	real(sp) :: qs, rm, cpm, rh
	qs=eps1*svp_liq(t)/(p-svp_liq(t))
	rm=ra+rv*q
	cpm=cp+cpv*qs !+cpw*max(q-qs,0._sp)
	rh=q/qs
	calc_theta_q3=t*(1.e5_sp/p)**(rm/cpm)*exp(lv*q/(t*(cpm))) * &
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
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: svp_liq
		svp_liq = 100._sp*6.1121_sp* &
			  exp((18.678_sp - (t-ttr)/ 234.5_sp)* &
			  (t-ttr)/(257.14_sp + (t-ttr)))
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
		use nrtype
		implicit none
		real(sp), intent(in) :: t
		real(sp) :: svp_ice
		svp_ice = 100._sp*6.1115_sp* &
            exp((23.036_sp - (t-ttr)/ 333.7_sp)* &
            (t-ttr)/(279.82_sp + (t-ttr)))

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
      use nrtype
      implicit none
      real(sp), intent(in) :: t, p
      real(sp) :: dd, t1
      t1=max(t,200._sp)
      dd=2.11e-5_sp*(t1/ttr)**1.94_sp*(101325._sp/p)
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
      use nrtype
      implicit none
      real(sp), intent(in) :: t
      real(sp) :: ka, t1
      t1=max(t,200._sp)
      ka=(5.69_sp+0.017_sp*(t1-ttr))*1.e-3_sp*joules_in_a_cal
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
        use nrtype
        implicit none
        real(sp), intent(in) :: t
        real(sp) :: viscosity_air
        real(sp) :: tc

        tc = t-ttr
        tc = max(tc,-200._sp)

        if( tc.ge.0._sp) then
            viscosity_air = (1.718_sp+0.0049_sp*tc) * 1e-5_sp ! the 1d-5 converts from poise to si units
        else
            viscosity_air = (1.718_sp+0.0049_sp*tc-1.2e-5_sp*tc**2) * 1e-5_sp
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
        use nrtype
        real(sp), intent(in) :: t
        real(sp) :: tc, surface_tension

        tc=t-ttr
        tc = max(tc,-40._sp)

        ! pruppacher and klett pg 130 
        surface_tension = 75.93_sp + 0.115_sp * tc + 6.818e-2_sp * tc**2 + &
                          6.511e-3_sp * tc**3 + 2.933e-4_sp * tc**4 + &
                          6.283e-6_sp * tc**5 + 5.285e-8_sp * tc**6
        if(tc.ge.0._sp) then
            surface_tension = 76.1_sp - 0.155_sp*tc
        end if
    
        surface_tension = surface_tension*joules_in_an_erg ! convert to j/cm2 
        surface_tension = surface_tension*1.e4_sp ! convert to j/m2 

    !    surface_tension=72d-3
        !sigma = 75.93_sp * joules_in_an_erg*1d4
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
      use nrtype
      implicit none
      real(sp), dimension(:), intent(in) :: mwat
      real(sp), dimension(:,:), intent(in) :: mbin,rhobin
      integer(i4b), intent(in) :: sz
      real(sp), dimension(:),intent(inout) :: dw
      
      real(sp), dimension(sz) :: rhoat

      ! calculate the diameter and radius
      rhoat(:)=mwat(:)/rhow+sum(mbin(:,1:n_comps)/rhobin(:,:),2)
      rhoat(:)=(mwat(:)+sum(mbin(:,1:n_comps),2))/rhoat(:);
  
      ! wet diameter:
      dw(:)=((mwat(:)+sum(mbin(:,1:n_comps),2))*6._sp/(pi*rhoat(:)))**(1._sp/3._sp)
      
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
      use nrtype
      implicit none
      real(sp), dimension(:), intent(in) :: mwat
      real(sp), dimension(:,:), intent(in) :: mbin,rhobin,nubin,molwbin
      integer(i4b), intent(in) :: sz
      real(sp), dimension(sz) :: nw
      real(sp), dimension(:),intent(inout) :: rh_eq,rhoat, dw
      real(sp), intent(in) :: t
      real(sp) :: sigma

      ! calculate the diameter and radius
      nw(:)=mwat(:)/molw_water
      rhoat(:)=mwat(:)/rhow+sum(mbin(:,1:n_comps)/rhobin(:,:),2)
      rhoat(:)=(mwat(:)+sum(mbin(:,1:n_comps),2))/rhoat(:);
  
      ! wet diameter:
      dw(:)=((mwat(:)+sum(mbin(:,1:n_comps),2))*6._sp/(pi*rhoat(:)))**(1._sp/3._sp)
  
      ! calculate surface tension
      sigma=surface_tension(t)

      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      rh_eq(:)=exp(4._sp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))* &
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
      use nrtype
      implicit none
      real(sp), dimension(:), intent(in) :: mwat
      real(sp), dimension(:,:), intent(in) :: mbin,rhobin,kappabin,molwbin
      integer(i4b), intent(in) :: sz
      real(sp), dimension(sz) :: nw,dd,kappa
      real(sp), dimension(:),intent(inout) :: rh_eq,rhoat, dw
      real(sp), intent(in) :: t
      real(sp) :: sigma
      ! calculate the diameter and radius
      nw(:)=mwat(:)/molw_water
      rhoat(:)=mwat(:)/rhow+sum(mbin(:,1:n_comps)/rhobin(:,:),2)
      rhoat(:)=(mwat(:)+sum(mbin(:,1:n_comps),2))/rhoat(:);
  
      ! wet diameter:
      dw(:)=((mwat(:)+sum(mbin(:,1:n_comps),2))* 6._sp/(pi*rhoat(:)))**(1._sp/3._sp)
  
      dd(:)=((sum(mbin(:,1:n_comps)/rhobin(:,:),2))*6._sp/(pi))**(1._sp/3._sp) ! dry diameter
                                  ! needed for eqn 6, petters and kreidenweis (2007)
  
      kappa(:)=sum(mbin(:,1:n_comps)/rhobin(:,:)*kappabin(:,:),2) &
               / sum(mbin(:,1:n_comps)/rhobin(:,:),2)
               ! equation 7, petters and kreidenweis (2007)

      ! calculate surface tension
      sigma=surface_tension(t)

      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      rh_eq(:)=exp(4._sp*molw_water*sigma/r_gas/t/rhoat(:)/dw(:))* &
           (dw(:)**3-dd(:)**3)/(dw(:)**3-dd(:)**3*(1._sp-kappa))
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
      use nrtype
      implicit none
      real(sp), intent(in) :: nw
      real(sp) :: massw
      real(sp) :: rhoat, dw,koehler02
      real(sp) :: sigma
      
      

      ! wet diameter:
      massw=nw*molw_water
      rhoat=massw/rhow+sum(parcel1%mbin(n_sel,1:n_comps) / &
            parcel1%rhobin(n_sel,1:n_comps))
      rhoat=(massw+parcel1%maer(n_sel))/rhoat;
      dw=((massw+parcel1%maer(n_sel))* 6._sp/(pi*rhoat))**(1._sp/3._sp)
  
      ! calculate surface tension
      sigma=surface_tension(parcel1%t)
  
      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      koehler02=mult*(exp(4._sp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
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
      use nrtype
      implicit none
      real(sp), intent(in) :: nw
      real(sp) :: massw
      real(sp) :: rhoat, dw,dd,kappa,kkoehler02
      real(sp) :: sigma

      ! wet diameter:
      massw=nw*molw_water
      rhoat=massw/rhow+sum(parcel1%mbin(n_sel,1:n_comps) / &
            parcel1%rhobin(n_sel,1:n_comps))
      rhoat=(massw+parcel1%maer(n_sel))/rhoat;
      dw=((massw+parcel1%maer(n_sel))* 6._sp/(pi*rhoat))**(1._sp/3._sp)
  
      ! calculate surface tension
      sigma=surface_tension(parcel1%t)
  
      dd=(sum(parcel1%mbin(n_sel,1:n_comps) / parcel1%rhobin(n_sel,:))* &
          6._sp/(pi))**(1._sp/3._sp) ! dry diameter
                                  ! needed for eqn 6, petters and kreidenweis (2007)
  
      kappa=sum(parcel1%mbin(n_sel,1:n_comps) / parcel1%rhobin(n_sel,1:n_comps)* &
               parcel1%kappabin(n_sel,:)) &
               / sum(parcel1%mbin(n_sel,1:n_comps) / parcel1%rhobin(n_sel,1:n_comps))
               ! equation 7, petters and kreidenweis (2007)

      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      kkoehler02=mult*(exp(4._sp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           (dw**3-dd**3)/(dw**3-dd**3*(1._sp-kappa)))-rh_act
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
      use nrtype
      implicit none
      real(sp), intent(in) :: mbin ! mass of aerosol particle
      real(sp) :: massw
      real(sp) :: rhoat, dw,nw,koehler03
      real(sp) :: sigma
      
      ! calculate the diameter and radius
      nw=d_dummy/molw_water ! moles of water
      rhoat=d_dummy/rhow+mbin* sum(mass_frac_aer1(n_sel,1:n_comps)/ &
                           density_core1(1:n_comps))
      rhoat=(d_dummy+mbin)/rhoat;
      dw=((d_dummy+mbin)* 6._sp/(pi*rhoat))**(1._sp/3._sp)
  
      ! calculate surface tension
      sigma=surface_tension(parcel1%t)
  
      ! equilibrium rh over particle - nb rh_act set to zero if not root-finding
      koehler03=mult*(exp(4._sp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
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
      use nrtype
      implicit none
      real(sp), intent(in) :: mbin ! mass of aerosol particle
      real(sp) :: massw
      real(sp) :: rhoat, dw,nw,dd,kappa, kkoehler03
      real(sp) :: sigma
      ! calculate the diameter and radius
      nw=d_dummy/molw_water ! moles of water
      rhoat=d_dummy/rhow+mbin* sum(mass_frac_aer1(n_sel,1:n_comps)/ &
                           density_core1(1:n_comps))
      rhoat=(d_dummy+mbin)/rhoat;
      dw=((d_dummy+mbin)* 6._sp/(pi*rhoat))**(1._sp/3._sp)
  
      dd=((sum(mbin*mass_frac_aer1(n_sel,1:n_comps)/ density_core1(1:n_comps),1))* &
          6._sp/(pi))**(1._sp/3._sp) ! dry diameter
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
      kkoehler03=mult*(exp(4._sp*molw_water*sigma/r_gas/parcel1%t/rhoat/dw)* &
           (dw**3-dd**3)/(dw**3-dd**3*(1._sp-kappa)))-rh_act
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
      use nrtype
      implicit none
      real(sp), intent(in) :: t, p
      real(sp), dimension(:), intent(in) :: diam, rhoat
      integer(i4b), intent(in) :: sz
      real(sp), dimension(sz), intent(inout) :: nre,cd,vel 
      real(sp) :: eta, sigma, physnum, phys6, mfpath, tc, rhoa
      real(sp), dimension(sz) :: diam2, mass, bondnum, bestnm, x,y


        tc=t-ttr
        vel = 0._sp ! zero array
        rhoa = p / (ra * t) ! density of air
        diam2=diam ! temporary array that can be changed

        eta = viscosity_air(t)

        nre = 0._sp ! zero array
        where(diam2.gt.7000.e-6_sp)
            diam2=7000.e-6_sp
        end where
        mass=pi/6._sp*diam2**3*rhoat
    
        sigma = surface_tension(t)
        
        ! regime 3:  eqns 5-12, 10-146 & 10-148 from p & k 
        physnum = (sigma**3._sp) * (rhoa**2_sp) / ((eta**4._sp) * grav * (rhow - rhoa))		
        phys6 = physnum**(1._sp / 6._sp)
        where(diam2.gt.1070.e-6_sp) 
            bondnum = (4._sp/3._sp)*grav * (rhow - rhoa) * (diam2**2) / sigma

            x = log(bondnum*phys6)
            y = -5.00015_sp + 5.23778_sp * x - 2.04914_sp * x * x + 0.475294_sp * (x**3) &
                - 0.542819e-1_sp * (x**4._sp) + 0.238449e-2_sp * (x**5)

            nre = phys6 * exp(y)

            vel = eta * (nre)/ (rhoa * diam2)

            cd = 8._sp * mass * grav * rhoa/(pi * ((diam2 / 2._sp)* eta)**2)
            cd = cd	/ (nre**2) 
        end where

        ! regime 2:  eqns 10-142, 10-145 & 10-146 from p & k 
        where(diam2.le.1070.e-6_sp.and.diam2.gt.20.e-6_sp)
            bestnm = 32._sp * ((diam2 / 2._sp)**3) * (rhow - rhoa) * rhoa * &
                      grav / (3._sp * eta**2)
            x = log(bestnm)
            y = -3.18657_sp + 0.992696_sp * x - 0.153193e-2_sp * x * x &
                -0.987059e-3_sp * (x**3) - 0.578878e-3_sp * (x**4) &
                + 0.855176e-4_sp * (x**5) - 0.327815e-5_sp * (x**6)
            nre =  exp(y)
            vel = eta * nre / (2._sp * rhoa * (diam2 / 2._sp))
            cd = bestnm/(nre**2)
        end where

        ! regime 1:  eqns 10-138, 10-139 & 10-140 from p & k 
        mfpath = 6.6e-8_sp * (101325_sp / p) * (t / 293.15_sp)
        where(diam2.le.20.e-6_sp) 
            vel = 2._sp * ((diam2 / 2._sp)**2) * grav * (rhow - rhoa) / (9._sp * eta)
            vel = vel * (1._sp + 1.26_sp * mfpath / (diam2 / 2._sp))
            nre = vel * rhoa * diam2 / eta

            cd = 8._sp * mass * grav * rhoa/(pi * ((diam2 / 2._sp)* eta)**2)
            cd = cd	/ (nre**2) 
         end where

        where(isnan(vel))
          vel=0._sp
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
      use nrtype
      implicit none
      real(sp), intent(in) :: t, p
      real(sp), dimension(:), intent(in) :: diam, rhoat
      real(sp), dimension(sz), intent(inout) :: fv,fh
      integer(i4b), intent(in) :: sz
      real(sp), dimension(sz) :: nre,cd,vel,calc
      real(sp) :: d1,k1,rhoa, eta, nu, nsc1,nsc2
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
      calc = (nsc1**(1._sp/3._sp)) * sqrt(nre)
      where(calc.gt.51.4_sp)
        calc=51.4_sp
      end where

      where(calc.lt.1.4_sp)
        fv=1.00_sp+0.108_sp*calc**2
      elsewhere
        fv=0.78_sp+0.308_sp*calc
      end where
      !-----------------------------------
    
      ! heat ventilation - use ka---------
      calc = (nsc2**(1._sp/3._sp)) * sqrt(nre)
      where(calc.gt.51.4_sp)
        calc=51.4_sp
      end where

      where(calc.lt.1.4_sp)
        fh=1.00_sp+0.108_sp*calc**2
      elsewhere
        fh=0.78_sp+0.308_sp*calc
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
      use nrtype
      implicit none
      real(sp), intent(in) :: t, p
      real(sp), dimension(:), intent(in) :: mwat, phi, rhoi, nump, rime
      integer(i4b), intent(in) :: sz
      real(sp), dimension(sz), intent(inout) :: nre,vel 
      real(sp) :: eta, rhoa
      real(sp), dimension(sz) :: dmax,drime,area,ar, x
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
!       ar=area/(pi/4._sp* (dmax**2._sp))
!       ar=min(max(ar,0.1_sp),1._sp)
      ar=1._sp
  
      ! heymsfield and westbrook
      x=rhoa*8._sp*mwat*grav/( (eta**2._sp)*pi*(ar**0.5_sp))
      nre=(8.0_sp**2._sp)/4._sp* &
          ( (sqrt(1._sp+(4._sp*sqrt(x))/( (8._sp**2._sp)*sqrt(0.35_sp)))-1._sp)**2._sp)
      vel=eta*(nre)/(rhoa*dmax)
    
      ! viscous regime
      where(nre.lt.1._sp) 
        vel = grav*mwat / (6._sp*pi*eta*0.465_sp*dmax*(ar**0.5_sp))
      end where

      where(isnan(vel)) 
            vel=0._sp
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
      use nrtype

      implicit none
      real(sp), intent(in) :: t, p
      real(sp), dimension(:), intent(in) :: mwat, phi, rhoi, nump, rime
      real(sp), dimension(sz), intent(inout) :: fv,fh
      integer(i4b), intent(in) :: sz
      real(sp), dimension(sz) :: nre,vel,nre2,x
      real(sp) :: d1,k1,rhoa, eta, nu, nsc1,nsc2, calc1, calc2
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
      calc1 = nsc1**(1._sp/3._sp)
      calc2 = nsc2**(1._sp/3._sp)
  
      ! columns
      nre2=min(nre,20._sp)
      where(phi.gt.1.0_sp)  
        x = calc1*sqrt(nre2)	
        fv = 1.0_sp - 0.000668_sp*x/4._sp + 2.39402_sp*((x/4._sp)**2._sp) + &
             0.73409_sp*((x/4._sp)**3._sp)-0.73911_sp*((x/4._sp)**4._sp)
        x = calc2*sqrt(nre2);	
        fh = 1.0_sp - 0.000668_sp*x/4._sp + 2.39402_sp*((x/4._sp)**2._sp) + &
             0.73409_sp*((x/4._sp)**3._sp)-0.73911_sp*((x/4._sp)**4._sp)
      end where
      !--------
  
      ! plates
      nre2=min(nre,120._sp)
      where(phi.le.1._sp) 
        x = calc1*sqrt(nre2)	
        fv = 1.0_sp - 0.06042_sp*x/10._sp + 2.79820_sp*((x/10._sp)**2._sp) - &
             0.31933_sp*((x/10._sp)**3._sp)-0.06247_sp*((x/10._sp)**4._sp)
        x = calc2*sqrt(nre2)	
        fh = 1.0_sp - 0.06042_sp*x/10._sp + 2.79820_sp*((x/10._sp)**2._sp) - &
             0.31933_sp*((x/10._sp)**3._sp)-0.06247_sp*((x/10._sp)**4._sp)
      end where
      !-------
  
      ! broad-branched crystals
      !nre2=min(nre,120d0) ! already done above
      where(phi.lt.0.2_sp.and.rhoi.le.500._sp) 
        x = calc1*sqrt(nre2)	
        fv = 1.0_sp + 0.35463_sp*x/10._sp + 3.55333_sp*((x/10._sp)**2._sp)
        x = calc2*sqrt(nre2)
        fh = 1.0_sp + 0.35463_sp*x/10._sp + 3.55333_sp*((x/10._sp)**2._sp)
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
      use nrtype
      implicit none
      real(sp), intent(in) :: t, p, rh
      real(sp), dimension(:), intent(in) :: rh_eq,rhoat, diam
      integer(i4b), intent(in) :: sz
      real(sp), dimension(sz) :: dropgrowthrate01
      real(sp), dimension(sz) :: rad, dstar,kstar,fv,fh
      real(sp) :: d1,k1,rhoa
  
      rad=diam/2._sp
      ! density of air
      rhoa=p/ra/t
      ! diffusivity of water vapour in air
      d1=dd(t,p)
      ! thermal conductivity of air
      k1=ka(t)
      ! ventilation coefficient
      fv=1._sp
      fh=1._sp
      if(vent_flag.eq.1) then
        call ventilation01(diam, rhoat,t, p, fv, fh,sz)
      end if
  
      ! modify diffusivity and conductivity
      dstar=d1*fv/(rad/(rad+0.7_sp*8.e-8_sp)+d1*fv/rad/alpha_cond*sqrt(2._sp*pi/rv/t))
      kstar=k1*fh/(rad/(rad+2.16e-7_sp)+k1*fh/rad/alpha_therm/cp/rhoa*sqrt(2._sp*pi/ra/t))
  
      ! 455 jacobson and 511 pruppacher and klett
      dropgrowthrate01=dstar*lv*rh_eq*svp_liq(t)* &
                       rhoat/kstar/t*(lv*molw_water/t/r_gas-1._sp) 
      dropgrowthrate01=dropgrowthrate01+rhoat*r_gas*t/molw_water  
      dropgrowthrate01=dstar*(rh-rh_eq)*svp_liq(t)/rad/dropgrowthrate01
                 
    end function dropgrowthrate01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate growth rate of an ice crystal					  				   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function icegrowthrate01(t,p,rh_ice,rh_eq,mwat,mbin,rhobin,phi, rhoi, nump,rime,sz) 
      use nrtype

      implicit none
      real(sp), intent(in) :: t, p, rh_ice
      real(sp), dimension(:), intent(in) :: rh_eq,mwat, phi, rhoi, nump, rime
      real(sp), dimension(:,:), intent(in) :: mbin,rhobin
      integer(i4b), intent(in) :: sz
      real(sp), dimension(sz) :: icegrowthrate01
      real(sp), dimension(sz) :: rad, dstar,kstar,rhoat,diam,fv,fh
      real(sp) :: d1,k1,rhoa
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
      fv=1._sp
      fh=1._sp
      if(vent_flag.eq.1) then
        call ventilation02(mwat, t, p,phi, &
             rhoi,nump, rime,fv, fh,sz)
      end if
      ! modify diffusivity and conductivity
      dstar=d1*fv/(rad/(rad+0.7_sp*8e-8_sp)+d1*fv/rad/alpha_dep*sqrt(2._sp*pi/rv/t))
      kstar=k1*fh/(rad/(rad+2.16e-7_sp)+k1*fh/rad/alpha_therm_ice/cp/rhoa*sqrt(2._sp*pi/ra/t))
  
      ! 473 jacobson 
      icegrowthrate01=dstar*ls*rh_eq*svp_ice(t)/ &
                       kstar/t*(ls*molw_water/t/r_gas-1._sp) 
      icegrowthrate01=icegrowthrate01+r_gas*t/molw_water  
      icegrowthrate01=4._sp*pi*rad*dstar*(rh_ice-rh_eq)*svp_ice(t)/icegrowthrate01
                 
    end function icegrowthrate01
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the capacitance of an ice crystal				  				   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculates the capacitance of oblate (longer a) and prolate (longer c) 
    ! spheroids. a is the length of the prism axis (half of) and c the basal 
    ! (half of). see page 1214 of chen and lamb, jas, 1994.
    function capacitance01(mwat,phi,rhoice,numi,rimemass,sz)
      use nrtype
      implicit none
      real(sp), dimension(:), intent(in) :: mwat, phi, rhoice, numi, rimemass
      real(sp), dimension(sz) :: capacitance01,vol,a,c,ecc,dmax,drime
  
      integer(i4b), intent(in) :: sz
  
      vol=mwat/rhoice
  
      a=( 3._sp*vol/(4._sp*pi*phi) )**(1._sp/3._sp)
      c=a*phi
  
      where(phi.lt.1._sp)
        ecc=sqrt(1._sp-phi**2._sp)
        capacitance01=a*ecc/asin(ecc)
      elsewhere
        ecc=sqrt(1._sp-phi**(-2._sp))
        capacitance01=c*ecc/log((1._sp+ecc)*phi)
      end where
    
      where(abs(phi-1._sp).lt.1.e-4_sp)
        capacitance01=a
      end where
      ! westbrook et al. (2008, jas): capacitance of aggregates is 0.25 times the 
      ! maximum dimension of the ice particle
!       call maxdimension01(mwat-rimemass,rhoice,phi,numi,rimemass,dmax,drime,sz)
!       where(numi.ge.2._sp)
!          capacitance01=0.25_sp*dmax
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
          use nrtype
          implicit none
          integer(i4b), intent(in) :: sz
          real(sp), parameter :: k_water_amb=1.6_sp, &
                    dk_water_dp=-8.8_sp,k_ice_amb=0.22_sp,dk_ice_dp=-0.17_sp
          real(sp) :: pg,t,p, & 
                      integral3, mudiff0
          real(sp), intent(in), dimension(:) :: aw
          real(sp), dimension(sz) :: koopnucrate,deltaaw,logj
      
          pg=p/1.e9_sp

      
          integral3=(-230.76_sp - 0.1478_sp * t + 4099.2_sp * t**(-1) + &
                     48.8341_sp * log(t) ) * &
            (pg - 0.5_sp * (k_water_amb + dk_water_dp * pg)* pg**2 &
             - (1._sp/6._sp) * dk_water_dp * pg**3 ) &
             - (19.43_sp - 2.2e-3_sp * t + 1.08e-5_sp * t**2 ) * &
            (pg - 0.5_sp * (k_ice_amb + dk_ice_dp * pg) * pg**2 - &
            (1._sp/6._sp) * dk_ice_dp * pg**3 )

          mudiff0 = 210368._sp + 131.438_sp * t - 3.32373e6_sp * t**(-1)  &
                   - 41729.1_sp * log(t)
    
          ! delta activity
          deltaaw = aw * exp(integral3 / r_gas / t) - exp(mudiff0 / r_gas / t)
        

          ! nucleation rate
          logj = -906.7_sp + 8502._sp * deltaaw - 26924._sp * deltaaw**2 + 29180._sp * &
                deltaaw**3

          koopnucrate = (10._sp**logj) * 1.e6_sp;	! nucleation rate in m^-3 s^-1 
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
        use nrtype
        use nr, only : dfridr,locate

        implicit none
        real(sp), intent(inout) :: tt
        real(sp), intent(inout), dimension(neq) :: y, ydot
        integer(i4b), intent(inout) :: neq
        real(sp), intent(inout) :: rpar
        integer(i4b), intent(inout) :: ipar

        ! local variables
        real(sp) :: wv=0._sp, wl=0._sp, wi=0._sp, rm, cpm, &
                  drv=0._sp, dri=0._sp,dri2=0._sp, &
                  rh,t,p,err,sl, w, &
                  te, qve, pe, var, dummy, rhoe, rhop, b

        integer(i4b) :: i, j,iloc, ipart, ipr, ite, irh, iz,iw

        ipart=parcel1%n_bin_mode
        ipr=parcel1%ipr
        ite=parcel1%ite
        irh=parcel1%irh
        iz =parcel1%iz
        iw =parcel1%iw

        rh=y(irh)
        t=y(ite)
        p=y(ipr)
        w=y(iw)
    

        ! check there are no negative values
        where(y(1:ipart).le.0.e1_sp)
            y(1:ipart)=1.e-22_sp
        end where


        ! calculate mixing ratios from rh, etc
        sl=svp_liq(t)*rh/(p-svp_liq(t)) ! saturation ratio
        sl=(sl*p/(1._sp+sl))/svp_liq(t)
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
            if(isnan(parcel1%da_dt(i)).or.parcel1%npart(i).le. 1.e-9_sp) then
              parcel1%da_dt(i)=0._sp
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
            iloc=locate(parcel1%z_sound(1:n_levels_s),y(iz))
            iloc=min(n_levels_s-1,iloc)
            iloc=max(1,iloc)
            ! linear interp p
            call polint(parcel1%z_sound(iloc:iloc+1), parcel1%p_sound(iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            pe=var
            ! linear interp qv
            call polint(parcel1%z_sound(iloc:iloc+1), parcel1%q_sound(1,iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            qve=var
            ! linear interp te
            call polint(parcel1%z_sound(iloc:iloc+1), parcel1%t_sound(iloc:iloc+1), &
                        min(y(iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            te=var
            ! env density:
            rhoe=pe/(rm*te)
            ! parcel density:
            rhop=p/(rm*t)
            !buoyancy
            if((parcel1%z_sound(n_levels_s) .lt. y(iz)) .or. &
                (parcel1%z_sound(1) .gt. y(iz))) then
                b=0._sp
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
!             parcel1%x_ent=max(0._sp,parcel1%x_ent)
!             parcel1%x_ent=min(1._sp,parcel1%x_ent)
!             
!             ydot(ite)=ydot(ite)+&
!                 parcel1%x_ent/dt*(parcel1%t_ctop-y(ite) + &
!                 lv/cpm*(parcel1%q_ctop-wv)) +&
!                 (1._sp-parcel1%x_ent)/dt*(parcel1%t_cbase-y(ite) + &
!                 lv/cpm*(parcel1%q_cbase-wv))
!         endif 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! change in rh of parcel                                                 !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ydot(irh)=(p-svp_liq(t))*svp_liq(t)*drv
        ydot(irh)=ydot(irh)+svp_liq(t)*wv*ydot(ipr)
        ydot(irh)=ydot(irh)-wv*p*dfridr(svp_liq,t,1.e0_sp,err)*ydot(ite)
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
        use nrtype
        use nr, only : dfridr,locate

        implicit none
        real(sp), intent(inout) :: tt
        real(sp), intent(inout), dimension(neq) :: y, ydot
        integer(i4b), intent(inout) :: neq
        real(sp), intent(inout) :: rpar
        integer(i4b), intent(inout) :: ipar

        ! local variables
        real(sp) :: wv=0._sp, wl=0._sp, wi=0._sp, rm, cpm, &
                  drv=0._sp, dri=0._sp,dri2=0._sp, &
                  rh,t,p,err,sl, w, &
                  te, qve, pe, var, dummy, rhoe, rhop, b, rh_ice

        integer(i4b) :: i, j,iloc, ipartice, ipr, ite, irh, iz,iw

        ipartice=parcel1%n_bin_mode
        ipr=parcel1%ipri
        ite=parcel1%itei
        irh=parcel1%irhi
        iw =parcel1%iwi
        
        ydot(iw)=0._sp

        rh=y(irh)
        t=y(ite)
        p=y(ipr)
        w=y(iw)
    

        ! check there are no negative values
        where(y(1:ipartice).le.0.e1_sp)
            y(1:ipartice)=1.e-22_sp
        end where


        ! calculate mixing ratios from rh, etc
        sl=svp_liq(t)*rh/(p-svp_liq(t)) ! saturation ratio
        sl=(sl*p/(1._sp+sl))/svp_liq(t)
        wv=eps1*rh*svp_liq(t) / (p-svp_liq(t)) ! vapour mixing ratio
        wl=sum(parcel1%npart*parcel1%y(1:ipartice))          ! liquid mixing ratio
        wi=sum(parcel1%npartice*y(1:ipartice))             ! liquid mixing ratio
        rh_ice = wv / ( eps1*svp_ice(t) / (p-svp_ice(t) ) ) ! rh over ice

        ! calculate the moist gas constants and specific heats
        rm=ra+wv*rv
        cpm=cp+wv*cpv+wl*cpw+wi*cpi

        ! now calculate derivatives
        ! adiabatic parcel model
        ydot(ipr)=0._sp      ! hydrostatic equation

        
        ! particle growth rate - mass growth rate
        parcel1%rh_eq=1._sp
        ydot(1:ipartice)=icegrowthrate01(t,p,rh_ice,parcel1%rh_eq,y(1:ipartice), &
            parcel1%mbinice(:,1:n_comps),parcel1%rhobinice,&
            parcel1%phi,parcel1%rhoi,parcel1%nump,parcel1%rime,ipartice) 
        
        ! do not bother if number concentration too small
        do i=1,ipartice
            if(isnan(ydot(i)).or.parcel1%npartice(i).le. 1.e-9_sp) then
              ydot(i)=0._sp
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
        ydot(irh)=ydot(irh)-wv*p*dfridr(svp_liq,t,1.e0_sp,err)*ydot(ite)
        ydot(irh)=ydot(irh) / (eps1*svp_liq(t)**2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
    end subroutine fparcelcold
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! jacobian for a warm parcel model : dummy subroutine                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine jparcelwarm(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
          use nrtype
          implicit none
          real(sp) :: t
          real(sp), dimension(neq) :: y
          real(sp), dimension(nrpd, neq) :: pd
          integer(i4b) :: neq, ml, mu, nrpd
          real(sp) :: rpar
          integer(i4b) :: ipar
      
    end subroutine jparcelwarm
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ice nucleation                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine icenucleation(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin,t,p,sz,sz2,yice,rh,dt) 
      use nrtype
      implicit none
      real(sp), intent(inout) :: t
      real(sp), intent(in) :: p,rh,dt
      real(sp), dimension(sz2), intent(inout) :: npart,npartice
      real(sp), dimension(:), intent(in) :: mwat
      real(sp), dimension(:,:), intent(in) :: mbin2, &
                                              rhobin,nubin,kappabin,molwbin
      integer(i4b), intent(in) :: sz,sz2
      real(sp), dimension(sz2) :: nw,aw,jw,dn01,m01,ns,dw,dd,kappa,rhoat
      real(sp), dimension(sz2,sz) :: dmaer01
      real(sp), dimension(sz2,sz), intent(inout) :: mbin2_ice
      
      real(sp), intent(inout), dimension(sz2) :: yice
      
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
  
          dw(:)=((mwat(:)+sum(mbin2(:,:),2))*6._sp/(pi*rhoat(:)))**(1._sp/3._sp)
  
          dd(:)=((sum(mbin2(:,:)/rhobin(:,:),2))* &
                 6._sp/(pi))**(1._sp/3._sp) ! dry diameter
                              ! needed for eqn 6, petters and kreidenweis (2007)
          kappa(:)=sum((mbin2(:,:)+1.e-60_sp)/rhobin(:,:)*kappabin(:,:),2) &
                 / sum((mbin2(:,:)+1.e-60_sp)/rhobin(:,:),2)
                 ! equation 7, petters and kreidenweis (2007)
          aw=(dw**3-dd**3)/(dw**3-dd**3*(1._sp-kappa)) ! from eq 6,p+k(acp,2007)
        case default
          print *,'error kappa_flag'
          stop
      end select
      ! koop et al. (2000) nucleation rate - due to homogeneous nucleation.
      jw(:)=koopnucrate(aw,t,p,sz2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! the number of ice crystals nucleated:
      dn01(:)=abs( npart(:)*(1._sp-exp(-jw(:)*mwat(:)/rhow*dt)) )

      if(t.gt.ttr) then
          dn01=0._sp
      endif
      !!!!
      ! total aerosol mass in each bin added together:
      dmaer01(:,:)=(mbin2_ice(:,:)*(spread(npartice(:),2,sz)+1.e-50_sp)+ &
                      mbin2(:,:)*spread(dn01(:),2,sz) ) 
      ! total water mass that will be in the ice bins:
      m01=(yice*npartice+mwat(:)*dn01(:)) 

      ! number conc. of liquid bins:
      npart(:)=npart(:)-dn01(:)
      ! number conc. of ice bins:
      npartice(:)=npartice(:)+dn01(:)
      ! new ice mass in bin:
      m01=m01/(npartice) 
      
      
      where(m01.gt.0._sp.and.npartice.gt.0._sp)
        yice=m01
      elsewhere
        yice=yice
      end where
      
      ! aerosol mass in ice bins
      mbin2_ice(:,:)=dmaer01(:,:)/(1.e-50_sp+spread(npartice,2,sz))

      ! latent heat of fusion:
      t=t+lf/cp*sum(mwat(:)*dn01(:))
    end subroutine icenucleation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! one time-step of the bin-microphysics                                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates one time-step of bin-microphysics
    subroutine bin_microphysics(func1,func2,func3)
    use nrtype
    use nr, only : zbrent
    implicit none
    real(sp) :: mass1, mass2, deltam, vapour_mass, liquid_mass, x1,x2 , cpm, &
        var, dummy
    integer(i4b) :: iloc
    
    interface
        subroutine func1(neq, tt, y, ydot, rpar, ipar)
            use nrtype
            use nr, only : dfridr,locate

            implicit none
            real(sp), intent(inout) :: tt
            real(sp), intent(inout), dimension(neq) :: y, ydot
            integer(i4b), intent(inout) :: neq
            real(sp), intent(inout) :: rpar
            integer(i4b), intent(inout) :: ipar
        end subroutine func1
    end interface
    interface
        subroutine func2(neq, tt, y, ydot, rpar, ipar)
            use nrtype
            use nr, only : dfridr,locate

            implicit none
            real(sp), intent(inout) :: tt
            real(sp), intent(inout), dimension(neq) :: y, ydot
            integer(i4b), intent(inout) :: neq
            real(sp), intent(inout) :: rpar
            integer(i4b), intent(inout) :: ipar
        end subroutine func2
    end interface
    interface
        subroutine func3(npart, npartice, mwat,mbin2,mbin2_ice, &
                         rhobin,nubin,kappabin,molwbin,t,p,sz,sz2,yice,rh,dt) 
            use nrtype
            implicit none
            real(sp), intent(inout) :: t
            real(sp), intent(in) :: p,rh,dt
            real(sp), dimension(sz2), intent(inout) :: npart,npartice
            real(sp), dimension(:), intent(in) :: mwat
            real(sp), dimension(:,:), intent(in) :: mbin2, &
                                                  rhobin,nubin,kappabin,molwbin
            integer(i4b), intent(in) :: sz,sz2
            real(sp), dimension(sz2,sz), intent(inout) :: mbin2_ice
            real(sp), intent(inout), dimension(sz2) :: yice
        end subroutine func3
    end interface
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mass balance                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(adiabatic_prof) then
        parcel1%yold=parcel1%y ! store old
        ! total water before:
        mass1=sum(parcel1%npart*parcel1%y(1:parcel1%n_bin_mode))+ &
            parcel1%y(parcel1%irh)*eps1* &
            svp_liq(parcel1%y(parcel1%ite)) / &
            (parcel1%y(parcel1%ipr)-svp_liq(parcel1%y(parcel1%ite)))
        if(ice_flag .eq. 1) then
            mass1=mass1+sum(parcel1%npartice*parcel1%yice(1:parcel1%n_bin_mode))
        endif
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
        ! check there are no negative values
        where(parcel1%y(1:parcel1%n_bin_mode).le.0.e1_sp)
            parcel1%y(1:parcel1%n_bin_mode)=1.e-22_sp
        end where
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    if(ice_flag .eq. 1) then
        ! ice part of the parcel model
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Ice nucleation                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call func3(parcel1%npart(1:parcel1%n_bin_mode), &
                parcel1%npartice(1:parcel1%n_bin_mode), &
                parcel1%y(1:parcel1%n_bin_mode), &
                parcel1%mbin(:,1:n_comps), &
                parcel1%mbinice(:,1:n_comps), &
                parcel1%rhobin(:,1:n_comps), &
                parcel1%nubin(:,1:n_comps), &
                parcel1%kappabin(:,1:n_comps), &
                parcel1%molwbin(:,1:n_comps), &
                parcel1%y(parcel1%ite), &
                parcel1%y(parcel1%ipr),&
                n_comps,parcel1%n_bin_mode, &
                parcel1%yice(1:parcel1%n_bin_mode), &
                parcel1%y(parcel1%irh), parcel1%dt) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ODE solver                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        parcel1%toutice=parcel1%ttice+parcel1%dt
        parcel1%yice(parcel1%ipri)=parcel1%y(parcel1%ipr)
        parcel1%yice(parcel1%itei)=parcel1%y(parcel1%ite)
        parcel1%yice(parcel1%irhi)=parcel1%y(parcel1%irh)
        do while (parcel1%ttice .lt. parcel1%toutice)
            parcel1%istateice=1
            call dvode(func2,parcel1%neqice,parcel1%yice,parcel1%ttice,parcel1%toutice,&
                           parcel1%itolice,parcel1%rtolice,parcel1%atolice,&
                           parcel1%itaskice,parcel1%istateice,parcel1%ioptice,&
                           parcel1%rworkice,parcel1%lrwice,&
                           parcel1%iworkice,parcel1%liwice,jparcelwarm, &
                           parcel1%mfice,parcel1%rparice,parcel1%iparice)
            ! check there are no negative values
            where(parcel1%yice(1:parcel1%n_bin_mode).le.0.e1_sp)
                parcel1%yice(1:parcel1%n_bin_mode)=1.e-22_sp
            end where
        enddo
        parcel1%y(parcel1%ipr)=parcel1%yice(parcel1%ipri)
        parcel1%y(parcel1%ite)=parcel1%yice(parcel1%itei)
        parcel1%y(parcel1%irh)=parcel1%yice(parcel1%irhi)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! stop the simulation if parcel is above cloud-top                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if((parcel1%y(parcel1%iz) .gt. parcel1%z_ctop)  .and. &
        vert_ent) parcel1%break_flag=.true.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! vertical entrainment outside of solver (see Sanchez et al, 2017, ACP)!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(vert_ent .and. (.not. adiabatic_prof) .and. &
        (parcel1%y(parcel1%iz) .gt. parcel1%z_cbase) ) then
        ! vapour mass:
        vapour_mass=parcel1%y(parcel1%irh)*eps1* &
            svp_liq(parcel1%y(parcel1%ite)) / &
            (parcel1%y(parcel1%ipr)-svp_liq(parcel1%y(parcel1%ite)))

        ! liquid mass:
        liquid_mass=sum(parcel1%npart*parcel1%y(1:parcel1%n_bin_mode))
        ! total water:
        mass2=liquid_mass+ vapour_mass

        if (set_theta_q_cb_flag) then
            parcel1%theta_q_cbase=calc_theta_q3(parcel1%y(parcel1%ite), &
                                    parcel1%y(parcel1%ipr), mass2)
            parcel1%q_cbase=mass2
            if(parcel1%y(parcel1%irh) .gt. 1._sp) then
                set_theta_q_cb_flag=.false.
            endif
        endif

        if(.not. set_theta_q_cb_flag) then
            cpm=cp+vapour_mass*cpv+liquid_mass*cpw
        
!             parcel1%theta_q=calc_theta_q3(parcel1%y(parcel1%ite), &
!                                         parcel1%y(parcel1%ipr), mass2)
            ! locate position
            iloc=locate(parcel1%z_sound(1:n_levels_s),parcel1%y(parcel1%iz))
            iloc=min(n_levels_s-1,iloc)
            iloc=max(1,iloc)
            ! linear interp theta_q
            call polint(parcel1%z_sound(iloc:iloc+1), parcel1%theta_q_sound(iloc:iloc+1), &
                    min(parcel1%y(parcel1%iz),parcel1%z_sound(n_levels_s)), var,dummy)        
            parcel1%theta_q=var
                                        
            ! equation 5 (Sanchez et al, 2017, acp)
            parcel1%x_ent=(parcel1%theta_q-parcel1%theta_q_cbase) / &
                    (parcel1%theta_q_ctop-parcel1%theta_q_cbase)
            x1=min(max(parcel1%x_ent,0._sp),1._sp)
            x2=1._sp-x1

            ! entrainment of vapour:
            vapour_mass= &!vapour_mass+&
                x1*(parcel1%q_ctop)+&
                x2*(parcel1%q_cbase)-liquid_mass
            ! entrainment of liquid:
            parcel1%npart=parcel1%npart*x2/x2old

            ! entrainment of theta_q (equation 4, Sanchez et al, 2017, ACP):
            parcel1%theta_q=x1*parcel1%theta_q_ctop + &
                x2*parcel1%theta_q_cbase
            
            p111=parcel1%y(parcel1%ipr)
            theta_q_sat=parcel1%theta_q
            parcel1%y(parcel1%ite)=zbrent(calc_theta_q,150._sp, &
                 theta_q_sat, 1.e-30_sp)
            print *,x1,x2

            ! rh
            parcel1%y(parcel1%irh)=vapour_mass / &
                 ( eps1*svp_liq(parcel1%y(parcel1%ite)) / &
                 (parcel1%y(parcel1%ipr)-svp_liq(parcel1%y(parcel1%ite))) )
            x2old=x2
        endif
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mass balance                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(adiabatic_prof) then
        ! total water after:
        mass2=sum(parcel1%npart*parcel1%y(1:parcel1%n_bin_mode))+ &
            parcel1%y(parcel1%irh)*eps1* &
            svp_liq(parcel1%y(parcel1%ite)) / &
            (parcel1%y(parcel1%ipr)-svp_liq(parcel1%y(parcel1%ite)))
        if(ice_flag .eq. 1) then
            mass2=mass2+sum(parcel1%npartice*parcel1%yice(1:parcel1%n_bin_mode))
        endif
        
        deltam=mass2-mass1
        vapour_mass=parcel1%y(parcel1%irh)*eps1* &
            svp_liq(parcel1%y(parcel1%ite)) / &
            (parcel1%y(parcel1%ipr)-svp_liq(parcel1%y(parcel1%ite)))
        ! adjust to conserve:
        vapour_mass=vapour_mass-deltam
        parcel1%y(parcel1%irh)=vapour_mass / &
            ( eps1*svp_liq(parcel1%y(parcel1%ite)) / &
            (parcel1%y(parcel1%ipr)-svp_liq(parcel1%y(parcel1%ite))) )
    endif    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    end subroutine bin_microphysics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! HELPER ROUTINE                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check(status)
    use netcdf
    use nrtype
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

    use nrtype
    use netcdf

    implicit none
    logical, intent(inout) :: new_file
    character (len=*),intent(in) :: outputfile
    ! output to netcdf file
    if(new_file) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! open / create the netcdf file                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call check( nf90_create(outputfile, NF90_CLOBBER, io1%ncid) )
        ! define dimensions (netcdf hands back a handle)
        call check( nf90_def_dim(io1%ncid, "times", NF90_UNLIMITED, io1%x_dimid) )
        call check( nf90_def_dim(io1%ncid, "nbins", n_bins, io1%bin_dimid) )
        call check( nf90_def_dim(io1%ncid, "nmodes", n_mode, io1%mode_dimid) )
        
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
        endif
        
        call check( nf90_enddef(io1%ncid) )
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
        sum(parcel1%y(1:parcel1%n_bin_mode)*parcel1%npart(1:parcel1%n_bin_mode)), &
                start = (/io1%icur/)))

    ! write variable: beta_ext
    call check( nf90_inq_varid(io1%ncid, "beta_ext", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        2._sp*sum((parcel1%y(1:parcel1%n_bin_mode)* &
            6._sp/(rhow*pi))**(2._sp/3._sp)* pi/4._sp* &
            parcel1%npart(1:parcel1%n_bin_mode)), &
                start = (/io1%icur/)))

    ! write variable: number > 2.5 microns (8.1812e-15 kg)
    parcel1%ndrop=0._sp
    where (parcel1%y(1:parcel1%n_bin_mode) > 6.5450e-14_sp)
        parcel1%ndrop=parcel1%npart(:)
    end where
    
    call check( nf90_inq_varid(io1%ncid, "ndrop", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        sum(parcel1%ndrop), start = (/io1%icur/)))

    ! write variable: effective radius
    call check( nf90_inq_varid(io1%ncid, "deff", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        sum((parcel1%y(1:parcel1%n_bin_mode)* &
            6._sp/(rhow*pi))**(3._sp/3._sp)*  &
            parcel1%npart(1:parcel1%n_bin_mode)) / &
        sum((parcel1%y(1:parcel1%n_bin_mode)* &
            6._sp/(rhow*pi))**(2._sp/3._sp)*  &
            parcel1%npart(1:parcel1%n_bin_mode)), &
                start = (/io1%icur/)))

    call check( nf90_inq_varid(io1%ncid, "mwat", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        reshape(parcel1%y(1:parcel1%n_bin_mode),(/n_bins,n_mode/)), start = (/1,1,io1%icur/)))

    if(ice_flag .eq. 1) then
        ! write variable: qi
        call check( nf90_inq_varid(io1%ncid, "qi", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            sum(parcel1%yice(1:parcel1%n_bin_mode)* &
                parcel1%npartice(1:parcel1%n_bin_mode)), &
                    start = (/io1%icur/)))

        ! write variable: number concentration of ice crystals
        parcel1%nice=parcel1%npartice
        call check( nf90_inq_varid(io1%ncid, "nice", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            sum(parcel1%nice), start = (/io1%icur/)))
    
    endif
    

    call check( nf90_close(io1%ncid) )


    io1%icur=io1%icur+1
    end subroutine output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! driver for bmm                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>driver for the bin-microphysics module
    subroutine bmm_driver()
    use nrtype
    implicit none
    integer(i4b) :: i, nt
    
    
    nt=ceiling(runtime / real(dt,kind=sp))
    do i=1,nt
        ! output to file
        call output(io1%new_file,outputfile)
        
        
        ! one time-step of model
        call bin_microphysics(fparcelwarm, fparcelcold, icenucleation)
        
        
             
        ! break-out if flag has been set 
        if(parcel1%break_flag) exit
    enddo
    ! output to file
    call output(io1%new_file,outputfile)
    
    
    
    
    end subroutine bmm_driver
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







	end module bmm	

