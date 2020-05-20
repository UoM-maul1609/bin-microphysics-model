	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>code to allocate arrays, and call activation 
	module sce
    use numerics_type
    use numerics, only : find_pos, poly_int
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>variables and types for the stochastic collection equation module

    implicit none
        ! constants for the sce module

        real(wp), parameter :: r_gas=8.314_wp, molw_a=29.e-3_wp,molw_water=18.e-3_wp, &
                                cp=1005.0_wp, cpv=1870._wp, cpw=4.27e3_wp, cpi=2104.6_wp, &
                                grav=9.81_wp, &
        						lv=2.5e6_wp, ls=2.837e6_wp, lf=ls-lv, ttr=273.15_wp, &
        						joules_in_an_erg=1.0e-7_wp,joules_in_a_cal=4.187e0_wp, &
        						rhow=1000._wp, ra=r_gas/molw_a,rv=r_gas/molw_water , &
        						eps1=ra/rv, rhoice=910._wp, qsmall=1.e-30_wp
        						

        type parcel
            ! variables for bin model
            integer(i4b) :: n_bins1,n_bins2,n_binst, &
                            n_modes,n_comps, n_bin_mode, n_bin_mode1, n_bin_mode2, &
                            n_bin_modea, n_bin_modea1, &
                            ice_flag
            real(wp) :: dt
            real(wp) :: z,p,t,w,rh
                        
                        
            ! liquid water
            real(wp), dimension(:), allocatable :: d, maer, npart, rho_core, &
                            rh_eq, rhoat, dw, da_dt, ndrop
            real(wp), dimension(:,:), allocatable :: mbin, rhobin, &
                            nubin,molwbin,kappabin, mbinedges,ecoll,ecoal, moments
                                ! all bins x all comps    
            integer(i4b), dimension(:), allocatable :: momenttype
            integer(i4b), dimension(:,:), allocatable :: indexc                            
            real(wp) :: tt, tout
            
            ! ice water
            real(wp), dimension(:), allocatable :: dice, maerice, npartice, rho_coreice, &
                            rh_eqice, rhoatice, dwice, da_dtice, nice, &
                            phi, rhoi, nump, rime
            real(wp), dimension(:,:), allocatable :: mbinice, rhobinice, &
                                        nubinice,molwbinice,kappabinice ! all bins x all comps                                

            
            
            logical :: break_flag=.false.
            
        end type parcel


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
        ! declare an io type
        type(io) :: io1
        

        ! some namelist variables
        logical :: micro_init=.true.
        real(wp) :: dmina,dmaxa
        real(wp) :: dminc,dmaxc
        
        real(wp) :: tinit,pinit,rhinit,alpha_therm, alpha_cond=1.0_wp, &
                    alpha_therm_ice=1._wp, alpha_dep=1._wp
        integer(i4b) :: microphysics_flag=0, kappa_flag,updraft_type, vent_flag=1, &
                        ice_flag=0, bin_scheme_flag=1
        real(wp) :: dz,dt, runtime
        integer(i4b) :: idum, n_sel
        real(wp) :: mult, rh_act

        ! aerosol setup
        integer(i4b) :: n_intern, n_mode,n_sv,sv_flag,n_bins,n_binsc,n_comps,kfac
        ! aerosol_spec
        real(wp), allocatable, dimension(:,:) :: n_aer1,d_aer1,sig_aer1, mass_frac_aer1
        real(wp), allocatable, dimension(:) ::  molw_core1,density_core1,nu_core1, &
                                        kappa_core1, ncloud

        ! cloud spec
        real(wp), allocatable, dimension(:,:) :: lwc, dbar
        
        

        ! variables for model
        real(wp) :: n_dummy, d_dummy


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
	!>@param[inout] n_aer1: number conc. in modes
	!>@param[inout] d_aer1: diameter in modes
	!>@param[inout] sig_aer1: geo std in modes
	!>@param[inout] mass_frac_aer1:mass_fraction of each component
	!>@param[inout] molw_core1:molw in core
	!>@param[inout] density_core1: solute density
	!>@param[inout] nu_core1: van hoff factor
	!>@param[inout] kappa_core1: kappa parameter
	!>@param[inout] lwc: liquid water content
	!>@param[inout] dbar: mean diameter
	!>@param[inout] ncloud: for initialisation
	subroutine allocate_arrays(n_intern,n_mode,n_sv,n_bins,n_comps, &
		                    n_aer1,d_aer1,sig_aer1,mass_frac_aer1, molw_core1, &
		                    density_core1, nu_core1, kappa_core1,lwc,dbar,ncloud)
		                    
		                    
		use numerics_type
		implicit none
		integer(i4b), intent(in) :: n_intern, n_mode, n_sv, n_bins,n_comps
		real(wp), dimension(:,:), allocatable, intent(inout) :: &
		                        n_aer1,d_aer1,sig_aer1,mass_frac_aer1,lwc,dbar
		real(wp), dimension(:), allocatable, intent(inout) :: molw_core1,density_core1, &
		                        nu_core1,kappa_core1,ncloud
		
		integer(i4b) :: AllocateStatus
		
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

		allocate( lwc(1:n_intern,1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( dbar(1:n_intern,1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
		allocate( ncloud(1:n_mode), STAT = AllocateStatus)
		if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
        lwc=0._wp
        dbar=10._wp

	end subroutine allocate_arrays
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! read in the namelist                                                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>read in the data from the namelists for the sce microphysics model
	!>@param[in] nmlfile
	subroutine read_in_sce_namelist(nmlfile)
		implicit none
        character (len=200), intent(in) :: nmlfile
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! namelists                                                            !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! define namelists for environment
        namelist /run_vars/ outputfile, runtime, dt, &
                    tinit,pinit,rhinit, &
                    microphysics_flag, ice_flag, bin_scheme_flag, &
                    kappa_flag
        namelist /aerosol_setup/ n_intern,n_mode,n_sv,sv_flag, n_bins,n_comps
        namelist /aerosol_spec/ n_aer1,d_aer1,sig_aer1, dmina,dmaxa, &
                                mass_frac_aer1, molw_core1, &
                                density_core1,nu_core1,kappa_core1
        namelist /cloud_setup/ n_binsc,kfac,dminc,dmaxc
        namelist /cloud_spec/ lwc,dbar
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists	and allocate arrays								   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
        read(8,nml=run_vars)
        read(8,nml=aerosol_setup)
        read(8,nml=cloud_setup)
        ! allocate memory / init
		call allocate_arrays(n_intern,n_mode,n_sv,n_bins,n_comps, &
		                    n_aer1,d_aer1,sig_aer1,mass_frac_aer1, molw_core1, &
		                    density_core1, nu_core1, kappa_core1, lwc, dbar,ncloud)
        
        read(8,nml=aerosol_spec)
        read(8,nml=cloud_spec)
        close(8)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine read_in_sce_namelist







	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! initialise arrays                                                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>interpolates the sounding to the grid
    subroutine initialise_sce_arrays()
    use numerics_type
    use numerics, only : find_pos, poly_int, zeroin, fmin,vode_integrate

    implicit none
    real(wp) :: num, ntot, number_per_bin, test, var1, &
                eps2, z1, z2, htry, hmin, var, dummy, mass, vol, rho
    real(wp), dimension(1) :: p1, z11
    real(wp) :: p11, p22, rm, cpm
    integer(i4b) :: i,j,k, AllocateStatus, iloc,i2, mode1,mode2,moden

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set variables and allocate arrays in parcel                                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    parcel1%n_bins1=n_bins 
    parcel1%n_bins2=n_binsc
    parcel1%n_binst=n_bins+n_binsc
    parcel1%n_modes=n_mode
    parcel1%n_comps=n_comps
    
    parcel1%n_bin_modea=parcel1%n_bins1*n_mode
    parcel1%n_bin_modea1=(parcel1%n_bins1+1)*n_mode ! extra due to bin edges
    parcel1%n_bin_mode=parcel1%n_binst*n_mode
    parcel1%n_bin_mode1=(parcel1%n_binst+1)*n_mode ! extra due bin edges
    parcel1%n_bin_mode2=(parcel1%n_binst+2)*n_mode ! extra due bin edges
    
    parcel1%p=pinit
    parcel1%t=tinit
    parcel1%rh=rhinit
    parcel1%dt=dt
    allocate( parcel1%d(1:parcel1%n_bin_mode1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%mbinedges(1:parcel1%n_binst+1,1:parcel1%n_modes), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%maer(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%npart(1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%mbin(1:parcel1%n_bin_mode,1:n_comps+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%rho_core(1:parcel1%n_modes), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	

    allocate( parcel1%moments(1:parcel1%n_bin_mode,1:n_comps), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%momenttype(1:n_comps), STAT = AllocateStatus)
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

    allocate( parcel1%ecoal(1:parcel1%n_bin_mode,1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%ecoll(1:parcel1%n_bin_mode,1:parcel1%n_bin_mode), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory ***"	
    allocate( parcel1%indexc(1:parcel1%n_bin_mode,1:parcel1%n_bin_mode), STAT = AllocateStatus)
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
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Aerosol                                                                  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        idum=k ! this is sent through to zbrent to select the correct mode
        ! find total number in mode between dmina and dmaxa:
        call lognormal_n_between_limits(n_aer1(:,k),d_aer1(:,k),sig_aer1(:,k), &
                                        n_intern,dmina,dmaxa, num)
        !print *,num
        ! set up variables for parcel model
        ntot=num
        number_per_bin=ntot/real(parcel1%n_bins1,wp)
        parcel1%npart(1+(k-1)*parcel1%n_binst:(k-1)*parcel1%n_binst+parcel1%n_bins1)=&
            number_per_bin
        parcel1%d(1+(k-1)*(parcel1%n_binst+1))=dmina
        do i=1,parcel1%n_bins1
            d_dummy=parcel1%d(i+(k-1)*(parcel1%n_binst+1))
            n_dummy=number_per_bin*(1._wp-1.e-15_wp)
            parcel1%d(i+1+(k-1)*(parcel1%n_binst+1))= zeroin(&
                        d_dummy*0.9_wp,dmaxa*2._wp,find_upper_diameter, 1.e-30_wp)
        enddo
        parcel1%d((parcel1%n_bins1+1)+ &
            (k-1)*(parcel1%n_binst+1))=dmaxa ! nail it to end point - round off
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! aerosol mass - total                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=1,parcel1%n_modes
        do j=1,parcel1%n_binst
            i=j+(k-1)*(parcel1%n_binst+1)
            if (j.le.parcel1%n_bins1) then
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! from aerosol distribution                                        !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                parcel1%maer(j+(k-1)*(parcel1%n_binst))= &
                    pi/6._wp*(0.5_wp*(parcel1%d(i+1)+parcel1%d(i)))**3 * &
                    parcel1%rho_core(k)
            else
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! just put some aerosol in cloud distribution                      !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                parcel1%maer(j+(k-1)*(parcel1%n_binst))= &
                    pi/6._wp*(1.e-7_wp)**3 * &
                    parcel1%rho_core(k)            
            endif
        enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the mass of each component in a bin, including water (Koehler eq)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,parcel1%n_binst
        do j=1,parcel1%n_modes
            do k=1,parcel1%n_comps
                parcel1%mbin(i+(j-1)*parcel1%n_binst,k)= &  
                        parcel1%maer(i+(j-1)*parcel1%n_binst)*mass_frac_aer1(j,k)
                ! density in each bin:
                parcel1%rhobin(i+(j-1)*parcel1%n_binst,k)=density_core1(k)
                ! nu in each bin:
                parcel1%nubin(i+(j-1)*parcel1%n_binst,k)=nu_core1(k)
                ! molw in each bin:
                parcel1%molwbin(i+(j-1)*parcel1%n_binst,k)=molw_core1(k)
                ! kappa in each bin:
                parcel1%kappabin(i+(j-1)*parcel1%n_binst,k)=kappa_core1(k)
            enddo
        enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! put water on bin, using koehler equation                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    select case(kappa_flag)
        case(0)
            do i=1,parcel1%n_bin_mode
                if ((modulo(i,parcel1%n_binst).le.parcel1%n_bins1).and.&
                    (modulo(i,parcel1%n_binst).gt.0)) then
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! Aerosol water                                                   !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    n_sel=i
                    rh_act=0._wp !min(parcel1%rh,0.999_wp)
                    mult=-1._wp
                    ! has to be less than the peak moles of water at activation
                    test=fmin(1.e-50_wp,1.e1_wp, koehler02,1.e-30_wp)
                    rh_act=min(parcel1%rh,0.999_wp)
                    mult=1._wp
                    d_dummy=zeroin(1.e-30_wp, test, koehler02,1.e-30_wp)*molw_water 
                    parcel1%mbin(i,n_comps+1)= d_dummy
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                else
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! cloud water                                                     !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! water mass bin centres
                    parcel1%mbin(i,n_comps+1)= &
                        2._wp**(1._wp/kfac)*parcel1%mbin(i-1,n_comps+1)                    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                endif
            enddo
        case(1)
            do i=1,parcel1%n_bin_mode
                if ((modulo(i,parcel1%n_binst).le.parcel1%n_bins1).and.&
                    (modulo(i,parcel1%n_binst).gt.0)) then
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! Aerosol water                                                   !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    n_sel=i
                    rh_act=0._wp !min(parcel1%rh,0.999_wp)
                    mult=-1._wp
                    ! has to be less than the peak moles of water at activation
                    test=fmin(1.e-50_wp,1.e1_wp, kkoehler02,1.e-30_wp)
                    rh_act=min(parcel1%rh,0.999_wp)
                    mult=1._wp
                    d_dummy=zeroin(1.e-30_wp, test, kkoehler02,1.e-30_wp)*molw_water 
                    parcel1%mbin(i,n_comps+1)= d_dummy
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                else
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! cloud water                                                     !
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! water mass bin centres
                    parcel1%mbin(i,n_comps+1)= &
                        2._wp**(1._wp/kfac)*parcel1%mbin(i-1,n_comps+1)                    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                endif
            enddo
            
        case default
            print *,'error kappa flag'
            stop
    end select
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set the water bin-edges                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=1,parcel1%n_modes
        parcel1%mbinedges(1,k)=0._wp
        parcel1%mbinedges(parcel1%n_binst+1,k)=pi/6._wp*rhow*dmaxc**3
        do j=2,parcel1%n_binst
            parcel1%mbinedges(j,k)= &
                10**(0.5_wp*(log10(parcel1%mbin(j,n_comps+1))+&
                     log10(parcel1%mbin(j-1,n_comps+1))))
        enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise the cloud PSD                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=1,parcel1%n_modes
        ! above the aerosol
        do j=1+parcel1%n_bins1,parcel1%n_binst
            ncloud=-1._wp/(4._wp/3._wp*pi*rhow*dbar(:,k)**3)*lwc(:,k)* &
                (exp(-parcel1%mbinedges(j+1,k)/(4._wp/3._wp*pi*rhow*dbar(:,k)**3))-&
                exp(-parcel1%mbinedges(j,k)/(4._wp/3._wp*pi*rhow*dbar(:,k)**3)))
            ncloud=max(ncloud,0._wp)


            parcel1%npart(j+(k-1)*parcel1%n_binst) = &
                sum(ncloud)
        enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise conserved moments                                                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j=1,parcel1%n_comps
        ! above the aerosol
        do i=1,parcel1%n_bin_mode
            ! aerosol moments
            parcel1%moments(i,j)=parcel1%npart(i)*parcel1%mbin(i,j)
        enddo
    enddo
    parcel1%momenttype(1:parcel1%n_comps)=1 ! 1 is mass, 2 is number
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the diameters                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call wetdiam(parcel1%mbin(:,n_comps+1),parcel1%mbin,&
        parcel1%rhobin,parcel1%n_bin_mode,parcel1%dw) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! collision efficiency - gravitational settling                                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call ecollision(parcel1%dw,parcel1%n_bin_mode,parcel1%ecoll) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! which mode should gain go into?                                              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    moden=parcel1%n_bin_mode/parcel1%n_binst
    do j=1,parcel1%n_bin_mode
        do i=1,parcel1%n_bin_mode
            ! which mode is first in?
            mode1=(i-1)/parcel1%n_binst+1
            ! which mode is second in?
            mode2=(j-1)/parcel1%n_binst+1
            ! if both in same mode, put them in same mode
            ! if different, put in last mode
            if (mode1.eq.mode2) then
                parcel1%indexc(i,j)=mode1
            else
                parcel1%indexc(i,j)=moden
            endif
        enddo
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end subroutine initialise_sce_arrays
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! collision efficiency                                                         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the collision efficiency between two droplets according to the 
	!> `Long' kernel
	!>@param[in] dw - diameters of droplets
	!>@param[in] sz - size of array
	!>@param[inout] ecoll: collision efficiency
    subroutine ecollision(dw,sz,ecoll) 
    integer(i4b), intent(in) :: sz
    real(wp), intent(in), dimension(sz) :: dw
    real(wp), intent(inout), dimension(sz,sz) :: ecoll
    
    integer(i4b) :: i,j
    real(wp) :: x1,x2,as,al,r,u,v
       
    do j=1,sz
        do i=1,sz ! left-most should vary quickest
            x1=min(dw(i),dw(j))
            x2=max(dw(i),dw(j))
            as=x1/2._wp
            al=x2/2._wp
            
            r = max(x2*1.e6_wp,x1*1.e6_wp) ! microns
            u = 4._wp/3._wp*pi*al**3*1.e6_wp    ! cm^3
            v = 4._wp/3._wp*pi*as**3*1.e6_wp    ! cm^3
            if(r<=50._wp) then
                ecoll(i,j) = 9.44e9_wp*(u**2+v**2)/1.e6_wp
            else
                ecoll(i,j) = 5.78e3_wp*(u+v)/1.e6_wp
            endif
        enddo
    enddo
    
    end subroutine ecollision
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
    ! one time-step of the bin-microphysics                                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates one time-step of sce-microphysics
    subroutine sce_microphysics(n_binst,n_bin_mode,n_moments,npart,moments,momtype, &
                                ecoll,indexc,xn)
    use numerics_type
    use numerics, only : zeroin, dvode
    implicit none
    integer(i4b), intent(in) :: n_binst,n_bin_mode, n_moments
    real(wp), dimension(n_bin_mode,n_bin_mode), intent(in) :: ecoll
    integer(i4b), dimension(n_bin_mode,n_bin_mode), intent(in) :: indexc
    real(wp), dimension(n_bin_mode), intent(inout) :: npart,xn
    real(wp), dimension(n_bin_mode,n_moments), intent(inout) :: moments
    integer(i4b), dimension(n_moments), intent(in) :: momtype
    
    real(wp) :: remove1,remove2,massn,massaddto,nnew,gk,beta1,cw,fk05, &
                frac1, frac2, fracl, fracadj1, fracadj2, totloss
    real(wp), dimension(n_moments) :: momtemp
    integer(i4b) :: i,j,k,l,il,ih,jl,jh, modeinto
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Find the min and max bins to do computations on                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    il=n_bin_mode
    ih=1
    do i=1,n_bin_mode
        if (npart(i).gt.qsmall) ih=i            
    enddo
    do i=n_bin_mode,1,-1
        if (npart(i).gt.qsmall) il=i            
    enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Main loop                                                                          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=il,ih
        if (npart(i).lt.qsmall) cycle
        do j=i+1,ih
            ! numbers removed from each bin:
            remove1=min(npart(i)*ecoll(j,i)*npart(j),npart(i))
            remove2=min(npart(i)*ecoll(j,i)*npart(j),npart(j))
            
            ! interaction creates a drops of this mass:
            massn = xn(i) + xn(j)
            ! total mass removed (i.e. added to new bin):
            massaddto=xn(i)*remove1+xn(j)*remove2
            ! number to add to new bin:
            nnew = massaddto/massn
            
            if (nnew.lt.qsmall) cycle
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! now, which bin does the new particle go into?
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            modeinto=indexc(j,i)
            jl=(modeinto-1)*n_binst+1
            jh=(modeinto)*n_binst
            do k=jl,jh
                if (xn(k).gt.massn) exit
            enddo
            l=k-1
            if (l.eq.jh) cycle
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



            !*****
            ! fraction to be removed from both bins
            frac1=remove1/npart(i)
            frac2=remove2/npart(j)
            !*****
            
            
            ! remove the particles from the colliding bins
            npart(i)=npart(i)-remove1
            npart(j)=npart(j)-remove2  
            
            
            
            !*****
            ! the moments that need to be transferred out
            momtemp=moments(i,:)*frac1+moments(j,:)*frac2
            ! remove the moments from the colliding bins
            moments(i,:)=moments(i,:)*(1._wp-frac1)
            moments(j,:)=moments(j,:)*(1._wp-frac2)
            !*****
            
            

            ! add the mass into the new bin:
            gk=npart(l)*xn(l)+massaddto
            
            ! now for the flux:
            beta1=log(npart(l+1)*xn(l+1)/gk+qsmall)
            ! courant number
            cw=(log(massn)-log(xn(l))) / (log(xn(l+1))-log(xn(l)))
            ! exponential flux
            fk05=massaddto/beta1*(exp(beta1*0.5_wp)-exp(beta1*(0.5_wp-cw)))
            fk05=min(fk05,massaddto)
            
            
            ! now apply the flux:
            fracadj1=(massaddto-fk05)/xn(l)
            fracadj2=(fk05)/xn(l+1)
            npart(l)=npart(l)+fracadj1
            npart(l+1)=npart(l+1)+fracadj2


            
            !*****
            ! the partitioning between bin l and l+1 is by mass fraction
            ! for mass variables. Fraction of total going into l:
            totloss=remove1+remove2
            fracadj1=fracadj1/(totloss)
            fracadj2=fracadj2/(totloss)
            fracl=(massaddto-fk05)/massaddto
            ! add the 'loss' moments to the new bin
            do k=1,n_moments
                if (momtemp(k).eq.2) then
                    moments(l,k)=moments(l,k)+momtemp(k)*fracadj1
                    moments(l+1,k)=moments(l+1,k)+momtemp(k)*fracadj2
                elseif(momtemp(k).eq.1) then
                    moments(l,k)=moments(l,k)+momtemp(k)*fracl
                    moments(l+1,k)=moments(l+1,k)+momtemp(k)*(1._wp-fracl)
                endif
            enddo
            !*****
        enddo
    enddo    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    



    
    end subroutine sce_microphysics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    







    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! driver for bmm                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>driver for the sce module
    subroutine sce_driver()
    use numerics_type
    implicit none
    integer(i4b) :: i, nt,j
    
    
    nt=ceiling(runtime / real(dt,kind=wp))
    do i=1,nt
        ! output to file
        call output(io1%new_file,outputfile)
        
        
        ! one time-step of model
        call sce_microphysics(parcel1%n_binst,parcel1%n_bin_mode,parcel1%n_comps,&
                            parcel1%npart,parcel1%moments,parcel1%momenttype, &
                            parcel1%ecoll,parcel1%indexc, &
                            parcel1%mbin(:,n_comps+1))
        
        
        ! redefine the mass of each component of aerosol
        do j=1,parcel1%n_comps
            where (parcel1%npart(:).gt.qsmall)
                parcel1%mbin(:,j)=parcel1%moments(:,j)/parcel1%npart(:)
            end where
        enddo    
             
        ! break-out if flag has been set 
        if(parcel1%break_flag) exit
    enddo
    ! output to file
    call output(io1%new_file,outputfile)
    
    
    
    end subroutine sce_driver
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
    ! output to netcdf file
    if(new_file) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! open / create the netcdf file                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call check( nf90_create(outputfile, NF90_CLOBBER, io1%ncid) )
        ! define dimensions (netcdf hands back a handle)
        call check( nf90_def_dim(io1%ncid, "times", NF90_UNLIMITED, io1%x_dimid) )
        call check( nf90_def_dim(io1%ncid, "nbins", n_bins, io1%bin_dimid) )
        call check( nf90_def_dim(io1%ncid, "nbinst", parcel1%n_binst, io1%bin2_dimid) )
        call check( nf90_def_dim(io1%ncid, "nbinsedge", parcel1%n_binst+1, io1%bin3_dimid) )
        call check( nf90_def_dim(io1%ncid, "nmodes", n_mode, io1%mode_dimid) )
        call check( nf90_def_dim(io1%ncid, "ncomps", n_comps, io1%comp_dimid) )
        
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
                   
        ! define variable: mwatedge
        call check( nf90_def_var(io1%ncid, "mwatedge", NF90_DOUBLE, &
                    (/io1%bin3_dimid,io1%mode_dimid/), io1%varid) )
        ! get id to a_dimid
        call check( nf90_inq_varid(io1%ncid, "mwatedge", io1%a_dimid) )
        ! units
        call check( nf90_put_att(io1%ncid, io1%a_dimid, &
                   "units", "kg") )

        ! define variable: mwat
        call check( nf90_def_var(io1%ncid, "mwat", NF90_DOUBLE, &
                    (/io1%bin2_dimid,io1%mode_dimid, io1%x_dimid/), io1%varid) )
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
!     ! write variable: z
!     call check( nf90_inq_varid(io1%ncid, "z", io1%varid ) )
!     call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%iz), &
!                 start = (/io1%icur/)))
! 
!     ! write variable: p
!     call check( nf90_inq_varid(io1%ncid, "p", io1%varid ) )
!     call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%ipr), &
!                 start = (/io1%icur/)))
! 
!     ! write variable: t
!     call check( nf90_inq_varid(io1%ncid, "t", io1%varid ) )
!     call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%ite), &
!                 start = (/io1%icur/)))
! 
!     ! write variable: rh
!     call check( nf90_inq_varid(io1%ncid, "rh", io1%varid ) )
!     call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%irh), &
!                 start = (/io1%icur/)))
! 
!     ! write variable: w
!     call check( nf90_inq_varid(io1%ncid, "w", io1%varid ) )
!     call check( nf90_put_var(io1%ncid, io1%varid, parcel1%y(parcel1%iw), &
!                 start = (/io1%icur/)))
! 
! 
!     ! write variable: ql
!     call check( nf90_inq_varid(io1%ncid, "ql", io1%varid ) )
!     call check( nf90_put_var(io1%ncid, io1%varid, &
!         sum(parcel1%y(1:parcel1%n_bin_mode)*parcel1%npart(1:parcel1%n_bin_mode)), &
!                 start = (/io1%icur/)))
! 
!     ! write variable: beta_ext
!     call check( nf90_inq_varid(io1%ncid, "beta_ext", io1%varid ) )
!     call check( nf90_put_var(io1%ncid, io1%varid, &
!         2._wp*sum((parcel1%y(1:parcel1%n_bin_mode)* &
!             6._wp/(rhow*pi))**(2._wp/3._wp)* pi/4._wp* &
!             parcel1%npart(1:parcel1%n_bin_mode)), &
!                 start = (/io1%icur/)))
! 
    ! write variable: number > 2.5 microns (8.1812e-15 kg)
    parcel1%ndrop=0._wp
    where (parcel1%mbin(1:parcel1%n_bin_mode,parcel1%n_comps+1) > 6.5450e-14_wp)
        parcel1%ndrop=parcel1%npart(:)
    end where
    
    call check( nf90_inq_varid(io1%ncid, "ndrop", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        sum(parcel1%ndrop), start = (/io1%icur/)))
! 
!     ! write variable: effective radius
!     call check( nf90_inq_varid(io1%ncid, "deff", io1%varid ) )
!     call check( nf90_put_var(io1%ncid, io1%varid, &
!         sum((parcel1%y(1:parcel1%n_bin_mode)* &
!             6._wp/(rhow*pi))**(3._wp/3._wp)*  &
!             parcel1%npart(1:parcel1%n_bin_mode)) / &
!         sum((parcel1%y(1:parcel1%n_bin_mode)* &
!             6._wp/(rhow*pi))**(2._wp/3._wp)*  &
!             parcel1%npart(1:parcel1%n_bin_mode)), &
!                 start = (/io1%icur/)))
! 
    if (io1%icur==1) then
        call check( nf90_inq_varid(io1%ncid, "mwatedge", io1%varid ) )
        call check( nf90_put_var(io1%ncid, io1%varid, &
            reshape(parcel1%mbinedges(1:parcel1%n_binst+1,1:parcel1%n_modes), &
            (/parcel1%n_binst+1,n_mode/)), start = (/1,1/)))
    endif

    call check( nf90_inq_varid(io1%ncid, "mwat", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        reshape(parcel1%mbin(1:parcel1%n_bin_mode,parcel1%n_comps+1), &
        (/parcel1%n_binst,n_mode/)), start = (/1,1,io1%icur/)))

    call check( nf90_inq_varid(io1%ncid, "nwat", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        reshape(parcel1%npart(1:parcel1%n_bin_mode), &
        (/parcel1%n_binst,n_mode/)), start = (/1,1,io1%icur/)))

    call check( nf90_inq_varid(io1%ncid, "maer", io1%varid ) )
    call check( nf90_put_var(io1%ncid, io1%varid, &
        reshape(parcel1%mbin(1:parcel1%n_bin_mode,1:parcel1%n_comps), &
        (/parcel1%n_binst,n_mode,n_comps/)), start = (/1,1,1,io1%icur/)))

! 
!     if(ice_flag .eq. 1) then
!         ! write variable: qi
!         call check( nf90_inq_varid(io1%ncid, "qi", io1%varid ) )
!         call check( nf90_put_var(io1%ncid, io1%varid, &
!             sum(parcel1%yice(1:parcel1%n_bin_mode)* &
!                 parcel1%npartice(1:parcel1%n_bin_mode)), &
!                     start = (/io1%icur/)))
! 
!         ! write variable: number concentration of ice crystals
!         parcel1%nice=parcel1%npartice
!         call check( nf90_inq_varid(io1%ncid, "nice", io1%varid ) )
!         call check( nf90_put_var(io1%ncid, io1%varid, &
!             sum(parcel1%nice), start = (/io1%icur/)))
!     
!         call check( nf90_inq_varid(io1%ncid, "mice", io1%varid ) )
!         call check( nf90_put_var(io1%ncid, io1%varid, &
!             reshape(parcel1%yice(1:parcel1%n_bin_mode),(/n_bins,n_mode/)), &
!              start = (/1,1,io1%icur/)))
! 
!     endif
!     

    call check( nf90_close(io1%ncid) )


    io1%icur=io1%icur+1
    end subroutine output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	end module sce	

