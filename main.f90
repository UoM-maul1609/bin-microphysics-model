	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2018
	!>@brief
	!>Bin Microphysics Model (BMM): 
	!>Bin cloud parcel model based on earlier ACPIM
	!> <br><br>
	!>\f$ F\left(t,z \right)
	!>   = initialisation,microphysics,etc \f$
	!> <br><br>
	!> compile using the Makefile (note requires netcdf) and then run using: <br>
	!> ./main.exe namelist.in
	!> <br><br>
	!> (namelist used for initialisation).



	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme reads in information, allocates arrays, then calls the model driver

    program main
        use numerics_type
        use bmm, only : read_in_bmm_namelist, initialise_bmm_arrays, bmm_driver, io1, &
                        write_sce_to_bmm, write_sce_grid_to_bmm, &
                        project_initial_bmm_to_fixed_grid, &
                        scefile, sce_flag, hm_flag, break_flag, mode1_flag,mode2_flag, &
                        psurf, tsurf, q_read, theta_read, rh_read, z_read, &
                        time_chamber, press_chamber, temp_chamber, qtot_chamber, &
                        runtime, dt, zinit, tpert, use_prof_for_tprh, &
                        winit, &
                        winit2, amplitude2, tinit, pinit, &
                        rhinit, radinit, bubble_flag, &
                        microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                        kappa_flag, updraft_type, t_thresh, &
                        adiabatic_prof, entrain_period, thresh_to_start_hom_mix, &
                        vert_ent, z_ctop, ent_rate, &
                        n_levels_s, n_levels_c, alpha_therm, alpha_cond, alpha_therm_ice, &
                        alpha_dep, n_intern, n_mode, n_sv, sv_flag, n_bins, n_comps, &
                        n_aer1,d_aer1,sig_aer1,dmina,dmaxa,mass_frac_aer1,molw_core1, &
                        density_core1, nu_core1, kappa_core1, org_content1, molw_org1, &
                        kappa_org1, density_org1, delta_h_vap1,nu_org1, log_c_star1, &
                        BIN_FULL_MOVING, BIN_MOVING_CENTRE, BIN_CHEN_LAMB, &
                        FIXED_GRID_LEGACY, FIXED_GRID_HYBRID, fixed_grid_mode
                        
        use sce, only : read_in_sce_namelist, initialise_sce_arrays, &
                        n_binsc, n_binst, &
                        kfac, dminc, dmaxc, lwc, dbar, iwc, dbari, parcel1

        implicit none
		logical :: fixed_grid_required
        integer(i4b) :: n_bins_aerosol, n_bins_init, sce_grid_mode
        character (len=200) :: nmlfile = ' '




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists													   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        call read_in_bmm_namelist(nmlfile)
        n_bins_aerosol=n_bins

        ! check if fixed_grid required
		fixed_grid_required = (sce_flag.gt.0) .or. &
			(bin_scheme_flag.eq.BIN_MOVING_CENTRE) .or. &
			(bin_scheme_flag.eq.BIN_CHEN_LAMB)
        if (fixed_grid_required) then
            call read_in_sce_namelist(scefile,nmlfile)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! allocate and initialise the grid                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (fixed_grid_required) then
            ! note, this initialises parcel1 arrays in sce module
			! For bin schemes 1 and 2 this may be used only to construct
			! the common fixed mass grid; it does not imply SCE is run.
            ! Full-moving also uses the hybrid construction for the appended
            ! zero-number receiving categories: the first n_bins aerosol cohorts
            ! retain their equal-number initialisation, while bins n_bins+1 onward
            ! continue geometrically in water mass to dmaxc.  The collision gain
            ! treatment remains full-moving; this only supplies sensible empty-bin
            ! reference pivots.
            if (bin_scheme_flag.eq.BIN_FULL_MOVING) then
                sce_grid_mode=FIXED_GRID_HYBRID
            else
                sce_grid_mode=fixed_grid_mode
            endif

            call initialise_sce_arrays(n_bins_aerosol, n_binsc,n_mode, n_comps, n_intern, &
                    ice_flag, &
                    pinit,tinit,rhinit,dt,dmina,dmaxa,dminc,dmaxc,&
                    mass_frac_aer1,density_core1,nu_core1,molw_core1, kappa_core1, &
                    n_aer1,d_aer1,sig_aer1, sce_grid_mode, kappa_flag)

            ! BMM arrays must span the complete fixed/SCE grid.
            !
            ! Full-moving always starts from the aerosol resolution requested
            ! in the BMM namelist.  Thus, for e.g. n_bins=60 and n_binsc=80,
            ! bins 1:60 per mode contain the 60 equal-number aerosol cohorts
            ! and bins 61:140 are initially zero-number receiving categories.
            ! This makes the initial aerosol discretisation independent of
            ! whether SCE is enabled.
            !
            ! Hybrid moving-centre/Chen-Lamb uses the same n_bins aerosol
            ! populations followed by empty geometric cloud bins.  Legacy
            ! fixed-grid mode retains its historical n_binst-population
            ! initialisation for reproducibility.
            n_bins=n_binst
            if (bin_scheme_flag.eq.BIN_FULL_MOVING) then
                n_bins_init=n_bins_aerosol
            elseif (((bin_scheme_flag.eq.BIN_MOVING_CENTRE) .or. &
                     (bin_scheme_flag.eq.BIN_CHEN_LAMB)) .and. &
                    (fixed_grid_mode.eq.FIXED_GRID_HYBRID)) then
                n_bins_init=n_bins_aerosol
            else
                n_bins_init=n_bins
            endif
        else
            n_bins_init=n_bins
        endif 
        ! initialise parcel1 arrays in bmm module       
        call initialise_bmm_arrays(psurf, tsurf, q_read, theta_read, rh_read, z_read, &
        			time_chamber, press_chamber, temp_chamber, qtot_chamber, &
                    runtime, dt, zinit, tpert, use_prof_for_tprh, &
                    winit, tinit, pinit, &
                    rhinit, radinit, bubble_flag, &
                    microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                    kappa_flag, updraft_type, adiabatic_prof, vert_ent, z_ctop, &
                    ent_rate, n_levels_s, n_levels_c, &
                    alpha_therm, alpha_cond, alpha_therm_ice, &
                    alpha_dep, n_intern, n_mode, n_sv, sv_flag, n_bins, n_bins_init, n_comps, &
                    n_aer1,d_aer1,sig_aer1,dmina,dmaxa,mass_frac_aer1,molw_core1, &
                    density_core1, nu_core1, kappa_core1, org_content1, molw_org1, &
                    kappa_org1, density_org1, delta_h_vap1,nu_org1, log_c_star1,sce_flag)
        
        ! This code writes the SCE variables to the BMM            
        if(fixed_grid_required) then
            ! send the SCE arrays, and use the local BMM arrays to map
            ! parcel1 here are the sce module vars. They are written to the bmm version of
            ! parcel1
			call write_sce_grid_to_bmm( &
				parcel1%n_bin_mode,parcel1%n_bin_modew,parcel1%n_binst, &
				parcel1%n_modes,parcel1%n_comps, &
				parcel1%indexc,parcel1%mbinedges)
			if (bin_scheme_flag.ne.BIN_FULL_MOVING) then
				call project_initial_bmm_to_fixed_grid()
			endif
        endif        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! run the model                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        io1%new_file=.true.
        call bmm_driver(sce_flag,hm_flag,break_flag,mode1_flag,mode2_flag) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main



