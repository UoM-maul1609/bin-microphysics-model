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
                        runtime, dt, zinit, tpert, use_prof_for_tprh, chamber_override, &
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
                        BIN_FULL_MOVING, BIN_MOVING_CENTRE, BIN_CHEN_LAMB
                        
        use sce, only : read_in_sce_namelist, initialise_sce_arrays, &
                        n_binsc, n_binst, &
                        kfac, dminc, dmaxc, lwc, dbar, iwc, dbari, parcel1

        implicit none
		logical :: fixed_grid_required
        character (len=200) :: nmlfile = ' '




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists													   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        call read_in_bmm_namelist(nmlfile)
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
            call initialise_sce_arrays(n_bins, n_binsc,n_mode, n_comps, n_intern, &
                    ice_flag, &
                    pinit,tinit,rhinit,dt,dmina,dmaxa,dminc,dmaxc,&
                    mass_frac_aer1,density_core1,nu_core1,molw_core1, kappa_core1, &
                    n_aer1,d_aer1,sig_aer1)
            n_bins=n_binst ! n_bins from parcel is increased so allocations are the 
                           ! same as SCE code
        endif 
        ! initialise parcel1 arrays in bmm module       
        call initialise_bmm_arrays(psurf, tsurf, q_read, theta_read, rh_read, z_read, &
        			time_chamber, press_chamber, temp_chamber, qtot_chamber, &
                    runtime, dt, zinit, tpert, use_prof_for_tprh, chamber_override, &
                    winit, tinit, pinit, &
                    rhinit, radinit, bubble_flag, &
                    microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                    kappa_flag, updraft_type, adiabatic_prof, vert_ent, z_ctop, &
                    ent_rate, n_levels_s, n_levels_c, &
                    alpha_therm, alpha_cond, alpha_therm_ice, &
                    alpha_dep, n_intern, n_mode, n_sv, sv_flag, n_bins, n_comps, &
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



