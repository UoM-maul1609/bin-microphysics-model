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
                        write_sce_to_bmm, &
                        scefile, sce_flag, hm_flag, break_flag, mode1_flag,mode2_flag, &
                        psurf, tsurf, q_read, theta_read, rh_read, z_read, &
                        runtime, dt, zinit, tpert, use_prof_for_tprh, winit, &
                        winit2, amplitude2, tinit, pinit, &
                        rhinit, radinit, bubble_flag, &
                        microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                        kappa_flag, updraft_type, t_thresh, &
                        adiabatic_prof, entrain_period, thresh_to_start_hom_mix, &
                        vert_ent, z_ctop, &
                        ent_rate, n_levels_s, alpha_therm, alpha_cond, alpha_therm_ice, &
                        alpha_dep, n_intern, n_mode, n_sv, sv_flag, n_bins, n_comps, &
                        n_aer1,d_aer1,sig_aer1,dmina,dmaxa,mass_frac_aer1,molw_core1, &
                        density_core1, nu_core1, kappa_core1, org_content1, molw_org1, &
                        kappa_org1, density_org1, delta_h_vap1,nu_org1, log_c_star1
                        
        use sce, only : read_in_sce_namelist, initialise_sce_arrays, &
                        n_binsc, n_binst, &
                        kfac, dminc, dmaxc, lwc, dbar, iwc, dbari, parcel1

        implicit none

        character (len=200) :: nmlfile = ' '




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! read in namelists													   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getarg(1,nmlfile)
        call read_in_bmm_namelist(nmlfile)
        if (sce_flag.gt.0) then
            call read_in_sce_namelist(scefile)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! allocate and initialise the grid                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (sce_flag.gt.0) then
            ! note, this initialises parcel1 arrays in sce module
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
                    runtime, dt, zinit, tpert, use_prof_for_tprh, winit, tinit, pinit, &
                    rhinit, radinit, bubble_flag, &
                    microphysics_flag, ice_flag, bin_scheme_flag, vent_flag, &
                    kappa_flag, updraft_type, adiabatic_prof, vert_ent, z_ctop, &
                    ent_rate, n_levels_s, alpha_therm, alpha_cond, alpha_therm_ice, &
                    alpha_dep, n_intern, n_mode, n_sv, sv_flag, n_bins, n_comps, &
                    n_aer1,d_aer1,sig_aer1,dmina,dmaxa,mass_frac_aer1,molw_core1, &
                    density_core1, nu_core1, kappa_core1, org_content1, molw_org1, &
                    kappa_org1, density_org1, delta_h_vap1,nu_org1, log_c_star1,sce_flag)
        
        ! This code writes the SCE variables to the BMM            
        if(sce_flag.gt.0) then
            ! send the SCE arrays, and use the local BMM arrays to map
            ! parcel1 here are the sce module vars. They are written to the bmm version of
            ! parcel1
            call write_sce_to_bmm(parcel1%n_bin_mode,parcel1%n_bin_modew,parcel1%n_binst,&
                    parcel1%n_modes, parcel1%n_comps, parcel1%n_comps+parcel1%imoms, &
                    parcel1%ice_flag, &
                    parcel1%npart, parcel1%moments, parcel1%mbin, parcel1%vel, &
                    parcel1%indexc, parcel1%ecoll, &
                    parcel1%mbinedges,adiabatic_prof)
        endif        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! run the model                                                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        io1%new_file=.true.
        call bmm_driver(sce_flag,hm_flag,break_flag,mode1_flag,mode2_flag) 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    end program main



