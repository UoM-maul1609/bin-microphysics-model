    !>@author
    !>Paul Connolly, The University of Manchester
    !>@brief
    !>Modified anomalous diffraction approximation (ADA / "ADT") of
    !>Mitchell (2000, J. Atmos. Sci., 57, 1311-1326), "Parameterization of the
    !>Mie extinction and absorption coefficients for water clouds".
    !>
    !>Unlike the paper's Eqs (19)-(23), which integrate analytically over a
    !>gamma size distribution, this module evaluates the single-particle
    !>efficiencies Qext(D) and Qabs(D) from Eqs (1)-(18) directly for one
    !>diameter at a time. That's the natural fit for a bin microphysics
    !>scheme: sum Qext(D_i)*(pi/4)*D_i^2*N_i over the bins instead of doing
    !>the closed-form gamma-distribution integral.
    module adt_scattering_mod
    use numerics_type
    use numerics, only : find_pos, poly_int
    use refractive_indices_mod
    implicit none

    private
    public :: qext_qabs_adt, refractive_index_lookup, beta_ext_adt_bin, &
              cloud_albedo_conservative

    ! default asymmetry parameter for water cloud droplets (Mitchell 2000
    ! uses g=0.90 in his own Table 1/2 comparisons; Slingo (1989) and others
    ! commonly use ~0.85-0.86). Treated as a constant since ADT does not
    ! provide g - see discussion after Eq (1) in the paper.
    real(wp), parameter :: g_asym_default=0.90_wp

    ! dispersion parameter for the tunnelling gamma-function form, Eq (7)-(9).
    ! m=0.5 was found by Mitchell (2000) to best reproduce the Mie dispersion.
    real(wp), parameter :: m_disp=0.5_wp
    ! recommended edge-effect coefficient, Eq (17) (a6=1 gives the better fit
    ! - see discussion after Eq (17)/Fig 7 in the paper).
    real(wp), parameter :: a6=1._wp
    ! local copies so this module doesn't depend on bmm's own parameter list
    real(wp), parameter :: onequarter=0.25_wp, twothirds=2._wp/3._wp

    contains

    !>@author
    !>Paul Connolly, The University of Manchester
    !>@brief
    !>look up the real and imaginary refractive index of water or ice at a
    !>given wavelength (m), by linear interpolation of the tables in
    !>refractive_indices_mod (values are clipped at the table edges, matching
    !>the behaviour of adt_scattering.py / refractive_indices_s.py).
    !>@param[in] wvl: wavelength (m)
    !>@param[in] is_ice: .true. for the ice table, .false. for water
    !>@param[out] nr, ni: real and imaginary parts of the refractive index
    !> (ni is returned positive, i.e. n = nr - i*ni)
    subroutine refractive_index_lookup(wvl,is_ice,nr,ni)
    implicit none
    real(wp), intent(in) :: wvl
    logical, intent(in) :: is_ice
    real(wp), intent(out) :: nr, ni

    integer(i4b) :: iloc, n
    real(wp) :: dy

        if(is_ice) then
            n=size(lam_ice_tab)
            if(wvl <= lam_ice_tab(1)) then
                nr=nr_ice_tab(1); ni=ni_ice_tab(1)
            else if(wvl >= lam_ice_tab(n)) then
                nr=nr_ice_tab(n); ni=ni_ice_tab(n)
            else
                iloc=find_pos(lam_ice_tab,wvl)
                iloc=max(1,min(iloc,n-1))
                call poly_int(lam_ice_tab(iloc:iloc+1),nr_ice_tab(iloc:iloc+1),wvl,nr,dy)
                call poly_int(lam_ice_tab(iloc:iloc+1),ni_ice_tab(iloc:iloc+1),wvl,ni,dy)
            endif
        else
            n=size(lam_h2o_tab)
            if(wvl <= lam_h2o_tab(1)) then
                nr=nr_h2o_tab(1); ni=ni_h2o_tab(1)
            else if(wvl >= lam_h2o_tab(n)) then
                nr=nr_h2o_tab(n); ni=ni_h2o_tab(n)
            else
                iloc=find_pos(lam_h2o_tab,wvl)
                iloc=max(1,min(iloc,n-1))
                call poly_int(lam_h2o_tab(iloc:iloc+1),nr_h2o_tab(iloc:iloc+1),wvl,nr,dy)
                call poly_int(lam_h2o_tab(iloc:iloc+1),ni_h2o_tab(iloc:iloc+1),wvl,ni,dy)
            endif
        endif
        ni=abs(ni)
    end subroutine refractive_index_lookup



    !>@author
    !>Paul Connolly, The University of Manchester
    !>@brief
    !>single-particle extinction and absorption efficiencies from the
    !>modified ADA of Mitchell (2000), Eqs (1)-(18).
    !>@param[in] d: particle diameter (m)
    !>@param[in] wvl: wavelength (m)
    !>@param[in] nr, ni: real and imaginary refractive index (n=nr-i*ni)
    !>@param[out] qext, qabs: extinction and absorption efficiencies
    subroutine qext_qabs_adt(d,wvl,nr,ni,qext,qabs)
    implicit none
    real(wp), intent(in) :: d, wvl, nr, ni
    real(wp), intent(out) :: qext, qabs

    real(wp) :: x, k, a1, qabs_ada, c1, ra_max, eps0, kmax, c2, &
                r_ext, c3, qext_ada, qedge, tun
    complex(wp) :: t, kt, ncom

        ! Mie size parameter and the "k=D/lambda" variable used for the
        ! resonance/tunnelling terms (distinct quantities - see Eqs (7),(14))
        x=pi*d/wvl
        k=d/wvl

        ! --- Qabs,ADA (Eq 2), using V/P=(2/3)*D for a sphere ---
        qabs_ada=1._wp-exp(-4._wp*pi*ni*(twothirds*d)/wvl)

        ! --- internal reflection/refraction correction, Eqs (3),(5) ---
        a1=0.25_wp+0.25_wp*exp(-1167._wp*ni)
        c1=a1*(1._wp-qabs_ada)

        ! --- tunnelling contribution, Eqs (6)-(9) ---
        ra_max=0.7393_wp*nr-0.6069_wp
        eps0=0.25_wp+0.6_wp*(1._wp-exp(-8._wp*pi*ni/3._wp))**2
        kmax=m_disp/eps0
        if(kmax>0._wp .and. k>0._wp) then
            tun=(k**m_disp*exp(-eps0*k))/(kmax**m_disp*exp(-m_disp))
        else
            tun=0._wp
        endif
        c2=ra_max*tun

        ! --- absorption efficiency, Eq (10) ---
        qabs=(1._wp+c1+c2)*qabs_ada
        qabs=max(qabs,0._wp)

        ! --- extinction: tunnelling term for Qext, Eqs (11)-(12) ---
        r_ext=ra_max/2._wp
        c3=r_ext*tun

        ! --- Qext,ADA via the van de Hulst K(t) function, Eqs (14)-(15) ---
        ncom=cmplx(nr,-ni,kind=wp)
        t=cmplx(0._wp,1._wp,kind=wp)*2._wp*pi*d*(ncom-cmplx(1._wp,0._wp,kind=wp))/wvl
        kt=0.5_wp+exp(-t)/t+(exp(-t)-1._wp)/t**2
        qext_ada=4._wp*real(kt,kind=wp)

        ! --- edge-effect term, Eq (17), with a6=1 ---
        qedge=a6*(1._wp-exp(-0.06_wp*x))*x**(-twothirds)

        ! --- total extinction efficiency, Eq (18) ---
        qext=(1._wp+c3)*qext_ada+qedge
        qext=max(qext,0._wp)

    end subroutine qext_qabs_adt



    !>@author
    !>Paul Connolly, The University of Manchester
    !>@brief
    !>bulk volume-extinction coefficient (m-1) for a bin/sectional size
    !>distribution, using the single-particle ADT Qext(D_i) in place of the
    !>fixed geometric-optics value Qext=2 used elsewhere in the model:
    !>    beta_ext = rhod * sum_i Qext(D_i) * (pi/4) * D_i^2 * N_i
    !>where N_i has units of #/kg (as in parcel1%npart) and rhod is the
    !>dry-air density (kg/m3), giving beta_ext in m-1, directly comparable
    !>with the existing beta_ext = 2*pi*rhod*sum(0.25*dw**2*npart) diagnostic.
    !>@param[in] d: array of bin diameters (m)
    !>@param[in] npart: array of bin number concentrations (#/kg)
    !>@param[in] rhod: dry air density (kg/m3)
    !>@param[in] wvl: wavelength (m)
    !>@param[in] is_ice: .true. to use the ice refractive index table
    !>@param[out] beta_ext: bulk extinction coefficient (m-1)
    !>@param[out] beta_abs: bulk absorption coefficient (m-1), optional-use
	subroutine beta_ext_adt_bin(d,npart,rhod,wvl,is_ice, &
								beta_ext,beta_abs,area)
	
		implicit none
		real(wp), dimension(:), intent(in) :: d,npart
		real(wp), intent(in) :: rhod,wvl
		logical, intent(in) :: is_ice
		real(wp), intent(out) :: beta_ext,beta_abs
		real(wp), dimension(:), intent(in), optional :: area
		integer(i4b) :: i
		real(wp) :: nr,ni,qext,qabs,aproj
	
		call refractive_index_lookup(wvl,is_ice,nr,ni)
	
		beta_ext=0._wp
		beta_abs=0._wp
		do i=1,size(d)
			if(d(i) <= 0._wp .or. npart(i) <= 0._wp) cycle
			call qext_qabs_adt(d(i),wvl,nr,ni,qext,qabs)
			if(present(area)) then
				aproj=area(i)
			else
				aproj=onequarter*pi*d(i)**2
			endif
			beta_ext=beta_ext + qext*aproj*npart(i)
			beta_abs=beta_abs + qabs*aproj*npart(i)
		enddo
		beta_ext=rhod*beta_ext
		beta_abs=rhod*beta_abs
	end subroutine beta_ext_adt_bin


    !>@author
    !>Paul Connolly, The University of Manchester
    !>@brief
    !>plane-parallel cloud albedo for conservative (non-absorbing) scattering
    !>with overhead sun, from Seinfeld & Pandis ("Atmospheric Chemistry and
    !>Physics"), attributed there to Lacis and Hansen (1974):
    !>    R = tau*(1-g) / (2 + tau*(1-g))
    !>where tau*(1-g) is the similarity-scaled optical depth. Since this
    !>ignores absorption, tau should come from the ADT beta_ext (not
    !>beta_abs) integrated over the cloud depth: tau = beta_ext * depth.
    !>@param[in] tau: visible-wavelength cloud optical depth (-)
    !>@param[in] g: asymmetry parameter (-); omit to use g_asym_default=0.90
    !>@param[out] albedo: plane-parallel cloud albedo (-)
    subroutine cloud_albedo_conservative(tau,albedo,g)
    implicit none
    real(wp), intent(in) :: tau
    real(wp), intent(out) :: albedo
    real(wp), intent(in), optional :: g

    real(wp) :: gg, tau_scaled

        gg=g_asym_default
        if(present(g)) gg=g

        tau_scaled=tau*(1._wp-gg)
        albedo=tau_scaled/(2._wp+tau_scaled)

    end subroutine cloud_albedo_conservative

    end module adt_scattering_mod
