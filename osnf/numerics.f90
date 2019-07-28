module numerics
    use numerics_type
    implicit none
    real(wp), parameter :: LN2 = 6.931471805599453094172321214581e-1_wp
    real(wp), parameter :: A0 = 1.1975323115670912564578e0_wp
    real(wp), parameter :: A1 = 4.7072688112383978012285e1_wp
    real(wp), parameter :: A2 = 6.9706266534389598238465e2_wp
    real(wp), parameter :: A3 = 4.8548868893843886794648e3_wp
    real(wp), parameter :: A4 = 1.6235862515167575384252e4_wp
    real(wp), parameter :: A5 = 2.3782041382114385731252e4_wp
    real(wp), parameter :: A6 = 1.1819493347062294404278e4_wp
    real(wp), parameter :: A7 = 8.8709406962545514830200e2_wp
    real(wp), parameter :: B0 = 1.0000000000000000000e0_wp
    real(wp), parameter :: B1 = 4.2313330701600911252e1_wp
    real(wp), parameter :: B2 = 6.8718700749205790830e2_wp
    real(wp), parameter :: B3 = 5.3941960214247511077e3_wp
    real(wp), parameter :: B4 = 2.1213794301586595867e4_wp
    real(wp), parameter :: B5 = 3.9307895800092710610e4_wp
    real(wp), parameter :: B6 = 2.8729085735721942674e4_wp
    real(wp), parameter :: B7 = 5.2264952788528545610e3_wp
    real(wp), parameter :: C0 = 1.42343711074968357734e0_wp
    real(wp), parameter :: C1 = 4.63033784615654529590e0_wp
    real(wp), parameter :: C2 = 5.76949722146069140550e0_wp
    real(wp), parameter :: C3 = 3.64784832476320460504e0_wp
    real(wp), parameter :: C4 = 1.27045825245236838258e0_wp
    real(wp), parameter :: C5 = 2.41780725177450611770e-1_wp
    real(wp), parameter :: C6 = 2.27238449892691845833e-2_wp
    real(wp), parameter :: C7 = 7.74545014278341407640e-4_wp
    real(wp), parameter :: D0 = 1.4142135623730950488016887e0_wp
    real(wp), parameter :: D1 = 2.9036514445419946173133295e0_wp
    real(wp), parameter :: D2 = 2.3707661626024532365971225e0_wp
    real(wp), parameter :: D3 = 9.7547832001787427186894837e-1_wp
    real(wp), parameter :: D4 = 2.0945065210512749128288442e-1_wp
    real(wp), parameter :: D5 = 2.1494160384252876777097297e-2_wp
    real(wp), parameter :: D6 = 7.7441459065157709165577218e-4_wp
    real(wp), parameter :: D7 = 1.4859850019840355905497876e-9_wp
    real(wp), parameter :: E0 = 6.65790464350110377720e0_wp
    real(wp), parameter :: E1 = 5.46378491116411436990e0_wp
    real(wp), parameter :: E2 = 1.78482653991729133580e0_wp
    real(wp), parameter :: E3 = 2.96560571828504891230e-1_wp
    real(wp), parameter :: E4 = 2.65321895265761230930e-2_wp
    real(wp), parameter :: E5 = 1.24266094738807843860e-3_wp
    real(wp), parameter :: E6 = 2.71155556874348757815e-5_wp
    real(wp), parameter :: E7 = 2.01033439929228813265e-7_wp
    real(wp), parameter :: F0 = 1.414213562373095048801689e0_wp
    real(wp), parameter :: F1 = 8.482908416595164588112026e-1_wp
    real(wp), parameter :: F2 = 1.936480946950659106176712e-1_wp
    real(wp), parameter :: F3 = 2.103693768272068968719679e-2_wp
    real(wp), parameter :: F4 = 1.112800997078859844711555e-3_wp
    real(wp), parameter :: F5 = 2.611088405080593625138020e-5_wp
    real(wp), parameter :: F6 = 2.010321207683943062279931e-7_wp
    real(wp), parameter :: F7 = 2.891024605872965461538222e-15_wp    
	interface
		function find_pos(xarr,x)
		    use numerics_type
            real(wp), dimension(:), intent(in) :: xarr
            real(wp), intent(in) :: x
            integer(i4b) :: find_pos
		end function find_pos
	end interface
!
	interface
		subroutine poly_int(xarr,yarr,x,y,dy)
		use numerics_type
		real(wp), dimension(:), intent(in) :: xarr,yarr
		real(wp), intent(in) :: x
		real(wp), intent(out) :: y,dy
		end subroutine poly_int
	end interface
!
    interface
        subroutine tridiagonal(a,b,c,r,x)
            use numerics_type
            implicit none
            real(wp), dimension(:), intent(in) :: a,b,c,r
            real(wp), dimension(:), intent(inout) :: x
        end subroutine tridiagonal
    end interface
!
	interface
		function zeroin(ax,bx,f,tol)
		use numerics_type
		real(wp), intent(in) :: ax,bx,tol
		real(wp) :: zeroin
		interface 
		    function f(x)
		    use numerics_type
		    real(wp), intent(in) :: x
		    real(wp) :: f
		    end function f
		end interface
		end function zeroin
	end interface
!
	interface assert_eq
		module procedure assert_eq2,assert_eq3,assert_eq4,assert_eqn
	end interface
	interface imaxloc
		module procedure imaxloc_r,imaxloc_i
	end interface
!
    interface
        function r1mach(I)
            use numerics_type
            real(sp) :: r1mach
            integer :: I
        end function r1mach
    end interface
!
    interface
        function d1mach(I)
            use numerics_type
            real(dp) :: d1mach
            integer :: I
        end function d1mach
    end interface
!
    interface
      function dfsid1(func,x,h0,delta,err)
	    use numerics_type
        implicit none
        real(wp), intent(in) :: x,h0,delta
        real(wp), intent(inout) :: err
        real(wp) :: dfsid1
        interface 
            function func(x)
            use numerics_type
            real(wp), intent(in) :: x
            real(wp) :: func
            end function func       
        end interface        
      end function dfsid1 
    end interface
!
    interface
        subroutine erfinv(x,ans) 
            use, intrinsic :: IEEE_ARITHMETIC
            use numerics_type
            implicit none
            real(wp), intent(in) :: x
            real(wp), intent(inout) :: ans
        end subroutine erfinv
    end interface 
!
    interface
        subroutine vode_integrate(ystart,x1,x2,eps,h1,hmin,fdash)
            use numerics_type
            real(wp), intent(in) :: x1,x2,eps,h1,hmin
            real(wp), intent(inout), dimension(:) :: ystart
            interface 
                subroutine fdash(p,z,dzdp)
                    use numerics_type
                    implicit none
                    real(wp), intent(in) :: p
                    real(wp), dimension(:), intent(in) :: z
                    real(wp), dimension(:), intent(out) :: dzdp
                end subroutine fdash        
            end interface
    	end subroutine vode_integrate
    end interface
!
    interface
          subroutine xsetf (mflag)
            integer mflag
          end subroutine xsetf
    end interface
!	
	interface dvode
    	! single precision version
        subroutine svode (f, neq, y, t, tout, itol, rtol, atol, itask, &
                 istate, iopt, rwork, lrw, iwork, liw, jac, mf,&
                 rpar, ipar)
        external f, jac
        real y, t, tout, rtol, atol, rwork, rpar
        integer neq, itol, itask, istate, iopt, lrw, iwork, liw,&
             mf, ipar
        dimension y(*), rtol(*), atol(*), rwork(lrw), iwork(liw),&
               rpar(*), ipar(*)
        end subroutine svode
    
        ! double precision version
        subroutine dvode (f, neq, y, t, tout, itol, rtol, atol, itask, &
                 istate, iopt, rwork, lrw, iwork, liw, jac, mf,&
                 rpar, ipar)
        external f, jac
        double precision y, t, tout, rtol, atol, rwork, rpar
        integer neq, itol, itask, istate, iopt, lrw, iwork, liw,&
             mf, ipar
        dimension y(*), rtol(*), atol(*), rwork(lrw), iwork(liw),&
               rpar(*), ipar(*)
        end subroutine dvode
	end interface
!	
	interface fmin
    	! single precision version
        real function sfmin(ax,bx,f,tol)
        real ax,bx,f,tol
        external f
        end function sfmin
    
        ! double precision version
        double precision function fmin(ax,bx,f,tol)
        double precision ax,bx,f,tol
        external f
        end function fmin
	end interface
!
    contains
!
	function assert_eq2(n1,n2,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1,n2
        integer :: assert_eq2
        if (n1 == n2) then
            assert_eq2=n1
        else
            write (*,*) 'numerics: an assert_eq failed:', &
                string
            stop 'program terminated by assert_eq2'
        end if
	end function assert_eq2
!
	function assert_eq3(n1,n2,n3,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1,n2,n3
        integer :: assert_eq3
        if (n1 == n2 .and. n2 == n3) then
            assert_eq3=n1
        else
            write (*,*) 'numerics: an assert_eq failed:', &
                string
            stop 'program terminated by assert_eq3'
        end if
	end function assert_eq3
!
	function assert_eq4(n1,n2,n3,n4,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1,n2,n3,n4
        integer :: assert_eq4
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
            assert_eq4=n1
        else
            write (*,*) 'numerics: an assert_eq failed:', &
                string
            stop 'program terminated by assert_eq4'
        end if
	end function assert_eq4
!
	function assert_eqn(nn,string)
        character(len=*), intent(in) :: string
        integer, dimension(:), intent(in) :: nn
        integer :: assert_eqn
        if (all(nn(2:) == nn(1))) then
            assert_eqn=nn(1)
        else
            write (*,*) 'numerics: an assert_eq failed:', &
                string
            stop 'program terminated by assert_eqn'
        end if
	end function assert_eqn
!
	function imaxloc_r(arr)
	    use numerics_type
        real(wp), dimension(:), intent(in) :: arr
        integer(i4b) :: imaxloc_r
        integer(i4b), dimension(1) :: imax
        imax=maxloc(arr(:))
        imaxloc_r=imax(1)
	end function imaxloc_r
!
	function imaxloc_i(iarr)
	    use numerics_type
        integer(i4b), dimension(:), intent(in) :: iarr
        integer(i4b), dimension(1) :: imax
        integer(i4b) :: imaxloc_i
        imax=maxloc(iarr(:))
        imaxloc_i=imax(1)
	end function imaxloc_i
!
	function iminloc(arr)
	    use numerics_type
        real(wp), dimension(:), intent(in) :: arr
        integer(i4b), dimension(1) :: imin
        integer(i4b) :: iminloc
        imin=minloc(arr(:))
        iminloc=imin(1)
	end function iminloc 
!
	subroutine numerics_error(string)
        character(len=*), intent(in) :: string
        write (*,*) 'numerics_error: ',string
        stop 'program terminated by numerics_error'
	end subroutine numerics_error
	


	
end module numerics
