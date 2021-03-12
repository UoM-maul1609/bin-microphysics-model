module numerics
    use numerics_type
    implicit none
    private 
    public find_pos, poly_int, tridiagonal, zeroin,assert_eq, r1mach, d1mach, &
        dfsid1, erfinv, vode_integrate, xsetf, dvode, fmin, assert_eq2, assert_eq3, &
        assert_eq4, assert_eqn, imaxloc, iminloc, numerics_error, &
        quad2d_romb, gammainc, gammainc_scal, invgammainc, &
        LN2, &
        A0, A1, A2, A3, A4, A5, A6, A7, &
        B0, B1, B2, B3, B4, B5, B6, B7, &
        C0, C1, C2, C3, C4, C5, C6, C7, &
        D0, D1, D2, D3, D4, D5, D6, D7, &
        E0, E1, E2, E3, E4, E5, E6, E7, &
        F0, F1, F2, F3, F4, F5, F6, F7
    real(wp) :: xsav

    integer(i4b), parameter :: npar_arith=16,npar_arith2=8
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

    interface gammainc
        module procedure gammainc_scal, gammainc_vec
    end interface
    
    interface gcf
        module procedure gcf_scal, gcf_vec
    end interface
    interface gser
        module procedure gser_scal, gser_vec
    end interface

    contains
    ! series solutions for the incomplete gamma function
	function gser_scal(a,x,gln)
	implicit none
	real(wp), intent(in) :: a,x
	real(wp), optional, intent(out) :: gln
	real(wp) :: gser_scal
	integer(i4b), parameter :: itmax=100
	real(wp), parameter :: eps=epsilon(x)
	integer(i4b) :: n
	real(wp) :: ap,del,summ
	if (x == 0.0_wp) then
		gser_scal=0.0_wp
		return
	end if
	ap=a
	summ=1.0_wp/a
	del=summ
	do n=1,itmax
		ap=ap+1.0_wp
		del=del*x/ap
		summ=summ+del
		if (abs(del) < abs(summ)*eps) exit
	end do
	if (n > itmax) call numerics_error('a too large, itmax too small in gser_scal')
	if (present(gln)) then
		gln=log_gamma(a)
		gser_scal=summ*exp(-x+a*log(x)-gln)
	else
		gser_scal=summ*exp(-x+a*log(x)-log_gamma(a))
	end if
	end function gser_scal


    ! series solutions for the incomplete gamma function
	function gser_vec(a,x,gln)
	implicit none
	real(wp), dimension(:), intent(in) :: a,x
	real(wp), dimension(:), optional, intent(out) :: gln
	real(wp), dimension(size(a)) :: gser_vec
	integer(i4b), parameter :: itmax=100
	real(wp), parameter :: eps=epsilon(x)
	integer(i4b) :: n
	real(wp), dimension(size(a)) :: ap,del,summ
	logical(lgt), dimension(size(a)) :: converged,zero
	n=assert_eq2(size(a),size(x),'gser_vec')
	zero=(x == 0.0_wp)
	where (zero) gser_vec=0.0_wp
	ap=a
	summ=1.0_wp/a
	del=summ
	converged=zero
	do n=1,itmax
		where (.not. converged)
			ap=ap+1.0_wp
			del=del*x/ap
			summ=summ+del
			converged = (abs(del) < abs(summ)*eps)
		end where
		if (all(converged)) exit
	end do
	if (n > itmax) call numerics_error('a too large, itmax too small in gser_vec')
	if (present(gln)) then
		if (size(gln) < size(a)) call &
			numerics_error('gser: not enough space for gln')
		gln=log_gamma(a)
		where (.not. zero) gser_vec=summ*exp(-x+a*log(x)-gln)
	else
		where (.not. zero) gser_vec=summ*exp(-x+a*log(x)-log_gamma(a))
	end if
	end function gser_vec


    ! continued fraction solution of incomplete gamma
	function gcf_scal(a,x,gln)
	implicit none
	real(wp), intent(in) :: a,x
	real(wp), optional, intent(out) :: gln
	real(wp) :: gcf_scal
	integer(i4b), parameter :: itmax=100
	real(wp), parameter :: eps=epsilon(x),fpmin=tiny(x)/eps
	integer(i4b) :: i
	real(wp) :: an,b,c,d,del,h
	if (x == 0.0_wp) then
		gcf_scal=1.0_wp
		return
	end if
	b=x+1.0_wp-a
	c=1.0_wp/fpmin
	d=1.0_wp/b
	h=d
	do i=1,itmax
		an=-i*(i-a)
		b=b+2.0_wp
		d=an*d+b
		if (abs(d) < fpmin) d=fpmin
		c=b+an/c
		if (abs(c) < fpmin) c=fpmin
		d=1.0_wp/d
		del=d*c
		h=h*del
		if (abs(del-1.0_wp) <= eps) exit
	end do
	if (i > itmax) call numerics_error('a too large, itmax too small in gcf_scal')
	if (present(gln)) then
		gln=log_gamma(a)
		gcf_scal=exp(-x+a*log(x)-gln)*h
	else
		gcf_scal=exp(-x+a*log(x)-log_gamma(a))*h
	end if
	end function gcf_scal


    ! continued fraction solution of incomplete gamma
	function gcf_vec(a,x,gln)
	implicit none
	real(wp), dimension(:), intent(in) :: a,x
	real(wp), dimension(:), optional, intent(out) :: gln
	real(wp), dimension(size(a)) :: gcf_vec
	integer(i4b), parameter :: itmax=100
	real(wp), parameter :: eps=epsilon(x),fpmin=tiny(x)/eps
	integer(i4b) :: i
	real(wp), dimension(size(a)) :: an,b,c,d,del,h
	logical(lgt), dimension(size(a)) :: converged,zero
	i=assert_eq2(size(a),size(x),'gcf_vec')
	zero=(x == 0.0_wp)
	where (zero)
		gcf_vec=1.0_wp
	elsewhere
		b=x+1.0_wp-a
		c=1.0_wp/fpmin
		d=1.0_wp/b
		h=d
	end where
	converged=zero
	do i=1,itmax
		where (.not. converged)
			an=-i*(i-a)
			b=b+2.0_wp
			d=an*d+b
			d=merge(fpmin,d, abs(d)<fpmin )
			c=b+an/c
			c=merge(fpmin,c, abs(c)<fpmin )
			d=1.0_wp/d
			del=d*c
			h=h*del
			converged = (abs(del-1.0_wp)<=eps)
		end where
		if (all(converged)) exit
	end do
	if (i > itmax) call numerics_error('a too large, itmax too small in gcf_vec')
	if (present(gln)) then
		if (size(gln) < size(a)) call &
			numerics_error('gser: not enough space for gln')
		gln=log_gamma(a)
		where (.not. zero) gcf_vec=exp(-x+a*log(x)-gln)*h
	else
		where (.not. zero) gcf_vec=exp(-x+a*log(x)-log_gamma(a))*h
	end if
	end function gcf_vec


    ! calculates this incomplete gamma function
	function gammainc_scal(a,x)
	implicit none
	real(wp), intent(in) :: a,x
	real(wp) :: gammainc_scal
	call assert1( x >= 0.0_wp,  a > 0.0_wp, 'gamminc_scal args')
	if (x<a+1.0_wp) then
		gammainc_scal=gser(a,x)
	else
		gammainc_scal=1.0_wp-gcf(a,x)
	end if
	end function gammainc_scal


	function gammainc_vec(a,x)
	implicit none
	real(wp), dimension(:), intent(in) :: a,x
	real(wp), dimension(size(x)) :: gammainc_vec
	logical(lgt), dimension(size(x)) :: mask
	integer(i4b) :: ndum
	ndum=assert_eq2(size(a),size(x),'gamminc_vec')
	call assert1( all(x >= 0.0_wp),  all(a > 0.0_wp), 'gamminc_vec args')
	mask = (x<a+1.0_wp)
	gammainc_vec=merge(gser(a,merge(x,0.0_wp,mask)), &
		1.0_wp-gcf(a,merge(x,0.0_wp,.not. mask)),mask)
	end function gammainc_vec



    function invgammainc(p,a) 
        implicit none
        real(wp), intent(in) :: a
        real(wp), intent(inout) :: p
        real(wp) :: invgammainc

        integer(i4b) :: j
        real(wp) :: x,err,t,u,pp,lna1,afac,gln,a1
        real(wp), parameter :: eps=1.e-8_wp
    
        a1=a-1.0_wp
        gln=log_gamma(a)
        if(a<=0.0_wp) print *, "a must be positive"
        if(p>=1.0_wp) then
            invgammainc=max(100._wp,a+100._wp*sqrt(a))
            return
        endif
        if(p<=0.0_wp) then
            invgammainc=0.0_wp
            return
        endif
        if(a>1.0_wp) then
            lna1=log(a1)
            afac=exp(a1*(lna1-1._wp)-gln)
            if(p<0.5_wp) then
                pp=p
            else
                pp=1._wp-p
            endif
            t=sqrt(-2._wp*log(pp))
            x=(2.30753_wp+t*0.27061_wp)/(1._wp+t*(0.99229_wp+t*0.04481_wp))-t
            if(p<0.5_wp) x=-x
            x=max(1.e-3_wp,a*(1._wp-1._wp/(9._wp*a)-x/(3._wp*sqrt(a)))**3)
        else
            t=1.0_wp-a*(0.253_wp+a*0.12_wp)
            if(p<t) then 
                x=(p/t)**(1._wp/a)
            else
                x=1._wp-log(1._wp-(p-t)/(1._wp-t))
            endif
        
        endif
    
        do j=0,11
            if(x<=0.0_wp) then
                invgammainc=0.0_wp
                return
            endif
            err=gammainc(a,x)-p
            if(a>1.) then
                t=afac*exp(-(x-a1)+a1*(log(x)-lna1))
            else
                t=exp(-x+a1*log(x)-gln)
            endif
            u=err/t
            t=u/(1._wp-0.5_wp*min(1._wp,u*((a-1._wp)/x-1._wp)))
            x=x-t
            if(x<=0._wp) x=0.5_wp*(x+t)
            if(abs(t) < eps*x) exit
        enddo
    
        invgammainc=x
        return 
    end function invgammainc

	recursive function romb(func1,a,b)
		use numerics_type

	implicit none
	real(wp), intent(in) :: a,b
	real(wp) :: romb
	interface
		function func1(x)
		use numerics_type
		real(wp), dimension(:), intent(in) :: x
		real(wp), dimension(size(x)) :: func1
		end function func1
	end interface
	integer(i4b), parameter :: jmax=20,jmaxp=jmax+1,k=5,km=k-1
	real(wp), parameter :: eps=3.0e-6_wp
	real(wp), dimension(jmaxp) :: h,s
	real(wp) :: dromb
	integer(i4b) :: j
	h(1)=1.0_wp
	do j=1,jmax
		call trapezoid(func1,a,b,s(j),j)
		if (j >= k) then
			call poly_int(h(j-km:j),s(j-km:j),0.0_wp,romb,dromb)
			if (abs(dromb) <= eps*abs(romb)) return
		end if
		s(j+1)=s(j)
		h(j+1)=0.25_wp*h(j)
	end do
	call numerics_error('romb: too many steps')
	end function romb

	recursive subroutine trapezoid(func1,a,b,s,n)
	use numerics_type
	implicit none
	real(wp), intent(in) :: a,b
	real(wp), intent(inout) :: s
	integer(i4b), intent(in) :: n
	interface
		function func1(x)
		use numerics_type
		real(wp), dimension(:), intent(in) :: x
		real(wp), dimension(size(x)) :: func1
		end function func1
	end interface
	real(wp) :: del,fsum
	integer(i4b) :: it
	if (n == 1) then
		s=0.5_wp*(b-a)*sum(func1( (/ a,b /) ))
	else
		it=2**(n-2)
		del=(b-a)/it
		fsum=sum(func1(arithmetic_prog(a+0.5_wp*del,del,it)))
		s=0.5_wp*(s+del*fsum)
	end if
	end subroutine trapezoid




	subroutine quad2d_romb(func1,y1,y2,x1,x2,ss)
		use numerics_type
	implicit none
	real(wp), intent(in) :: x1,x2
	real(wp), intent(out) :: ss

	interface
		function func1(x,y)
		use numerics_type, only : wp
		implicit none
		real(wp), intent(in) :: x
		real(wp), dimension(:), intent(in) :: y
		real(wp), dimension(size(y)) :: func1
		end function func1
!

		function y1(x)
        use numerics_type, only : wp
        implicit none
		real(wp), intent(in) :: x
		real(wp) :: y1
		end function y1
!
		function y2(x)
        use numerics_type, only : wp
        implicit none
		real(wp), intent(in) :: x
		real(wp) :: y2
		end function y2

    end interface

	ss=romb(h,x1,x2)

    ! internal routine
    contains
         function f(y)
        use numerics_type, only : wp
        implicit none
         real(wp), dimension(:), intent(in) :: y
         real(wp), dimension(size(y)) :: f
         f=func1(xsav,y)
         end function f
        !bl
        function h(x)
        use numerics_type, only : wp
        implicit none
        real(wp), dimension(:), intent(in) :: x
        real(wp), dimension(size(x)) :: h
        integer(i4b) :: i
        do i=1,size(x)
            xsav=x(i)
            h(i)=romb(f,y1(xsav),y2(xsav))
        end do
        end function h

	end subroutine quad2d_romb

    function arithmetic_prog(first,increment,n)
    real(wp), intent(in) :: first,increment
    integer(i4b), intent(in) :: n
    real(wp), dimension(n) :: arithmetic_prog
    integer(i4b) :: k,k2
    real(wp) :: temp
    if (n > 0) arithmetic_prog(1)=first
    if (n <= npar_arith) then
        do k=2,n
            arithmetic_prog(k)=arithmetic_prog(k-1)+increment
        end do
    else ! if it is greater than 16, we do the first 8 as usual
        do k=2,npar_arith2
            arithmetic_prog(k)=arithmetic_prog(k-1)+increment
        end do
        ! then replicate in sizes of 8, 16, 32, and so on
        temp=increment*npar_arith2
        k=npar_arith2
        do
            if (k >= n) exit
            k2=k+k
            arithmetic_prog(k+1:min(k2,n))=temp+arithmetic_prog(1:min(k,n-k))
            temp=temp+temp
            k=k2
        end do
    end if
    end function arithmetic_prog

	subroutine assert1(l1,l2,string)
        character(len=*), intent(in) :: string
        logical, intent(in) :: l1,l2
        if (.not.((l1).and.(l2)) ) then
            write (*,*) 'numerics: an assert1 failed:', &
                string
            stop 'program terminated by assert1'
        end if
	end subroutine assert1
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
