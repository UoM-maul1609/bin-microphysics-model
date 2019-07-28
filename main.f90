	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2019
	!>@brief
	!>A collection of numerical routines for use in models
	!>that can be distributed and modified


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main programme for testing 

    program main
        use numerics_type
        use numerics, only : dvode, dfsid1, zeroin,fmin, vode_integrate, &
                tridiagonal, erfinv
        use random, only : random_normal
        use hypergeo, only : hygfx
        implicit none
        real(wp) :: machep
        real(wp) :: ax,bx,tol,f1,f2,func,h0,err
        
        real(wp), dimension(3) :: atol, y,rtol
        real(wp) :: t,tout
        real(wp), dimension(67) :: rwork
        integer(i4b), dimension(33) :: iwork
        integer(i4b) :: neq, itol, itask,istate,iopt,lrw,liw,mf, &
                                    iout,i3,i4
        real(wp), allocatable,dimension(:) :: rpar
        integer(i4b), allocatable,dimension(:) :: ipar
        
        real(wp), dimension(1) :: ystart
        real(wp) :: x1,x2,eps,h1,hmin,f,ff1,aa,bb
        real(wp), dimension(10000001) :: b,x,r
        real(wp), dimension(10000000) :: a,c
        integer(i4b), allocatable, dimension(:) :: seed
        integer(i4b) :: l
        real(wp) :: rr

        external func 
        external func2
        external fex
        external jex

        machep=epsilon(tol)
        
        tol=1.e-8_wp
        ax=2._wp
        bx=6._wp
        f1=zeroin(ax,bx,func,tol) ! see zbrent
        print *,f1
        
        h0=1.e-3_wp
        f2=f1
        f1=dfsid1(func,f2,h0,tol,err) ! see dfridr
        print *,f1 

        f1=fmin(ax,bx,func,tol) ! see brent
        print *,f1,ax,bx
        
        call xsetf(0)

        neq = 3
        y(1) = 1.0e0_wp
        y(2) = 0.0e0_wp
        y(3) = 0.0e0_wp
        t = 0.0e0_wp
        tout = 0.4e0_wp
        itol = 2
        rtol = 1.e-4_wp
        atol(1) = 1.e-8_wp
        atol(2) = 1.e-14_wp
        atol(3) = 1.e-6_wp
        itask = 1
        istate = 1
        iopt = 0
        lrw = 67
        liw = 33
        mf = 21
        
        do iout = 1,12
            call dvode(fex,neq,y,t,tout,itol,rtol,atol,itask,istate, &
                    iopt,rwork,lrw,iwork,liw,jex,mf,rpar,ipar)
            write(6,20)t,y(1),y(2),y(3)
        20  format(' at t =',e12.4,'   y =',3e14.6)
            
            if (istate .lt. 0) then
                write(6,90) istate
            90  format(///' error halt: istate =',i3)
                stop            
            endif
            
            tout = tout*10._wp
        enddo
        write(6,60) iwork(11),iwork(12),iwork(13),iwork(19), &
                iwork(20),iwork(21),iwork(22)           
    60  format(/' no. steps =',i4,'   no. f-s =',i4, &
           '   no. j-s =',i4,'   no. lu-s =',i4/ &
           '  no. nonlinear iterations =',i4/ &
           '  no. nonlinear convergence failures =',i4/ &
           '  no. error test failures =',i4/)



        ! now test integration by vode_integrate
        ystart(1)=0._wp
        x1=0._sp
        x2=pi/2._wp
        eps=1.e-4_wp
        h1=0.01_wp
        hmin=0.0_wp
        call vode_integrate(ystart,x1,x2,eps,h1,hmin,func2)
        print *,ystart
        
        ! tridiagonal solver
        a=1._wp
        b=2._wp
        c=1._wp
        r=1._wp
        call tridiagonal(a,b,c,r,x)
        
        ! inverse erf
        call erfinv(0.5_wp,f)
        print *,f,erf(f)
        
        ! hypergeometric function
        aa=1._wp
        bb=4._wp+2._wp*2.5_wp+0.25_wp
        call hygfx(aa, bb, real(1,sp)+2.5_wp+1.0_wp,0.5_wp,ff1)
        print *,ff1
        
        ! random number
		call random_seed(size=l)
		allocate(seed(1:l))
		seed(:)=2
		call random_seed(put=seed)
        
        do l=1,10
            rr=random_normal()
            print *,rr
        enddo
        
        
    end program main



    function func(x)
        use numerics_type
        use numerics
        implicit none
        real(wp) :: func
        real(wp), intent(in) :: x
        
        func=cos(x)
    end function func


    subroutine fex (neq, t, y, ydot, rpar, ipar)
        use numerics_type
        real(wp) :: rpar, t
        real(wp), dimension(neq) :: y,ydot
        ydot(1) = -.04e0_wp*y(1) + 1.e4_wp*y(2)*y(3)
        ydot(3) = 3.e7_wp*y(2)*y(2)
        ydot(2) = -ydot(1) - ydot(3)
        return
    end subroutine fex

    subroutine func2 (t, y, ydot)
        use numerics_type
        implicit none
        real(wp), intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: ydot
        ydot(1) = cos(t)
    end subroutine func2

    subroutine jex (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
        use numerics_type
        real(wp) :: rpar, t
        real(wp), dimension(neq) :: y
        real(wp), dimension(nrpd,neq) :: pd
        pd(1,1) = -.04e0_wp
        pd(1,2) = 1.e4_wp*y(3)
        pd(1,3) = 1.e4_wp*y(2)
        pd(2,1) = .04e0_wp
        pd(2,3) = -pd(1,3)
        pd(3,2) = 6.e7_wp*y(2)
        pd(2,2) = -pd(1,2) - pd(3,2)
        return
    end subroutine jex    

