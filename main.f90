    
    ! module to test double integration
    module test
    
    use numerics_type, only : i4b, wp
    integer (i4b), parameter :: dim_num=2
    real(wp), dimension(dim_num) :: a,b
    integer (i4b), dimension(dim_num) :: sub_num
    integer (i4b) :: it_max, eval_num, ind, i
    real(wp) :: tol, result
    
    real(wp), parameter :: lf=333550._wp, ttr=273.15_wp, cw=4200._wp, &
                gamma_liq=0.072_wp, DEcrit=0.2_wp, alphar=2.5_wp, alphai=2.5_wp, &
                pi=4.0_wp*atan(1.0_wp), cr=pi/6.0_wp*1000.0_wp, dr=3.0_wp, ar=836.0_wp, &
                br=0.8_wp, ci=0.04_wp, di=2.0_wp, ai=8.97_wp, bi=0.42_wp, &
                oneoversix=1._wp/6._wp, dtt=10.e-6_wp, rhow=1000._wp, rhoice=920._wp
    real(wp) :: phi, T, f, nr1, qr, ni, qi, n0r, lambda0r, n0i, lambda0i, &
        mrthresh, mrupper, miupper, vx, vy, lam_freeze, n0_freeze
    
    contains            
        
        
    subroutine calculate_parameters(n,q,c,d,alpha, n0, lambda)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: n,q,c,d,alpha
        real(wp), intent(inout) :: n0, lambda
        
        lambda=((c*gamma(1._wp+d+alpha)/gamma(1._wp+alpha))*N/q)**(1.0_wp/d)
        n0=n*lambda**(1.0_wp+alpha) / gamma(1.0_wp+alpha)
        
    end subroutine calculate_parameters
    
    
    function collisions_precip_mass(dim_num, x)
        use numerics_type, only : wp, i4b
        implicit none
        integer(i4b), intent(in) :: dim_num
        real(wp), dimension(dim_num), intent(in) :: x
        real(wp) :: collisions_precip_mass
        real(wp) :: diamr, diami, mr, mi


        mr=x(1)
        mi=x(2)
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci)**(1.0_wp/di)
        
        collisions_precip_mass=pi/4.0_wp*(diamr+diami)**2* &
            abs(ar*diamr**br-ai*diami**bi)*n0r*diamr**alphar* &
            exp(-lambda0r*diamr)*n0i*diami**alphai*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci*di)
    end function collisions_precip_mass
        


    function dintegral(x,y)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral
        real(wp) :: diamr, mr, vr
        real(wp), dimension(size(y)) :: mi, diami, delv, vi, k0, de, nfrag, nfrag_freeze1, &
            nfrag_freeze2


        mr=x
        mi=y
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci)**(1.0_wp/di)
        ! fall-speeds
        vr=ar*diamr**br
        vi=ai*diami**bi
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral=pi/4.0_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alphar* &
            exp(-lambda0r*diamr)*n0i*diami**alphai*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci*di)
        
        ! cke from equation 6
        k0=0.5_wp*(mr*mi/(mr+mi))*(vr-vi)**2
        ! de parameter
        de=k0/(gamma_liq*pi*diamr**2)
        ! number of fragments in spalsh
        nfrag=3.0_wp*max(de-decrit, 0.0_wp)
        ! number of fragments in splash that freeze due to mode 1
        nfrag_freeze1=nfrag*f
        ! number of fragments in splash that freeze due to mode 2
        nfrag_freeze2=nfrag*(1.0_wp-f)*phi
        
        dintegral=dintegral*nfrag_freeze2
    
    end function dintegral
    
    function dintegral2(x,y)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(size(y)) :: dintegral2
        real(wp) :: diamr, mr, vr, n,nt,nb,mb,mt
        real(wp), dimension(size(y)) :: mi, diami, delv, vi
        integer(i4b) :: i

        mr=x
        mi=y
        diamr=(mr/cr)**(1.0_wp/dr)
        diami=(mi/ci)**(1.0_wp/di)
        ! fall-speeds
        vr=ar*diamr**br
        vi=ai*diami**bi
        delv=abs(vr-vi)
        !delv=max((vx+vy)/8.0,abs(vx-vy))
        ! last bit is to convert to integral over m
        dintegral2=pi/4.0_wp*(diamr+diami)**2* &
            delv*n0r*diamr**alphar* &
            exp(-lambda0r*diamr)*n0i*diami**alphai*exp(-lambda0i*diami)* &
            (diamr**(1.0_wp-dr)) / (cr*dr)*(diami**(1.0_wp-di)) / (ci*di)
        
        do i=1,size(y)
            call calculate_mode1(mr,mi(i),260._wp,n,nt,nb,mb,mt)
            dintegral2(i)=dintegral2(i)*n
        enddo
        
    
    end function dintegral2
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mode 1 fragmentation integral over size distribution                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the number of fragments and their mass
    function integral_m1(x)
        use numerics_type, only : wp, i4b
        implicit none
		real(wp), dimension(:), intent(in) :: x
		real(wp), dimension(size(x)) :: integral_m1
		
		real(wp), dimension(size(x)) :: nfrag
		real(wp) :: n,nt,nb,mb,mt
		integer(i4b) :: i
		
		do i=1,size(x)
		    call calculate_mode1(x(i),0._wp,260._wp,n,nt,nb,mb,mt)
		    nfrag(i)=n
		enddo
		
    
        integral_m1=n0_freeze*exp(-lam_freeze*x)*x**alphar !*nfrag
    end function integral_m1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mode 1 fragmentation                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>calculates the number of fragments and their mass
    subroutine calculate_mode1(min1,min2,t,n,nt,nb,mb,mt)
        use numerics_type
        implicit none
        real(wp), intent(in) :: min1,min2,t
        real(wp), intent(inout) :: n, nt, nb, mb, mt
        real(wp) :: tc, dthresh, x, beta1,log10zeta, log10nabla, t0, zetab, nablab, tb0, &
            sigma, omega, m,d, fac1
        
        if((min2>min1).or.(min1<=6.55e-11_wp)) then
            ! the ice is more massive than the drop or drop small, don't do it
            n=0._wp
            nt=0._wp
            nb=0._wp
            return
        endif
    
        d = (6._wp*min1/rhow)**(1._wp/3._wp)
        tc=t-ttr
        dthresh = min(d,1.6e-3)
        x = log10(dthresh*1000._wp)
        
        ! table 3, phillips et al.
        beta1 = 0.
        log10zeta = 2.4268_wp*x*x*x + 3.3274_wp*x*x + 2.0783_wp*x + 1.2927_wp
        log10nabla = 0.1242_wp*x*x*x - 0.2316_wp*x*x - 0.9874_wp*x - 0.0827_wp
        t0 = -1.3999_wp*x*x*x - 5.3285_wp*x*x - 3.9847_wp*x - 15.0332_wp
        
        ! table 4, phillips et al. 
        zetab = -0.4651_wp*x*x*x - 1.1072_wp*x*x - 0.4539_wp*x+0.5137_wp
        nablab = 28.5888*x*x*x + 49.8504_wp*x*x + 22.4873_wp*x + 8.0481_wp
        tb0 = 13.3588_wp*x*x*x + 15.7432_wp*x*x - 2.6545_wp*x - 18.4875_wp
        
        sigma = min(max((d-50.e-6_wp)/10.e-6_wp,0._wp), 1._wp)
        omega = min(max((-3._wp-tc)/3._wp,0._wp),1._wp)
        
        n = sigma*omega*(10._wp**log10zeta *(10**log10nabla)**2) / &
            ((tc-t0)**2+(10._wp*log10nabla)**2+beta1*tc)
        
        ! total number of fragments
        n=n*d/dthresh
        ! number of large fragments
        nb = min(sigma*omega*(zetab*nablab**2/((tc-tb0)**2+nablab**2)),n)
        ! number of small fragments
        nt = n-nb
        
        m=oneoversix*rhow*pi*d**3
        
        ! mass of large fragments
        mb=0.4_wp*m
        
        ! mass of small fragments
        mt=oneoversix*rhoice*pi*dtt**3
        
        fac1=min((mt*nt+mb*nb)/min1,1._wp)
        nt = nt *fac1
        nb = nb *fac1
        n=nt+nb
        
    end subroutine calculate_mode1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    

    function limit1(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit1
        limit1=x
    end function limit1
!
    function limit2(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit2
        limit2=miupper
    end function limit2
    
    function limit3(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit3
        limit3=mrthresh
    end function limit3
!
    function limit4(x)
        use numerics_type, only : wp
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: limit4
        limit4=x
    end function limit4
    
    
    end module test
    
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
                tridiagonal, erfinv, invgammainc, quad2d_qgaus, &
                quad2d_romb, romb, gammainc_scal, gammainc
        use random, only : random_normal
        use hypergeo, only : hygfx
        use test, only : phi, ni, qi, nr1, qr, ttr, vx,vy, lf, &
            miupper, mrupper, mrthresh, ar, ai, br, bi, cr, ci, dr, di, &
            n0i, n0r, lambda0r, lambda0i, alphar, alphai, cw, f,&
            limit1, limit2, limit3, limit4, dintegral, calculate_parameters, &
            dintegral2, integral_m1, lam_freeze, n0_freeze
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
        real(wp) :: x1,x2,eps,h1,hmin,ff1,aa,bb!, f
        real(wp), dimension(10000001) :: b,x,r
        real(wp), dimension(10000000) :: a,c
        integer(i4b), allocatable, dimension(:) :: seed
        integer(i4b) :: l
        real(wp) :: rr
        real(wp) :: ss, p

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
        x1=0._wp
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
        
        
        ! now test double integration - mode 2 fragmentation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Some inputs for model                                                   !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        phi=0.5_wp
        T=265._wp
        f=-cw*(t-ttr)/lf ! fraction frozen at the end of stage 1 of freezing
        nr1=2.e3_wp
        qr=2.e-5_wp
        ni=.001e3_wp
        qi=0.1e-3_wp
        mrthresh = cr*150.e-6_wp**dr    ! mass of rain with D>150 microns
        mrupper=cr*1.e-3_wp**dr         ! mass of rain with D>1 mm
        miupper=ci*5.e-3_wp**di         ! mass of ice with D>5 mm
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Calculate some things prior                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call calculate_parameters(nr1,qr,cr,dr,alphar, n0r, lambda0r)
        call calculate_parameters(ni,qi,ci,di,alphai, n0i, lambda0i)
        vx=ar*gamma(br+alphar+1.0_wp)/gamma(alphar+1.0_wp)/lambda0r**(br)
        vy=ai*gamma(bi+alphai+1.0_wp)/gamma(alphai+1.0_wp)/lambda0i**(bi)
    
        p=0.9999_wp
        p=invgammainc(p,alphar+1.0_wp)
!         p=14.938751954612494_wp
!         print *,p,gammainc(alphar+1.0,p), alphar
        mrupper=cr*(p/lambda0r)**dr     ! mass of rain that is larger than 99.99% of dist
        miupper=ci*(p/lambda0i)**di     ! mass of ice that is larger than 99.99% of dist
        mrupper=min(mrupper,miupper)    ! min of two
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
        x1=mrthresh
        x2=mrupper
    
    !     print *,n0r,lambda0r,n0i,lambda0i, &
    !         n0r*gamma(1.+alphar)/lambda0r**(1.+alphar), &
    !         cr*n0r*gamma(1.+alphar+dr)/lambda0r**(1.+alphar+dr), &
    !         n0i*gamma(1.+alphai)/lambda0i**(1.+alphai), &
    !         ci*n0i*gamma(1.+alphai+di)/lambda0i**(1.+alphai+di)
        ss=0._wp
    !     print *,cr*(p/lambda0r)**dr,ci*(p/lambda0i)**di,mrthresh
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! for mode-2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(mrupper.gt.mrthresh) then
            call quad2d_romb(dintegral,limit1,limit2,x1,x2,ss)
        endif
        
        print *,'Double Integral solution is: ', ss, miupper,gammainc(2.0_wp,0.99_wp)        
        ss=0._wp
        if(mrupper.gt.mrthresh) then
            call quad2d_qgaus(dintegral,limit1,limit2,x1,x2,ss)
        endif
        
        print *,'Double Integral solution 2 is: ', ss   
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! for mode-1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        mrthresh = cr*1.e-6_wp**dr    ! mass of rain with D>1 microns
        x1=mrthresh ! can happen for all sizes, so pick a small one
        x2=mrupper
        if(mrupper.gt.mrthresh) then
            call quad2d_romb(dintegral2,limit3,limit4,x1,x2,ss)
        endif
        
        print *,'Double Integral solution is: ', ss, miupper,gammainc(2.0_wp,0.99_wp)        
        ss=0._wp
        if(mrupper.gt.mrthresh) then
            call quad2d_qgaus(dintegral2,limit3,limit4,x1,x2,ss)
        endif
        
        print *,'Double Integral solution 2 is: ', ss   
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! for mode-1 - freezing, no collisions
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        mrthresh = cr*1.e-6_wp**dr    ! mass of rain with D>1 microns
        x1=mrthresh ! can happen for all sizes, so pick a small one
        x2=mrupper
        lam_freeze=lambda0r
        n0_freeze=n0r
        if(mrupper.gt.mrthresh) then
            ss=romb(integral_m1,x1,x2)
        endif
        
        print *,'Integral solution is: ', ss,x1,x2      
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        
        
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
    


