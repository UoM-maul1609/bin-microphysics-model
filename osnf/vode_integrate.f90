	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>module variables for integration using vode
    module ode_vars
        use numerics_type
        implicit none
        logical(lgt), save :: initialise_ode=.true.
        integer(i4b) :: neq,itol=2,itask=1,istate=1,iopt=1,lrw,liw,mf=22,iout
        !real(wp), dimension(:), allocatable :: rpar
        !integer(i4b), dimension(:), allocatable :: ipar
        real(wp), dimension(1) :: rpar
        integer(i4b), dimension(1) :: ipar
        
        
        
        contains
        
        subroutine jac (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
            use numerics_type
            implicit none
            real(wp) :: t
            real(wp), dimension(neq) :: y
            real(wp), dimension(nrpd, neq) :: pd
            integer(i4b) :: neq, ml, mu, nrpd
            real(wp) :: rpar
            integer(i4b) :: ipar

        end subroutine jac    
    end module ode_vars

	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>code to perform integration using vode
	!>@ref dvode
	!>param[inout] ystart: solution vector
	!>param[in] x1,x2,eps,h1,hmin
	!>param[in] fdash: function to calculate derivatives
	subroutine vode_integrate(ystart,x1,x2,eps,h1,hmin,fdash)
	    use numerics_type
	    use ode_vars
	    use numerics, only : dvode,xsetf
	    implicit none
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
        
	    real(wp) :: x
	    real(wp), dimension(size(ystart)) :: rtol,atol,y
	    real(wp), dimension(22+9*size(ystart)+2*size(ystart)**2) :: rwork
	    integer(i4b), dimension(30+size(ystart)) :: iwork

        

	    if(initialise_ode) then
	        neq=size(ystart)
            lrw=22+9*neq+2*neq**2
            liw=30+neq
            atol=1.e-10_wp
            rtol=eps
            x=x1
            y=ystart
            iwork=0
            rwork=0._wp
            call xsetf(0) ! do not print messages
            iwork(6)=100 ! max steps
            iwork(5)=5 ! order
            rwork(5)=h1 ! initial step
            rwork(7)=hmin ! minimum step
            rwork(14)=2._wp ! tolerance scale-factor
	    endif
	    
	    if(x1.lt.x2) then
            do while (x .lt. x2)
                istate=1
            	!>@callgraph 
            	!>@code call dvode(fex,neq,y,x,x2,itol,rtol,atol,itask,istate,iopt, & 
                !>            rwork,lrw,iwork,liw,jac,mf,rpar,ipar)
                !>@endcode
                call dvode(fex,neq,y,x,x2,itol,rtol,atol,itask,istate,iopt, & 
                            rwork,lrw,iwork,liw,jac,mf,rpar,ipar)
            enddo
        else
            do while (x .gt. x2)
                istate=1
                call dvode(fex,neq,y,x,x2,itol,rtol,atol,itask,istate,iopt, & 
                            rwork,lrw,iwork,liw,jac,mf,rpar,ipar)   
            enddo        
        endif
        ystart=y
        
        contains
            subroutine fex(neq,tt,y,ydot,rpar,ipar)
                integer(i4b), intent(in) :: neq
                real(wp), intent(inout) :: tt
                real(wp), intent(in) :: rpar
                integer(i4b), intent(in) :: ipar
                real(wp), dimension(neq), intent(inout) :: y,ydot
            
                call fdash(tt,y,ydot)
            end subroutine fex
        
	
	end subroutine vode_integrate
