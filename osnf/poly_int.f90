	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>code to do a linear interpolation based on neville's algorithm
	!> https://en.wikipedia.org/wiki/Neville%27s_algorithm
	!>param[in] xarr,yarr,x: arrays and position
	!>param[inout] y,dy: value of y and error in y
	subroutine poly_int(xarr,yarr,x,y,dy)
	    use numerics_type
        use numerics, only : assert_eq,iminloc,numerics_error
        implicit none
        real(wp), dimension(:), intent(in) :: xarr,yarr
        real(wp), intent(in) :: x
        real(wp), intent(out) :: y,dy
        integer(i4b) :: m,n,ns
        real(wp), dimension(size(xarr)) :: c,d,factor,dx
        n=assert_eq(size(xarr),size(yarr),'poly_int')
        c=yarr
        d=yarr
        dx=xarr-x
        ns=iminloc(abs(dx))
        y=yarr(ns) ! nearest neighbour
        ns=ns-1
        do m=1,n-1 ! there are n-1 additional terms
            factor(1:n-m)=dx(1:n-m)-dx(1+m:n)
            if (any(factor(1:n-m) == 0.0)) &
                call numerics_error('poly_int: calculation failure')
            factor(1:n-m)=(c(2:n-m+1)-d(1:n-m))/factor(1:n-m) ! 3.1.3 recurrence relation
            d(1:n-m)=dx(1+m:n)*factor(1:n-m) ! 3.1.5 formulae
            c(1:n-m)=dx(1:n-m)*factor(1:n-m)
            if (n-m > 2*ns) then
                dy=c(ns+1)
            else
                dy=d(ns)
                ns=ns-1
            end if
            y=y+dy
        end do
	end subroutine poly_int
