    !>@author
    !>Paul Connolly, The University of Manchester
    !>@brief
    !>code to solve a tridiagonal system of equations based 
    !> based on the Thomas algorithm: 
    !> https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    !>@param[in] a,b,c,r
    !>@param[inout] x - solution vector
    subroutine tridiagonal(a,b,c,r,x)
        use numerics_type
        use numerics, only : assert_eq, numerics_error
        implicit none
        real(wp), dimension(:), intent(in) :: a,b,c,r
        real(wp), dimension(:), intent(inout) :: x
        real(wp), dimension(size(b)) :: cp,rp
        real(wp) :: m
        integer(i4b) :: n,i
        
        n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(x)/),'tridiagonal')
        m=b(1)
        if (m == 0.0) call numerics_error('tridiagonal: Error first divide')
        x(1)=r(1)/m

        ! initialize c-prime and r-prime
        cp(1) = c(1)/b(1)
        rp(1) = x(1)
        ! solve for vectors c-prime and r-prime
        do i = 2,n
            m = b(i)-cp(i-1)*a(i-1)
            if (m == 0.0) call numerics_error('tridiagonal: Error 2nd divide')
            cp(i) = c(i-1)/m
            rp(i) = (r(i)-rp(i-1)*a(i-1))/m
        enddo
        ! initialize x
        x(n) = rp(n)
        ! solve for x from the vectors c-prime and r-prime
        do i = n-1, 1, -1
          x(i) = rp(i)-cp(i)*x(i+1)
        end do        
    end subroutine tridiagonal