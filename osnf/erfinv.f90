	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@copyright 2019
	!>@brief
	!>calculate the inverse of erf
	!> This implementation is based on the rational approximation
    !> of percentage points of normal distribution available from
    !> https://www.jstor.org/stable/2347330.
	!>@param[in] x: input integer
	!>@param[inout] ans: the output
	subroutine erfinv(x,ans) 
	    use, intrinsic :: IEEE_ARITHMETIC
	    use numerics, only : LN2,A0,A1,A2,A3,A4,A5,A6,A7, &
	                    B0,B1,B2,B3,B4,B5,B6,B7, &
	                    C0,C1,C2,C3,C4,C5,C6,C7, &
	                    D0,D1,D2,D3,D4,D5,D6,D7, &
	                    E0,E1,E2,E3,E4,E5,E6,E7, &
	                    F0,F1,F2,F3,F4,F5,F6,F7
	    use numerics_type
        implicit none
        real(wp), intent(in) :: x
        real(wp), intent(inout) :: ans
        complex(wp) :: num, den
        real(wp) :: abs_x, r
        
        if((x>1._wp) .or. (x<-1._wp)) then
            ans=IEEE_VALUE(0._wp,IEEE_QUIET_NAN)
            return
        elseif((x==1._wp) ) then
            ans=IEEE_VALUE(0._wp,IEEE_POSITIVE_INF)
            return
        elseif((x==-1._wp) ) then
            ans=IEEE_VALUE(0._wp,IEEE_POSITIVE_INF)
        endif
        
        
        abs_x=abs(x)
        if(abs_x <= 0.85_wp) then
            r =  0.180625_wp - 0.25_wp * x * x
            num = (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) &
                * r + A2) * r + A1) * r + A0)
            den = (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) &
                * r + B2) * r + B1) * r + B0)
            ans = x * num / den
            return
        endif
        
        r = sqrt(LN2 - log(1.0_wp - abs_x))
        if (r <= 5.0_wp) then
            r = r - 1.6;
            num = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * &
                r + C2) * r + C1) * r + C0)
            den = (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * &
                r + D2) * r + D1) * r + D0)
        else 
            r = r - 5.0_wp
            num = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * &
                r + E2) * r + E1) * r + E0)
            den = (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * &
                r + F2) * r + F1) * r + F0)
        endif
        if (x < 0._wp) then
            ans= real(-num/den,wp)
            return
        else 
            ans = real(num/den,wp)
            return
        endif    
        
    end subroutine erfinv
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


