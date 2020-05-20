	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>code to find position in an array
	!>param[in] xarr: array of xs 
	!>param[inout] x: value where position wanted
	!>@return find_pos: integer pointing to the position in the array
	function find_pos(xarr,x)
	    use numerics_type
        implicit none
        real(wp), dimension(:), intent(in) :: xarr
        real(wp), intent(in) :: x
        integer(i4b) :: find_pos
        integer(i4b) :: n,ilow,imid,iupp
        logical(lgt) :: ascnd
        n=size(xarr)
        ascnd = (xarr(n) >= xarr(1))
        ilow=0
        iupp=n+1
        do
            if (iupp-ilow <= 1) exit
            imid=(iupp+ilow)/2 ! just bisects array until it finds position 
            if (ascnd .eqv. (x >= xarr(imid))) then
                ilow=imid
            else
                iupp=imid
            end if
        end do
        if (x == xarr(1)) then
            find_pos=1
        else if (x == xarr(n)) then
            find_pos=n-1
        else
            find_pos=ilow
        end if
	end function find_pos
