module numerics_type
    
	integer, parameter :: lgt = kind(.true.)
    integer, parameter :: i4b = selected_int_kind(9)
    integer, parameter :: i2b = selected_int_kind(4)
    integer, parameter :: i1b = selected_int_kind(2)
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)
#if VAR_TYPE == 1
    integer, parameter :: wp = kind(1.0d0)
	integer, parameter :: wpc = kind((1.0d0,1.0d0))
#endif
#if VAR_TYPE == 0
    integer, parameter :: wp = kind(1.0)
	integer, parameter :: wpc = kind((1.0,1.0))
#endif
	integer, parameter :: spc = kind((1.0,1.0))
	integer, parameter :: dpc = kind((1.0d0,1.0d0))
    real(wp), parameter :: pi=3.141592653589793238462643383279502884197_wp
    real(wp), parameter :: twopi=2.*wp*pi
end module numerics_type
