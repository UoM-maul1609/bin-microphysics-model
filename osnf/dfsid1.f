      !>@author
      !>Paul Connolly, The University of Manchester
      !>@brief
      !>code to solve a calculate the 1st derivative of a function
      !> richardson extrapolation: 
      !> see https://www.researchgate.net/profile/Toshio_Fukushima/publication/320008726_Efficient_numerical_differentiation_by_Richardson_extrapolation_applied_to_forward_difference_formulas/links/59c73513a6fdccc7191edae2/Efficient-numerical-differentiation-by-Richardson-extrapolation-applied-to-forward-difference-formulas.pdf
      !>@param[in] x,h0,delta
      !>@param[inout] err
      !>@return dfsid1: gradient
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
        integer JMAX1; 
        real(wp) :: BETA,ERRMIN,LOG2,fx
        parameter (JMAX1=8,BETA=8.e0_wp,ERRMIN=1.e-35_wp)
        parameter (LOG2=0.69314718055994530941723212145818e0_wp)
        integer j,kmin,k; 
        real(wp),dimension(JMAX1,JMAX1) :: T
        real(wp) :: hj,errM,errP,factor,errT,errX 
        
        fx=func(x)
        hj=2.e0_wp**(floor(log(abs(h0))/LOG2))
        if(h0.lt.0.e0_wp) then
            hj=-hj 
        endif
        T(1,1)=(func(x+hj)-fx)/hj; dfsid1=T(1,1) 
        errM=0.5e0_wp*delta*abs(dfsid1); errP=1.e-38_wp; err=1.e38_wp !e-66 and e99
        do j=2,JMAX1
            hj=hj/BETA; factor=BETA; kmin=1; T(1,j)=(func(x+hj)-fx)/hj 
            do k=2,j
                T(k,j)=T(k-1,j)+(T(k-1,j)-T(k-1,j-1))/(factor-1.e0_wp) 
                factor=BETA*factor 
                errT=max(abs(T(k,j)-T(k-1,j)),abs(T(k,j)-T(k-1,j-1))) 
                if(errT.le.err) then
                    kmin=k; err=errT; dfsid1=T(k,j) 
                endif
            enddo 
            errX=err*err/errP 
            if(errX.le.errM) then
                err=max(ERRMIN,errX); return
            endif
            if(err.le.errM) return 
            if(kmin.eq.1) then
                err=errX; return 
            endif
            errP=err 
        enddo
        err=max(ERRMIN,errX) 
        return; 
      end function dfsid1