    !>@author
    !>Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
    !>Reference:
    !>    Shanjie Zhang, Jianming Jin,
    !>    Computation of Special Functions,
    !>    Wiley, 1996,
    !>    ISBN: 0-471-11963-6,
    !>    LC: QA351.C45.
    !>@copyright Public Domain
    !>@brief
    !> module for 
    !> calculating hypergeometric functions (arbitrary precision)
    !> Changes by Paul Connolly, 2019
	module hypergeo
		private
		public :: hygfx
		contains
		
        SUBROUTINE HYGFX(A,B,C,X,HF)
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI for computing psi function
C       ====================================================
C
        use numerics_type, only : wp
        IMPLICIT REAL(WP) (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793_wp
        EL=.5772156649015329_wp
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0_wp-X.LT.1.0e-15_wp.AND.C-A-B.LE.0.0_wp
        L2=A.EQ.INT(A).AND.A.LT.0.0_wp
        L3=B.EQ.INT(B).AND.B.LT.0.0_wp
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0_wp
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0_wp
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0e-15_wp
        IF (X.GT.0.95) EPS=1.0e-8_wp
        IF (X.EQ.0.0_wp.OR.A.EQ.0.0_wp.OR.B.EQ.0.0_wp) THEN
           HF=1.0_wp
           RETURN
        ELSE IF (1.0_wp-X.EQ.EPS.AND.C-A-B.GT.0.0_wp) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0_wp+X.LE.EPS.AND.ABS(C-A+B-1.0_wp).LE.EPS) THEN
           G0=sqrt(PI)*2.0_wp**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0_wp+A/2.0_wp-B,G2)
           CALL GAMMA(0.5_wp+0.5_wp*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0_wp
           R=1.0_wp
           DO 10 K=1,NM
              R=R*(A+K-1.0_wp)*(B+K-1.0_wp)/(K*(C+K-1.0_wp))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0_wp
           R=1.0_wp
           DO 15 K=1,NM
              R=R*(C-A+K-1.0_wp)*(C-B+K-1.0_wp)/(K*(C+K-1.0_wp))*X
15            HF=HF+R
           HF=(1.0_wp-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0_wp) THEN
           X=X/(X-1.0_wp)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0_wp) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75_wp) THEN
           GM=0.0_wp
           IF (abs(C-A-B-INT(C-A-B)).LT.1.0e-15_wp) THEN
              M=INT(C-A-B)
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A+M,GAM)
              CALL GAMMA(B+M,GBM)
              CALL PSI(A,PA)
              CALL PSI(B,PB)
              IF (M.NE.0) GM=1.0_wp
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0_wp
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0_wp
              R0=1.0_wp
              R1=1.0_wp
              SP0=0._wp
              SP=0.0_wp
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0_wp)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0_wp)*(B+K-1.0_wp)/(K*(K-M))*(1.0_wp-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0_wp/(A+K-1.0_wp)+1.0_wp/(B+K-1.0_wp)
     &                  -1.0_wp/K
                 F1=PA+PB+SP0+2.0_wp*EL+log(1.0_wp-X)
                 DO 55 K=1,250
                    SP=SP+(1.0_wp-A)/(K*(A+K-1.0_wp))+(1.0_wp-B)/
     &                    (K*(B+K-1.0_wp))
                    SM=0.0_wp
                    DO 50 J=1,M
50                     SM=SM+(1.0_wp-A)/((J+K)*(A+J+K-1.0_wp))+1.0_wp/
     &                    (B+J+K-1.0_wp)
                    RP=PA+PB+2.0_wp*EL+SP+SM+log(1.0_wp-X)
                    R1=R1*(A+M+K-1.0_wp)*(B+M+K-1.0_wp)/
     &                    (K*(M+K))*(1.0_wp-X)
                    F1=F1+R1*RP
                    IF (abs(F1-HW).LT.abs(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0_wp-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0_wp)*(B-M+K-1.0_wp)/
     &                    (K*(K-M))*(1.0_wp-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0_wp/K
                 F1=PA+PB-SP0+2.0_wp*EL+log(1.0_wp-X)
                 DO 80 K=1,250
                    SP=SP+(1.0_wp-A)/(K*(A+K-1.0_wp))+(1.0_wp-B)/
     &                    (K*(B+K-1.0_wp))
                    SM=0.0_wp
                    DO 75 J=1,M
75                     SM=SM+1.0_wp/(J+K)
                    RP=PA+PB+2.0_wp*EL+SP-SM+log(1.0_wp-X)
                    R1=R1*(A+K-1.0_wp)*(B+K-1.0_wp)/(K*(M+K))*(1.0_wp-X)
                    F1=F1+R1*RP
                    IF (ABS(F1-HW).LT.ABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              CALL GAMMA(C-A-B,GCAB)
              CALL GAMMA(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0_wp-X)**(C-A-B)
              HF=0.0_wp
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0_wp)*(B+K-1.0_wp)/
     &                 (K*(A+B-C+K))*(1.0_wp-X)
                 R1=R1*(C-A+K-1.0_wp)*(C-B+K-1.0_wp)/(K*(C-A-B+K))
     &              *(1.0_wp-X)
                 HF=HF+R0+R1
                 IF (ABS(HF-HW).LT.ABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0_wp
           IF (C.GT.A.AND.C.LT.2.0_wp*A.AND.
     &         C.GT.B.AND.C.LT.2.0_wp*B) THEN
              A0=(1.0_wp-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0_wp
           R=1.0_wp
           DO 100 K=1,250
              R=R*(A+K-1.0_wp)*(B+K-1.0_wp)/(K*(C+K-1.0_wp))*X
              HF=HF+R
              IF (ABS(HF-HW).LE.ABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0_wp) THEN
           X=X1
           C0=1.0_wp/(1.0_wp-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function ג(x)
C       Input :  x  --- Argument of ג(x)
C                       ( x is not equal to 0,-1,-2,תתת)
C       Output:  GA --- ג(x)
C       ==================================================
C
        use numerics_type, only : wp
        IMPLICIT REAL(WP) (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793_wp
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0_wp) THEN
              GA=1.0_wp
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0_wp+300
           ENDIF
        ELSE
           IF (abs(X).GT.1.0_wp) THEN
              Z=ABS(X)
              M=INT(Z)
              R=1.0_wp
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0_wp,0.5772156649015329_wp,
     &          -0.6558780715202538_wp, -0.420026350340952e-1_wp,
     &          0.1665386113822915_wp,-.421977345555443e-1_wp,
     &          -.96219715278770e-2_wp, .72189432466630e-2_wp,
     &          -.11651675918591e-2_wp, -.2152416741149e-3_wp,
     &          .1280502823882e-3_wp, -.201348547807e-4_wp,
     &          -.12504934821e-5_wp, .11330272320e-5_wp,
     &          -.2056338417e-6_wp, .61160950e-8_wp,
     &          .50020075e-8_wp, -.11812746e-8_wp,
     &          .1043427e-9_wp, .77823e-11_wp,
     &          -.36968e-11_wp, .51e-12_wp,
     &          -.206e-13_wp, -.54e-14_wp, .14e-14_wp, .1e-15_wp/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0_wp/(GR*Z)
           IF (ABS(X).GT.1.0_wp) THEN
              GA=GA*R
              IF (X.LT.0.0_wp) GA=-PI/(X*GA*SIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        SUBROUTINE PSI(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        use numerics_type, only : wp
        IMPLICIT REAL(WP) (A-H,O-Z)
        XA=ABS(X)
        PI=3.141592653589793_wp
        EL=.5772156649015329_wp
        S=0.0_wp
        IF (X.EQ.INT(X).AND.X.LE.0.0_wp) THEN
           PS=1.0_wp+300_wp
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0_wp/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5_wp)) THEN
           N=XA-.5_wp
           DO 20 K=1,N
20            S=S+1.0_wp/(2.0_wp*K-1.0_wp)
           PS=-EL+2.0_wp*S-1.386294361119891_wp
        ELSE
           IF (XA.LT.10.0_wp) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0_wp/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0_wp/(XA*XA)
           A1=-.8333333333333e-01_wp
           A2=.83333333333333333e-02_wp
           A3=-.39682539682539683e-02_wp
           A4=.41666666666666667e-02_wp
           A5=-.75757575757575758e-02_wp
           A6=.21092796092796093e-01_wp
           A7=-.83333333333333333e-01_wp
           A8=.4432598039215686_wp
           PS=LOG(XA)-.5_wp/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0_wp) PS=PS-PI*cos(PI*X)/SIN(PI*X)-1.0_wp/X
        RETURN
        END
	end module hypergeo
