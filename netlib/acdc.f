C.... *****************************************************************
C.... Beginning Of Sample Problem Description.
C.... *****************************************************************

C.... *****************************************************************

C.... This Is A Sample Main Program For Acdc, An Automatic 
C.... Continuation Code For Solving `Stiff' Boundary Value Problems.
C.... The Code Is Fully Documented In The User's Guide, Available As 
C.... A Latex File From Netlib. Machine Precision Is Provided 
C.... Automatically Via Netlib (In The Routine D1mach). The Routine 
C.... D1mach May Need To Be Altered By The User Depending On Which
C.... Machine Is Being Used.

C.... Any Questions, Reports Of Bugs, Or Comments (Good Or Bad) Will 
C.... Gladly Be Received By Ross Wright Or Jeff Cash, Department Of 
C.... Mathematics, Imperial College, London Sw7 2bz, United Kingdom. 
C.... [ E-Mail r.wright@ic.ac.uk Or j.cash@ic.ac.uk ]

C.... *****************************************************************

      Program Runacdc
      Implicit Double Precision(A-H,O-Z)
      Parameter(Lwrkfl = 500000)
      Parameter(Lwrkin = 50000)
      Dimension Wrk(Lwrkfl),Iwrk(Lwrkin)
      Dimension U(2,12400),Tol(2),Ltol(2),Xx(12400),Fixpnt(1)
      Logical Linear, Giveu, Givmsh, Giveps
      External Fsub, Dfsub, Gsub, Dgsub
      Common /Prob/ Linornon
      Common /Valpi/ Pi

      Pi = 4.0d0*Atan(1.0d0)

C.... *****************************************************************

C.... This Driver Is Set Up To Run Two Problems, One Linear And One
C.... Nonlinear.

C.... The Linear Test Problem (Which Has Boundary Layers At X = -1
C.... And At X = 1) Is 

C.... Eps*Y'' - Y = -(Eps*Pi*Pi + 1)*Cos(Pi*X),   Y(-1) = Y(1) = 0.

C.... The Nonlinear Test Problem (Which Has Boundary Layers At X = 0
C.... And At X = 1) Is 

C.... Eps*Y'' + Y*Y' - Y = 0,   Y(0) = 1, Y(1) = -1/3.

C.... Acdc Solves First Order Systems Of Differential Equations, Thus
C.... Both Problems Are Converted To First Order Form In Subroutine
C.... Fsub Below.

C.... *****************************************************************

C.... Specify The Number Of Components And The Desired Final Eps Value.

      Ncomp = 2
      Epsmin = 1.d-7

C.... Ntol Is The Number Of Tolerances

      Ntol = 2
      Tol(1) = 1.d-8
      Tol(2) = 1.d-4
      Ltol(1) = 1
      Ltol(2) = 2

C.... The User May Choose A Linear Or Nonlinear Test Problem

 10   Write(*,1000) 
      Read*, Linornon

C.... Specify The Boundary Points

      If (Linornon .eq. 1) Then
        Aleft = 0.d0
        Aright = 1.0d0
      Else If (Linornon .eq. 0) Then
        Aleft = -1.0d0
        Aright = 1.0d0
      Else
         Write(*,2000)
         Goto 10
      Endif

      If (Linornon .eq. 1) Then
         Linear = .false.
      Else
         Linear = .true.
      Endif

C.... Initialise Variables 

      Nudim = 2
      Nlbc = 1
      Nfxpnt = 0
      Nucol = 12400
      Nmsh = 0

      Giveu = .false.
      Givmsh = .false.
      Giveps = .false.

      Call Acdc(Ncomp, Nlbc, Nucol, Aleft, Aright, Nfxpnt, Fixpnt,
     +            Ntol, Ltol, Tol, Linear, Givmsh, Giveu, Nmsh, Xx,
     +            Nudim, U, Nmax, Lwrkfl, Wrk, Lwrkin, Iwrk, Giveps,
     +            Eps, Epsmin, Fsub, Dfsub, Gsub, Dgsub, Iflbvp)

C-----------------------------------------------------------------------
 1000 Format(//,'Choose A Test Problem (0 = Linear, 1 = Nonlinear).')
 2000 Format(/, '** Only 0 or 1 is allowed.')
C-----------------------------------------------------------------------
      End

      Subroutine Fsub(Ncomp,X,Z,F,Eps)
      Implicit Double Precision(A-H,O-Z)
      Dimension Z(*),F(*)
      Common /Prob/ Linornon
      Common /Valpi/ Pi

      If (Linornon .eq. 0) Then
         Pix = Pi*X
         F(1) = Z(2)
         F(2) = (Z(1)-Eps*Pi*Pi*Cos(Pix)-Cos(Pix))/Eps
      Else
         F(1) = Z(2)
         F(2) = (Z(1)*(1.d0-Z(2)))/Eps
      Endif

      Return
      End 

      Subroutine Dfsub(Ncomp,X,Z,Df,Eps)
      Implicit Double Precision(A-H,O-Z)
      Dimension Z(*),Df(Ncomp,*)
      Common /Prob/ Linornon

      If (Linornon .eq. 0) Then
         Df(1,1) = 0.d0
         Df(1,2) = 1.0d0
         Df(2,1) = 1.0d0/Eps
         Df(2,2) = 0.0d0
      Else 
         Df(1,1) = 0.d0
         Df(1,2) = 1.0d0
         Df(2,1) = (1.d0-Z(2))/Eps
         Df(2,2) = -Z(1)/Eps
      Endif

      Return
      End 

      Subroutine Gsub(I,Ncomp,Z,G,Eps)
      Implicit Double Precision(A-H,O-Z)
      Dimension Z(*)
      Common /Prob/ Linornon

      If (Linornon .eq. 0) Then
         G = Z(1)
      Else 
         If (I .eq. 1) G = Z(1)-1.d0
         If (I .eq. 2) G = Z(1)+1.d0/3.d0
      Endif

      Return 
      End 

      Subroutine Dgsub(I,Ncomp,Z,Dg,Eps)
      Implicit Double Precision(A-H,O-Z)
      Dimension Z(*),Dg(*)

      Dg(1) = 1.d0
      Dg(2) = 0.d0

      Return
      End 

C.... *****************************************************************
C.... End Of Sample Problem Description.
C.... *****************************************************************

      Subroutine Acdc(Ncomp, Nlbc, Nucol, Aleft, Aright, Nfxpnt, Fixpnt,
     +            Ntol, Ltol, Tol, Linear, Givmsh, Giveu, Nmsh, Xx,
     +            Nudim, U, Nmax, Lwrkfl, Wrk, Lwrkin, Iwrk, Giveps,
     +            Eps, Epsmin, Fsub, Dfsub, Gsub, Dgsub, Iflbvp)

      Implicit Double Precision (A-H,O-Z)
      Dimension Fixpnt(*), Ltol(*), Tol(*)
      Dimension Xx(*), U(Nudim,*)
      Dimension Wrk(Lwrkfl), Iwrk(Lwrkin)
      Dimension Phi(3), E(3), Npr(2), Pmax(2), Hord(2)

      Logical Linear, Givmsh, Giveu, Giveps
      External Fsub, Dfsub, Gsub, Dgsub

      Common/Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0
      Common /Flags/ Ifinal,Iatt,Iback,Iprec
      Common /Convg/ Nits
      Common /Mshvar/ Hsml,Npr,Pmax,Hord

      Parameter ( Zero = 0.0d+0, One = 1.0d+0, Two = 2.0d+0 )
      Parameter ( Three = 3.0d+0, Third = 1.0d+0/3.0d+0, Huge = 1.d+30 )

C.... The First Part Of This Subroutine Is Used For Checking Input
C.... Parameters And Dividing Up The Floating Point And Integer
C.... Workspaces. The Rest Of This Subroutine Is Devoted To The 
C.... Continuation Algorithm.

C.... Output Details Of The Problem

      If (Iprint .ge. 0) Then
         If (Linear) Then
            Write (6,1000) Ncomp
         Else
            Write (6,1001) Ncomp
         Endif
         Write (6,1002) Aleft, Aright
         If (Nfxpnt .gt. 0) 
     +        Write (6,1003) Nfxpnt, (Fixpnt(I), I=1,Nfxpnt)
         Write (6,1004) (Ltol(Ip), Ip = 1,Ntol)
         Write (6,1005) (Tol(Ip), Ip = 1,Ntol)
      Endif

C.... Check For Invalid Input Parameters.  If Any Parameters Are
C.... Invalid, Exit With The Flag Iflbvp Set To -1.

      Iflbvp = -1
      If (Ncomp .le. 0)  Return
      If (Nlbc .lt. 0 .or. Nlbc .gt. Ncomp) Return
      If (Nucol .le. 0) Return
      If (Aleft .ge. Aright) Return
     
      If (Nfxpnt .lt. 0)  Return
      If (Maxcon .lt. 0)  Return
      If (.not. Givmsh) Nmsh = 0
      If (Givmsh .and. Nmsh .lt. Nfxpnt+2) Return
      If (Givmsh .and. Xx(1) .ne. Aleft) Return
      If (Givmsh .and. Xx(Nmsh) .ne. Aright) Return
      If (Nfxpnt .gt. 0) Then
         If (Fixpnt(1) .le. Aleft) Return
         If (Fixpnt(Nfxpnt) .ge. Aright) Return
         Do 10 I = 1, Nfxpnt-1
            If (Fixpnt(I+1) .le. Fixpnt(I)) Return
   10    Continue
      Endif

      If (Ntol .lt. 1) Return
      Do 20 I = 1, Ntol
         If (Ltol(I) .lt. 0 .or. Ltol(I) .gt. Ncomp) Return
         If (Tol(I) .le. Zero) Return
   20 Continue

      If (Giveu .and. .not. Givmsh) Return

      If (Nudim .le. 0) Return
      If (Lwrkfl .le. 0 .or. Lwrkin .le. 0) Return
      If (Itsaim .lt. 1) Return
      If (.not. Giveps) Eps = 0.5d0
      If (Eps .le. Zero .or. Eps .ge. One) Return
      If (Epsmin .gt. Eps .or. Epsmin .le. Zero) Return
      If (Iprint .ge. 0) Write (6,1006) Eps, Epsmin

C.... Calculate Maximum Number Of Mesh Points Possible With The
C.... Given Floating-Point And Integer Workspace.

      Isp = Lwrkfl - 2*Ntol - 13*Ncomp - 14*Ncomp*Ncomp

      If (Linear) Then
         Iden = 4*Ncomp*Ncomp + 10*Ncomp + 5
      Else
         Iden = 4*Ncomp*Ncomp + 11*Ncomp + 5
      Endif

      Nmax1 = Isp/Iden

      Isp = Lwrkin - 3*Ncomp
      Nmax2 = Isp/(Ncomp+2)

      Nmax = Min(Nmax1, Nmax2)
      If (Iprint .eq. 1) Write(6,1007) Nmax
      Nmax = Min(Nmax, Nucol)
      If (Iprint .eq. 1 .and. Nmax .eq. Nucol) Write(6,1008) Nmax

      If (Nmax .le. 1) Return

C.... Partition Floating Point Workspace.

      Irhs = 1
      Lrhs = Ncomp*Nmax

      Itpblk = Irhs + Lrhs
      Ltpblk = Ncomp*Ncomp

      Ibtblk = Itpblk + Ltpblk
      Lbtblk = Ncomp*Ncomp

      Iajac = Ibtblk + Lbtblk
      Lajac = 2*Ncomp*Ncomp*Nmax

      Ibhold = Iajac + Lajac
      Lbhold = Ncomp*Ncomp*Nmax

      Ichold = Ibhold + Lbhold
      Lchold = Ncomp*Ncomp*Nmax

      Ifval = Ichold + Lchold
      Lfval = Ncomp*Nmax

      Idef = Ifval + Lfval
      Ldef = Ncomp*(Nmax-1)

      Idef8 = Idef + Ldef
      Ldef8 = Ncomp*(Nmax-1)

      Iuold = Idef8 + Ldef8
      Luold = Ncomp*Nmax

      Itmrhs = Iuold + Luold
      Ltmrhs = Ncomp*Nmax

      Irhtri = Itmrhs + Ltmrhs
      Lrhtri = Ncomp*Nmax

      Idelu = Irhtri + Lrhtri
      Ldelu = Ncomp*Nmax

      Ixmer = Idelu + Ldelu
      Lxmer = Ncomp*Nmax

      Iutri = Ixmer + Lxmer
      Lutri = Ncomp*Nmax

      Iermx = Iutri + Lutri
      Lermx = Nmax

      Irtdc = Iermx + Lermx
      Lrtdc = Nmax

      Ixxold = Irtdc + Lrtdc
      Lxxold = Nmax
 
      Iuint = Ixxold + Lxxold
      Luint = Ncomp

      Iftmp = Iuint + Luint
      Lftmp = Ncomp

      Idgtm = Iftmp + Lftmp
      Ldgtm = Ncomp

      Idftm1 = Idgtm + Ldgtm
      Ldftm1 = Ncomp*Ncomp

      Idftm2 = Idftm1 + Ldftm1
      Ldftm2 = Ncomp*Ncomp

      Itmp = Idftm2 + Ldftm2
      Ltmp = Ncomp*12

      Idsq = Itmp + Ltmp
      Ldsq = Ncomp*Ncomp

      Ietst6 = Idsq + Ldsq
      Letst6 = Ntol

      Ietst8 = Ietst6 + Letst6
      Letst8 = Ntol

      Iextra = Ietst8 + Letst8
      Lextra = 9*Ncomp*Ncomp

      Iextra1 = Iextra + Lextra
      Lextra1 = Nmax

      Imsh = Iextra1 + Lextra1
      Lmsh = Nmax

      If (Linear) Then
         Ilast = Imsh + Lmsh
      Else
         Isol = Imsh + Lmsh
         Lsol = Ncomp*Nmax
         Ilast = Isol + Lsol
      Endif

C.... Partition Integer Workspace.
      
      Iiref = 1
      Liref = Nmax

      Iihcom = Iiref + Liref
      Lihcom = Nmax

      Iipvbk = Iihcom + Lihcom
      Lipvbk = Ncomp*Nmax

      Iipvlu = Iipvbk + Lipvbk
      Lipvlu = 3*Ncomp

C.... *****************************************************************
C.... Initialisation And Explanation Of Variables Used In The 
C.... Continuation Strategy Begins Here.
C.... *****************************************************************

C.... Initialise Extrapolation Variables

      Do 30 J = 1,3
        Phi(J) = Zero
        E(J) = Zero
 30   Continue

      Emin = One/Epsmin
      Ep = One/Eps

C.... The Maximum Monitor Function Values On The First Two Meshes
C.... For A Given Step Are Stored In Pmax(1) And Pmax(2).

      Pmax(1) = Zero
      Pmax(2) = Zero

C.... Nits Counts The Number Of Newton Iterations Required For 
C.... Convergence On The First Mesh For A Given Continuation Problem

      Nits = 0

C.... Ifinal Is A Flag Indicating Whether The Final Problem Is Being
C.... Attempted. If Iback = 1 Then We Do Not Backtrack. Iprec Is A 
C.... Flag That Is Set To One When We Believe That The A Calculation
C.... Is Being Attempted Which May Be Beyond The Bounds Imposed By
C.... The Precision Of The Machine.
C.... Epsp Is The Value Of Eps Beyond Which The Machine Precision Is
C.... (Possibly) Not Sufficient To Solve The Given Problem.

      Ifinal = 0
      Iback = 0
      Iprec = 0
      Epsp = Zero

C.... Istep Is A Flag Which Limits The Size Of The Continuation Steps
C.... After Experiencing A Failure For Some Step. 

      Istep = 0

C.... Nc Counts The Total Number Of Continuation Steps Taken.
C.... Nss Counts The Number Of Successful Continuation Steps Taken.

      Nc = 0
      Nss = 0

C.... Idc Counts The Number Of Consecutive Steps For Which The Maximum
C.... Value Of The Monitor Function Is Decreasing.
C.... Iextrap Is A Flag That Indicates Whether Monitor Function
C.... Extrapolation Is Possible.
C.... Istuk Is A Flag Which Is Set To One When The Continuation 
C.... Algorithm Is Backtracking.

      Idc = 0
      Iextrap = 0
      Istuk = 0

C.... Ilin Is A Flag Indicating Whether Linear Or Quadratic
C.... Extrapolation Should Be Used.

      Ilin = -1

C.... Hrat Is A Variable Utilised When We Have A Moving Layer.

      Hrat = One

C.... The Following Four Variables Are Used In The Parameter 
C.... Selection Process.

      Amax = Huge
      Bmax = Zero
      Cmax = Huge
      Estep = Huge

C.... Calculate Phimax, Which Is The Largest Monitor Function Value 
C.... That We Believe Is Permissible In Double Precision For The Given
C.... Tolerance.

      Tolmax = One
      Phimax = One/D1mach(3)
      Do 40 J = 1,Ntol
        Tolmax = Min(Tol(J),Tolmax)
 40   Continue
      Phimax = Phimax*Tolmax 
      Phiaim = Zero
      Phialt = Zero
      Epold = Zero
      Hsmlold = Zero
      Htpv = Zero

C.... We Predict That For Epsilon = 1 The Maximum Value Of The 
C.... Monitor Function Is Minimised By Phi(3). We Will Use This Value
C.... When Performing Extrapolation.
 
      E(3) = One
      Phi(3) = One/(Aright-Aleft)      
  
C.... *****************************************************************
C.... End Of Initialisation Of Continuation Variables.
C.... *****************************************************************

C.... *****************************************************************
C.... The Algorithm Loops Back To Line 50 For Every Continuation Step.
C.... *****************************************************************
 
 50   Eps = Max(One/Ep,Epsmin)
      Ep = One/Eps
 
C.... Solve For Epsmin Exactly
 
      If (Eps .lt. (1.00001d0*Epsmin)) Ifinal = 1
 
C.... Increase The Continuation Step Counter Nc.
 
      Nc = Nc+1
      If (Iprint .ge. 0) Write(6,1010) Nc, Eps 

C.... If We Have Reached The Maximum Number Of Continuation Steps
C.... Then This Will Be The Last Problem We Attempt.

      If (Nc .eq. Maxcon) Then
         If (Iprint .ge. 0) Write(6,1011) Maxcon  
         Iback = 1
         Ifinal = 1
      Endif

C.... Iatt Is A Flag That Governs Whether Mesh Selection Is Performed 
C.... Or Control Is Returned To The Driver To Select A New Parameter.
 
      Iatt = -1
      Iflbvp = 0
 
C.... Attempt To Solve The Latest Continuation Problem.
 
      Call Bvpsol(Ncomp, Nmsh, Nlbc, Aleft, Aright, 
     *   Nfxpnt, Fixpnt, Ntol, Ltol, Tol, Nmax, Linear, 
     *   Giveu, Givmsh, Xx, Nudim, U, 
     *   Wrk(Idef), Wrk(Idelu), 
     *   Wrk(Irhs), Wrk(Ifval),
     *   Wrk(Itpblk), Wrk(Ibtblk), Wrk(Iajac), Wrk(Ibhold), 
     *   Wrk(Ichold), Iwrk(Iipvbk), Iwrk(Iipvlu),
     *   Wrk(Iuint), Wrk(Iftmp), Wrk(Itmrhs), 
     *   Wrk(Idftm1), Wrk(Idftm2), Wrk(Idgtm), 
     *   Wrk(Iutri), Wrk(Irhtri), Wrk(Ixmer), 
     *   Wrk(Ixxold), Wrk(Iuold),
     *   Wrk(Itmp), Wrk(Idsq), Wrk(Irtdc), 
     *   Wrk(Ietst6), Wrk(Ietst8), Wrk(Iermx), 
     *   Iwrk(Iihcom), Iwrk(Iiref), Wrk(Idef8),
     *   Wrk(Iextra), Wrk(Iextra1), 
     *   Fsub, Dfsub, Gsub, Dgsub, Iflbvp, Eps)
 
C.... *****************************************************************
C.... The Logic For A Successful Continuation Step Starts Here
C.... *****************************************************************
 
      If (Iflbvp .eq. 0) Then

         Nss = Nss + 1
 
C.... If The Problem Eps = Epsmin Has Been Solved Succesfully Then
C.... We May Finish. Ifinal = 1 When Eps = Epsmin.
 
         If (Ifinal .eq. 1) Then 
            If (Iprint .ge. 0) Then
               Write(6,1012) Eps
               Write(6,1018) ('U',Ltol(J), J=1,Ntol)
               Write(6,*)
               Jstep = Max(Nmsh/30,1)
               Do 55 I = 1, Nmsh-1, Jstep
                  Write(6,1019) I,Xx(I),(U(Ltol(J),I),J=1,Ntol)
 55            Continue
               Write(6,1019) Nmsh,Xx(Nmsh),(U(Ltol(J),Nmsh),J=1,Ntol)
               Write(6,1020)
            Endif
            Return
         Endif       
 
C.... When Bactracking The Program Should Find A Problem That It Can
C.... Eventually Solve. The Relationship Between The Last Two
C.... Successful Epsilon Values Is Used To Restrict Future Epsilon 
C.... Values. This Is The Purpose Of Estep. If Istep = 1, Then We
C.... Must Restrict Future Values Of Ep.
 
         If (Istuk .eq. 1) Then
            Istep = 1
            Estep = Min(Ep/E(3),Estep)
         Endif
         Istuk = 0
 
C.... Calculate The Best Approximation To The Maximum Value Of The 
C.... Monitor Function By Extrapolating The Maximum Value On The 
C.... First Two Meshes. The Best Approximation Is Stored In Phit.
 
         If (Iextrap .eq. 0) Then
            H1 = Hord(1)/Pmax(1)
            H2 = Hord(2)/Pmax(2)
            If (H1/H2 .gt. 1.1d0) Then
               C1h = (Pmax(2) - Pmax(1))/(H2-H1)
               Phit = Pmax(1) - C1h*H1
            Else
               Phit = Pmax(2)
            Endif
         Endif
         
C.... It Is Possible That We Are Attempting To Solve Problems 
C.... That Are Beyond The Bounds Imposed By The Machine Precision.
C.... The User Is Warned Of This Possibility.
 
         If (Iprec .eq. 0 .and. Phit .gt. Phimax) Then
            If (Iprint .ge. 0) Write(6,1013) Eps
            Epsp = Eps
            Iprec = 1
         Endif
         
C.... Save Details Of Last Problem Solved Successfully In Case We
C.... Need To Backtrack
 
         Nnewold = Nmsh
         Call Meshdet(Wrk(Imsh), Xx, Nmsh)
         If (.not. Linear) Call Soldet(Wrk(Isol), U, Ncomp, Nmsh, 
     +        Ncomp, Nudim)
         
C.... If Iextrap Equals 1 Then We Have Abandoned Monitor Function
C.... Extrapolation As Our Parameter Selection Basis.
 
         If (Iextrap .eq. 1) Then
            E(1) = E(2)
            E(2) = E(3)
            E(3) = Ep
            If (Istep .eq. 1) Then
               Ep = Estep*Ep
            Else
               Ep = Emin
            Endif
            Goto 70
         Endif
 
C.... If We Have A Decreasing Monitor Function Twice Consecutively 
C.... Then Extrapolation Will Not Work. The Flag Iextrap Equals 1
C.... When Extrapolation Does Not Work.
 
         If (Phit .le. Phi(3)) Then
            If (Idc .eq. 1 .or. Phit .lt. Phi(2)) Then
               Iextrap = 1
               E(1) = E(2)
               E(2) = E(3)
               E(3) = Ep
               If (Istep .eq. 1) Then
                  Ep = Estep*Ep
               Else
                  Ep = Emin
               Endif
               Goto 70
            Else
               Idc = 1
               E(3) = Ep
               Phi(3) = Phit
            Endif
         
C.... Otherwise Update Extrapolation Data
 
         Else
            Idc = 0
            Do 60 J = 1,2
               E(J) = E(J+1)
               Phi(J) = Phi(J+1)
 60         Continue
            E(3) = Ep
            Phi(3) = Phit
         Endif
         
C.... The Variable Irest Is A Flag Which Informs Us Whether Our
C.... Extrapolation Procedure Is Accurate Or Not. Irest = 1 Means
C.... Extrapolation Is Not Reliable And So We Should Restrict The
C.... Next Parameter Step.
 
         Irest = 0
 
C.... Check Errors In Extrapolation Procedure And Calculate Hmult
C.... Which Restricts The Size Of Parameter Steps If Extrapolation
C.... Is Working Poorly.
 
         If (Nss .eq. 2 .and. Phit .gt. Two*Phiaim) Then
            Irest = 1
            Hh = Ep - Epold
            Errp = Abs(Phit - Phiaim)
            Fx = Errp/Hh**2
            Hmax = Sqrt(0.2d0*Phiaim/Fx)
            Hmult = (Epold+Hmax)/Epold
         Endif
         
         If (Nss .gt. 2 .and. Phit .gt. Two*Phiaim) Then
            Irest = 1
            Hh = Ep - Epold
            Errp = Abs(Phit - Phiaim)
            Errp1 = Abs(Phit - Phialt)
            If (Ilin .eq. -1) Then
               Fx = Errp/Hh**2
               Hmax = Sqrt(0.2d0*Phiaim/Fx)
               Hmult = (Epold+Hmax)/Epold
               Fx = Errp1/Hh**3
               Hmax = (0.2d0*Phiaim/Fx)**(Third)
               Hmult1 = (Epold+Hmax)/Epold
            Else
               Fx = Errp/Hh**3
               Hmax = (0.2d0*Phiaim/Fx)**(Third)
               Hmult = (Epold+Hmax)/Epold
               Fx = Errp1/Hh**2
               Hmax = Sqrt(0.2d0*Phiaim/Fx)
               Hmult1 = (Epold+Hmax)/Epold
            Endif
         Endif 
         Epold = Ep
         
C.... If The Extrapolation Procedure Is Particularly Inaccurate Then 
C.... We Need To Restrict It Permanently By Setting Istep = 1.

         Pinacc = Three*Max(Pmax(1),Phiaim,Phialt)
         If (Nss .gt. 2 .and. Phit .gt. Pinacc) Then
            Istep = 1
            Pestep = Max(Hmult,Hmult1)
            Estep = Min(Estep,Pestep)
         Endif        
         
C.... Decide Whether Linear Or Quadratic Extrapolation Is Best.
 
         If (Nss .gt. 2) Then
            Per1 = Abs(Phiaim-Phit)
            Per2 = Abs(Phialt-Phit)
            If (Per2 .lt. Per1) Then
               Ilin = -Ilin 
               If (Phit .gt. Two*Phialt) Then
                  Hmult = Hmult1
               Else
                  Irest = 0
               Endif
            Endif
         Endif  
         
C.... If The Continuation Steps Are Becoming Very Small Then We Have 
C.... Reached Our Final Problem
 
         Dele = E(3)-E(2)
         If (Dele .lt. 0.01d0*E(3) .or. (Nmsh .gt. (3*Nmax/4))) Then
            Ep = E(3)
            Epsmin = One/Ep
            If (Iprint .ge. 0) Then
               If (Dele .lt. 0.01d0*E(3)) Then
                  Write(6,1021) Epsmin
               Else
                  Write(6,1022) Epsmin
               Endif
            Endif
            Emin = Ep
            Iback = 1
            Goto 70      
         Endif
         
C.... The Following Section Of Code Calculates The Desired Value 
C.... (Phiaim) That We Would Like The Maximum Value Of The Monitor 
C.... Function To Take.
  
         Itru = 0
         If (Phit .gt. Pmax(1) .and. Pmax(2) .ne. Phit) Itru = 1
         If (Itru .eq. 1) Amax = -Phit*Phit/(Two*C1h)
         Bmax = Phit*H1
         If (Npr(1) .gt. 3*Npr(2)/2) Then
            Ccm = Max(Bmax-Three,Bmax/1.5d0)
            Cmax = Min(Cmax,Ccm)
         Endif
         If (Itru .eq. 0) Amax = Max(1.5d0*Bmax,Bmax+Three)
         If (Nss .ge. 2) Hrat = Max(H1/Hsmlold,One)
         Hsmlold = Hsml
         Bbm = Max(1.5d0*Bmax,Bmax+Three)
         If (.not. Linear) Then
            Fmax = Dble(Itsaim)*Bmax/Dble(Nits)
            If (Fmax .gt. Bmax) Fmax = Max(Fmax,Bbm)
            If (Fmax .lt. Bmax) Fmax = Max(Bmax-Three,Bmax/1.5d0,Fmax)
         Endif
         Htot = Min(Amax,Bbm,Cmax)
         If (Linear) Fmax = Htot
         If (Nss .gt. 1 .and. Htot .ne. Cmax) Then
            If (Itru .eq. 0 .or. Amax .gt. Htpv) Htot = Max(Htot,Htpv)
            Htot = Min(Htot,Fmax,Cmax)
         Endif
         Htpv = Htot
         Phiaim = Htot/(Hsml*Hrat)
         Phiaim = Max(1.05d0*Phit,Phiaim)
         
C.... Quadratic And Linear Extrapolation To Find The Value Of Ep That
C.... Corresponds To Phiaim
         
         If (Nss .ne. 1) Then
            
            F1 = (Phi(2)-Phi(1))/(E(2)-E(1))
            F2 = (Phi(3)-Phi(2))/(E(3)-E(2))
            Ff2 = (F2-F1)/(E(3)-E(1))
            
            If (Ff2 .gt. 0) Then
               Qa = Ff2
               Qb = F1-(E(1)+E(2))*Ff2
               Qc = Phi(1)-E(1)*(F1-E(2)*Ff2)-Phiaim
               Qd = Qb**2-4.d0*Qa*Qc
               If (Ilin .eq. 1) Then
                  Ep = Min(Emin,(-Qb+Sqrt(Qd))/(Two*Qa))
                  If (Ep .eq. Emin) Phiaim = Phi(1)+
     +                 (Ep-E(1))*(F1+(Ep-E(2))*Ff2)
                  Phialt = Phi(2)+F2*(Ep-E(2))
               Else
                  Ep = Min(Emin,E(2)+(Phiaim-Phi(2))/F2)
                  If (Ep .eq. Emin) Phiaim = Phi(2)+F2*(Ep-E(2))
                  Phialt = Phi(1)+(Ep-E(1))*(F1+(Ep-E(2))*Ff2)
               Endif
            Else
               G1 = (E(2)-E(1))/(Phi(2)-Phi(1))
               G2 = (E(3)-E(2))/(Phi(3)-Phi(2))
               Gg2 = (G2-G1)/(Phi(3)-Phi(1))
               If (Ilin .eq. 1) Then 
                  Ep = Min(Emin,
     +                 E(1)+(Phiaim-Phi(1))*(G1+(Phiaim-Phi(2))*Gg2))
                  If (Ep .eq. Emin) Then
                     Qa = Gg2
                     Qb = G1-(Phi(1)+Phi(2))*Gg2
                     Qc = E(1)-Phi(1)*(G1-Phi(2)*Gg2)-Ep
                     Qd = Qb**2-4.d0*Qa*Qc
                     Phiaim = (-Qb+Sqrt(Qd))/(Two*Qa)
                  Endif 
                  Phialt = Phi(2)+F2*(Ep-E(2))
               Else
                  Ep = Min(Emin,E(2)+G2*(Phiaim-Phi(2)))
                  If (Ep .eq. Emin) Phiaim = Phi(2)+F2*(Ep-E(2))
                  Qa = Gg2
                  Qb = G1-(Phi(1)+Phi(2))*Gg2
                  Qc = E(1)-Phi(1)*(G1-Phi(2)*Gg2)-Ep
                  Qd = Qb**2-4.d0*Qa*Qc
                  Phialt = (-Qb+Sqrt(Qd))/(Two*Qa)
               Endif 
            Endif
            
C.... Linear Extrapolation
            
         Else
            F2 = (Phi(3)-Phi(2))/(E(3)-E(2))
            Ep = Min(Emin,E(2)+(Phiaim-Phi(2))/F2)
         Endif
         
C.... Extrapolation May Not Be Working Very Well - If Not, Adjust
C.... Continuation Stepsize And Recalculate Phiaim.
         
         If ((Irest .eq. 1 .or. Istep .eq. 1) .and. Nss .ne. 1) Then
            If (Irest .eq. 1 .and. Istep .eq. 1) Then
               Ealt = Min(Epold*Hmult,Epold*Estep)
            Elseif (Irest .eq. 1) Then
               Ealt = Epold*Hmult
            Else
               Ealt = Epold*Estep
            Endif
            If (Ep .gt. Ealt) Then
               Ep = Min(Emin,Ealt)
               Philin = Phi(2)+F2*(Ep-E(2))
               If (Ff2 .gt. 0) Then
                  Phiqua = Phi(1)+(Ep-E(1))*(F1+(Ep-E(2))*Ff2)
               Else
                  Qa = Gg2
                  Qb = G1-(Phi(1)+Phi(2))*Gg2
                  Qc = E(1)-Phi(1)*(G1-Phi(2)*Gg2)-Ep
                  Qd = Qb**2-4.d0*Qa*Qc
                  Phiqua = (-Qb+Sqrt(Qd))/(Two*Qa)
               Endif
               If (Ilin .eq. 1) Then
                  Phiaim = Phiqua
                  Phialt = Philin
               Else
                  Phiaim = Philin
                  Phialt = Phiqua
               Endif
            Endif
         Endif
         
 70      Continue
         
C.... *****************************************************************
C.... End Of Logic For A Successful Continuation Step And 
C.... Beginning Of Logic For An Unsuccessful Continuation Step.
C.... *****************************************************************
 
      Else
 
C.... If A Problem Is Not Solved Successfully Then Print Details And
C.... Select An Alternative Value For The Continuation Parameter. 
C.... This Process Is Known As Backtracking. Istuk = 1 Indicates That 
C.... We Are Backtracking.
 
         Istuk = 1
 
C.... If Iback = 1 Then The Final Problem Has Not Been Solved.
C.... In This Case We Stop.
 
         If (Iback .eq. 1) Then
            If (Iprint .ge. 0) Then
               If (Epsp .ne. Zero) Then
                  Write(6,1014) Eps,Iflbvp,Eps,Epsp
               Else
                  Write(6,1015) Eps,Iflbvp,Eps
               Endif
            Endif
            Return
         Endif
         
C.... If Iprec = 2, Then We Know That We Cannot Define A Mesh On
C.... Which The Current Eps Value Can Be Solved To The Requested
C.... Tolerances. We Alter Epsmin Accordingly.

         If (Iprec .eq. 2) Then
            Epsmin = One/(Max((Ep+E(3))/Two,0.9d0*Ep))
         Endif

C.... Insert Details For Backtracking
         
         If (Iprint .ge. 0) Write(6,1016) 
         Ifinal = 0
         Ep = (Ep+E(3))/Two
         If ((Ep-E(3)) .lt. 0.01d0*E(3)) Then
            Ep = E(3)
            Epsmin = One/Ep
            If (Iprint .ge. 0) Write(6,1021) Epsmin
            Emin = Ep
            Iback = 1
            Iflbvp = 0
         Endif
         
         If (Nss .eq. 1) Then
            Phiaim = Phi(2)+F2*(Ep-E(2))
         Else If (Nss .gt. 1) Then
            If (Iextrap .ne. 1) Then
               Philin = Phi(2)+F2*(Ep-E(2))
               If (Ff2 .gt. 0) Then
                  Phiqua = Phi(1)+(Ep-E(1))*(F1+(Ep-E(2))*Ff2)
               Else
                  Qa = Gg2
                  Qb = G1-(Phi(1)+Phi(2))*Gg2
                  Qc = E(1)-Phi(1)*(G1-Phi(2)*Gg2)-Ep
                  Qd = Qb**2-4.d0*Qa*Qc
                  Phiqua = (-Qb+Sqrt(Qd))/(Two*Qa)
               Endif
               
               If (Ilin .eq. 1) Then
                  Phiaim = Phiqua
                  Phialt = Philin
               Else
                  Phiaim = Philin
                  Phialt = Phiqua
               Endif
            Endif
         Endif
         
C.... Re-Insert Details From Last Problem Successfully Solved.
         
         If (Nss .gt. 0) Then
            Nmsh = Nnewold
            Call Meshdet(Xx, Wrk(Imsh), Nmsh)
            If (.not. Linear) Call Soldet(U, Wrk(Isol), Ncomp, Nmsh, 
     +           Nudim, Ncomp)
         Endif
         
      Endif   
      
C.... *****************************************************************
C.... End Of Logic For An Unsuccessful Continuation Step.
C.... *****************************************************************
      
C.... Pass On Mesh And Solution Details To The Next Continuation 
C.... Problem.
      
      If (Nss .gt. 0) Then
         Givmsh = .true.
         If (.not. Linear) Giveu = .true.
      Endif
      
C.... *****************************************************************
C.... The Program Loops Back To Line 50 For Every Continuation Step.
C.... *****************************************************************
      
      Goto 50
      
C-----------------------------------------------------------------------
 1000 Format(/,' The Number Of (Linear) Differential Equations Is ',I3)
 1001 Format(/,' The Number Of (Nonlinear) Differential Equations Is ',
     +     I3)
 1002 Format(' Left Boundary Point =',D12.4,', Right Boundary Point =',
     +     D12.4)
 1003 Format(' There Are',I5,' Fixed Points In The Mesh - ',10(6d10.4/))
 1004 Format(' Components Of U Requiring Tolerances -',8(7x,I2,1x),
     +     4(/38x,8i10))
 1005 Format(' Corresponding Error Tolerances -',6x,8d10.2,
     +     4(/39x,8d10.2))
 1006 Format (' The Initial Value Of Epsilon Is ',D11.4,/,
     +     ' The Desired Final Value Of Epsilon Is ',D11.4)
 1007 Format(1h ,'The Largest Mesh Size From Workspace, Nmax =',I8)
 1008 Format(1h ,'Nmax = Nucol =', I8)
 1010 Format (/,1x,82('*'),/,' Continuation Step ',I2,
     +        ', Epsilon = ',D9.4,/,1x,82('*'))
 1011 Format (/,1x,'This Is The Final Continuation Step Since ',
     +        'Maxcon =',I4)
 1012 Format (/,1x,82('$'),//,' Tolerances Satisfied For Final Problem',
     +        ' Epsilon = ',D9.4,//,1x,82('$'),/)
 1013 Format(/,' ** Machine Precision (Possibly) Not Sufficient For',
     +        ' Epsilon Less Than ',D9.4) 
 1014 Format (/,1x,82('$'),//,' Final Problem Epsilon = ',D10.4,/,
     +        ' Not Solved, Iflag =',I3,/, ' Try Running Problem Again',
     +        ' With Epsmin Greater Than ',D9.4,/, ' Machine Precision',
     +        ' (Possibly) Not Sufficient For Eps Less Than ',D9.4,//,
     +        1x,82('$'),/) 
 1015 Format (/,1x,82('$'),//,' Final Problem Epsilon = ',D10.4,/,
     +        ' Not Solved, Iflag =',I3,/, ' Try Running Problem Again',
     +        ' With Epsmin Greater Than ',D9.4,//,1x,82('$'),/) 
 1016 Format (/,' ** Failed Step - Bactracking For Larger Value ',
     +          'Of Epsilon')
 1018 Format(' The Final Mesh And Solution Components Are:',//,
     +       5x,'I',10x,'X(I) ',40(14x,A,'(',I2,')'))
 1019 Format(I6,41d19.8)
 1020 Format(/,1x,82('$'),/)
 1021 Format (/,' Continuation Steps Too Small, Change Epsmin To ',D9.4)
 1022 Format (/,' Storage Limit Being Approached, Change Epsmin To ',
     +     D9.4)
C-----------------------------------------------------------------------
      End 



      Subroutine Meshdet(X, X1, Nmsh)
      
      Implicit Double Precision (A-H,O-Z)
      
C.... This Subroutine Is Used For Saving And Re-Inserting Meshes
      
      Dimension X(*), X1(*)
      
      Do 10, J = 1,Nmsh
         X(J) = X1(J)
 10   Continue
      
      Return
      End

      Subroutine Soldet(U, U1, Ncomp, Nmsh, Idim1, Idim2)
      
      Implicit Double Precision (A-H,O-Z)
      
C.... This Subroutine Is Used For Saving And Re-Inserting Solutions.
      
      Dimension U(Idim1,*), U1(Idim2,*)
      
      Do 14, J = 1,Ncomp
         Do 12 K = 1,Nmsh
            U(J,K) = U1(J,K)
 12      Continue
 14   Continue
      
      Return
      End


* File Bnob.f (Nob = 'No Blas')

      Subroutine Blkdcm (Ntop, Nbot, Nrblk, Nbloks, 
     *    Topblk, Array, Botblk, Ipivot, Iflag)

*  This Subroutine Decomposes The Almost Block-Diagonal Matrix A
*  Given By
*
*          ( Topblk                         )
*          (     Array(1)                   )
*          (         Array(2)               )
*          (               .                )
*          (                 .              )
*          (               Array(Nbloks)    )
*          (                         Botblk ),
*
*  Where Topblk Is Ntop By Ntop+Nbot,
*  Array(K) Is Nrblk By Nrblk+Ntop+Nbot, For K = 1,...,Nbloks
*  And Botblk Is Nbot By Ntop+Bnot.
*  The Matrix Thus Contains Nbloks*Nrblk + Ntop + Nbot ( =  Nrow) Rows.

*  The Matrix Is Decomposed Using Gaussian Elimination With
*  Alternate Row And Column Elimination With Partial Pivoting
*  Within Each Block.

*  Ntop Is The Number Of (Conceptual) Rows In Topblk (May Be Zero).

*  Nblk Is The Number Of Rows In Each Array Along The Diagonal.

*  Nbot Is The Number Of (Conceptual) Rows In Botblk (May Be Zero).

*  The Declared Row Dimensions Of Both Topblk And Botblk
*  Must Be Ntop+Nbot.
*  (This Is Slightly Wasteful Of Storage, But Avoids 2 Extra
*  Parameters And Much Complication In Storage And Coding.

*  Topblk Is An Array Of Conceptual Dimension (Ntop,Ntop+Nbot);
*  The Upper Left Block Of A.  Its Declared Row Dimension Must Be
*  Ntop+Nbot.

*  Botblk Is An Array Of Conceptual Dimension (Nbot,Ntop+Nbot);
*  The Lower Right Block Of A.  Its Declared Row Dimension Must Be
*  Ntop+Nbot.

*   Array(*,*,K) Contains The K-Th Matrix Lying Along The Diagonal;
*                Each Block Is Nrblk By Nrblk + Ntop + Nbot.
*   Nrblk Is The Number Of Rows In Each Block
*   The Number Of Columns In Each Block Is Nrblk + Ntop + Nbot.
*   Nbloks Is The Number Of Blocks.

*   Ipivot Is An Integer Array Of Dimension Nrow,
*   Used To Store Information About Interchanges Performed
*   During Elimination.

*   Iflag Is An Output Flag.  Iflag = 0 If Everything Is Ok,
*   = 1 If There Are Invalid Input Parameters, And = -1 If The 
*   Matrix Is Judged To Be Singular Because Of A Too-Small Pivot.

      Implicit Double Precision (A-H, O-Z)

      Dimension Ipivot (*)
      Dimension Topblk(Ntop+Nbot,*), 
     *             Array(Nrblk,Ntop+Nbot+Nrblk,*), 
     *             Botblk(Ntop+Nbot,*)

      Common/Mchprs/ Flmax, Epsmch

      Intrinsic Abs

      Parameter (Zero = 0.0d+0)
      Parameter (Ten = 10.0d+0)


      Iflag = 0
      Novrlp = Ntop + Nbot
      Ncblk = Nrblk + Novrlp
      Pivtol = Ten*Epsmch

*  In Topblk, Apply Ntop Elimination Steps, One For Each Row,
*  Using Column Interchanges.
*  After The I-Th Step Of Reduction, Row I Has Been Reduced To
*  Lower-Triangular Form.


      Do 160 I = 1, Ntop

*  In Row I, Jpvt Gives The Index Of The Element Of Maximum Magnitude 
*  On Or To The Right Of The Diagonal. 

         Colmax = Zero
         Do 100 J = I, Novrlp
            Tabs = Abs(Topblk(I,J))
            If (Tabs .gt. Colmax) Then
               Jpvt = J
               Colmax = Tabs
            Endif
  100    Continue

             
         If (Colmax .le. Pivtol) Then
            Iflag = -1
            Return
         Endif
         Ipivot(I) = Jpvt
         If (Jpvt .ne. I) Then

            Do 110 L = I, Ntop
               Swap = Topblk(L,Jpvt)
               Topblk(L,Jpvt) = Topblk(L,I)
               Topblk(L,I) = Swap
  110       Continue
            Do 120 L = 1, Nrblk
               Swap = Array(L,Jpvt,1)
               Array(L,Jpvt,1) = Array(L,I,1)
               Array(L,I,1) = Swap
  120       Continue

         Endif
         Colpiv = Topblk(I,I)
         Do 150 J = I+1, Novrlp
            Colmlt = Topblk(I,J)/Colpiv
            Topblk(I,J) = Colmlt
            Do 130 L = I+1, Ntop
               Topblk(L,J) = Topblk(L,J) - Colmlt*Topblk(L,I)
  130       Continue
            Do 140 L = 1, Nrblk
               Array(L,J,1) = Array(L,J,1) - Colmlt*Array(L,I,1)
  140       Continue
  150    Continue
            
  160 Continue



*  In Each Diagonal Block Array, K = 1,...,Nbloks

      Do 320 K = 1, Nbloks
         Incr = (K-1)*Nrblk

*  In Each Matrix Array, Apply Nrblk-Ntop Elimination Steps
*  In Columns Ntop+1,...,Nrblk, Using Row Interchanges.

         Do 210 J = Ntop+1, Nrblk
            Jminn = J - Ntop
            Ipvt = Jminn
            Rowmax = Abs(Array(Jminn,J,K))
            Do 170 I = Jminn+1, Nrblk
               Tabs = Abs(Array(I,J,K))
               If (Tabs .gt. Rowmax) Then
                  Ipvt = I
                  Rowmax = Tabs
               Endif
  170       Continue
            If (Rowmax .le. Pivtol) Then
                Iflag = -1
                Return
            Endif

            Ipivot(Incr+J) = Incr+Ipvt+Ntop

*  Interchange Rows If Necessary.

            If (Ipvt .ne. Jminn) Then
               Do 180 L = J, Ncblk
                  Swap = Array(Ipvt,L,K)
                  Array(Ipvt,L,K) = Array(Jminn,L,K)
                  Array(Jminn,L,K) = Swap
  180          Continue
            Endif


*  Compute Multipliers.

            Rowpiv = Array(Jminn,J,K)
            Do 190 I = Jminn+1, Nrblk
               Array(I,J,K) = Array(I,J,K)/Rowpiv
  190       Continue
            Do 200 L = J+1, Ncblk
               Rowmlt = Array(Jminn,L,K)
               Do 200 I = Jminn+1, Nrblk
                  Array(I,L,K) = Array(I,L,K) - Rowmlt*Array(I,J,K)
  200       Continue
  210    Continue


*  Next, Perform Eliminations In Rows Nrblk-Ntop+1 Through Nrblk,
*  With Column Interchanges.

         Do 310 I = Nrblk - Ntop + 1, Nrblk

*           (The Index Ipn Runs From Nrblk+1 To Nrblk+Ntop.)

            Ipn = I + Ntop

*  Find The Index Of The Column With Maximum Element In Row I, Starting
*  With Column Ntop + I.  The Value Of Jpv1 Lies Between 1 And Ncblk-J+1. 
*  The Column Index In Array (*,*,K) Is Jpv1+Ntop. 
            Jpvt = Ipn
            Colmax = Abs(Array(I,Jpvt,K))
            Do 220 J = Ipn+1, Ncblk
               Tabs = Abs(Array(I,J,K))
               If (Tabs .gt. Colmax) Then
                  Jpvt = J
                  Colmax = Tabs
               Endif
  220       Continue

            If (Colmax .le. Pivtol) Then
               Iflag = -1
               Return
            Endif
            Ipivot(Incr+Ipn) = Incr+Jpvt

*  Interchange Columns If Necessary.

            If (Jpvt .ne. Ipn)  Then
               Do 230 L = I, Nrblk
                  Swap = Array(L,Jpvt,K)
                  Array(L,Jpvt,K) = Array(L,Ipn,K)
                  Array(L,Ipn,K) = Swap
  230          Continue

*  Depending On The Value Of K, Interchange Columns In Block K+1 Of 
*  Array Or In Botblk.

               If (K .lt. Nbloks) Then
                  Do 240 L = 1, Nrblk
                     Swap = Array(L,Jpvt-Nrblk,K+1)
                     Array(L,Jpvt-Nrblk,K+1) = Array(L,Ipn-Nrblk,K+1)
                     Array(L,Ipn-Nrblk,K+1) = Swap
  240             Continue
               Else
                  Do 250 L = 1, Nbot
                     Swap = Botblk(L,Jpvt-Nrblk)
                     Botblk(L,Jpvt-Nrblk) = Botblk(L,Ipn-Nrblk)
                     Botblk(L,Ipn-Nrblk) = Swap
  250             Continue
               Endif
            Endif

*  Compute Multipliers And Perform Elimination.

            Colpiv = Array(I,Ipn,K)
            Do 300 J = Ipn+1, Ncblk
               Colmlt = Array(I,J,K)/Colpiv
               Array(I,J,K) = Colmlt
               Do 260 L = I+1, Nrblk
                  Array(L,J,K) = Array(L,J,K) - Colmlt*Array(L,Ipn,K)
  260          Continue

               If (K .lt. Nbloks) Then
                  Do 270 L = 1, Nrblk
                     Array(L,J-Nrblk,K+1) = Array(L,J-Nrblk,K+1)
     *                 - Colmlt*Array(L,Ipn-Nrblk,K+1)
  270             Continue
               Else
                  Do 280 L = 1, Nbot
                     Botblk(L,J-Nrblk) = Botblk(L,J-Nrblk)
     *                 - Colmlt*Botblk(L,Ipn-Nrblk)
  280             Continue
               Endif
  300       Continue

*  310 Ends The Loop Over I
  310    Continue

*     320 Ends The Loop Over K
  320 Continue

*   In Botblk, Apply Nbot Eliminations With Row Interchanges.

      Incr = Nbloks*Nrblk
      Do 390 Jj = 1, Nbot
         J = Jj + Ntop
         Ipvt = Jj
         Rowmax = Abs(Botblk(Ipvt,J))
         Do 330 I = Jj+1, Nbot
            Tabs = Abs(Botblk(I,J))
            If (Tabs .gt. Rowmax) Then
               Ipvt = I
               Rowmax = Tabs
            Endif
  330    Continue
         If (Rowmax .le. Pivtol) Then
            Iflag = -1
            Return
         Endif
         Ipivot(Incr+J) = Incr+Ntop+Ipvt
         If (Ipvt .ne. Jj) Then
            Do 340 L = J, Novrlp
               Swap = Botblk(Ipvt,L)
               Botblk(Ipvt,L) = Botblk(Jj,L)
               Botblk(Jj,L) = Swap
  340       Continue
         Endif

         Rowpiv = Botblk(Jj,J)
         Do 350 I = Jj+1, Nbot
            Botblk(I,J) = Botblk(I,J)/Rowpiv
  350    Continue
         
         Do 380 L = J+1, Novrlp
            Rowmlt = Botblk(Jj,L)
            Do 380 I = Jj+1, Nbot
               Botblk(I,L) = Botblk(I,L) - Rowmlt*Botblk(I,J)
  380    Continue
  390 Continue

      Return
      End 

      Subroutine Blkslv (Ntop, Nbot, Nrblk, Nbloks,
     *       Topblk, Array, Botblk, Ipivot,
     *       B, X)

      Implicit Double Precision (A-H,O-Z)

*  Special Version Of Block-Structured Solution Where We 
*  Assume That Nrblk = Ntop + Nbot.
*  The Right-Hand Side B Is A Vector Of Length (Nbloks+1)*Nrblk.
*  The Solution X Is Represented As A Two-Dimensional Array, 
*  Of Dimension (Nrblk, Nbloks+1).

      Dimension Ipivot(*)
      Dimension Topblk(Ntop+Nbot,*), Array(Nrblk, 2*Nrblk,*),
     *    Botblk(Ntop+Nbot,*)
      Dimension B(*), X(Nrblk,Nbloks+1) 

      Parameter (Zero = 0.0d+0)

      Novrlp = Ntop+Nbot
      Ncblk = Novrlp + Nrblk


      Do 30 J = 1, Ntop
         X(J,1) = B(J)/Topblk(J,J)
         Do 10 I = J+1, Ntop
            B(I) = B(I) - X(J,1)*Topblk(I,J)
   10    Continue
   30 Continue

      Do 120 K = 1, Nbloks
         Incr = (K-1)*Nrblk

*  Forward Modification.

         Do 50 J = 1, Ntop
         Do 50 I = 1, Nrblk
            B(Incr+Ntop+I) = B(Incr+Ntop+I) - X(J,K)*Array(I,J,K)
   50    Continue

*  Forward Elimination.

         Do 80 J = Ntop+1, Nrblk
            Jpivot = Ipivot(Incr+J)
            If (Jpivot .ne. Incr+J) Then
               Swap = B(Incr+J)
               B(Incr+J) = B(Jpivot)
               B(Jpivot) = Swap
            Endif
            Do 70 I = J-Ntop+1, Nrblk
               Ii = I + Incr + Ntop
               B(Ii) = B(Ii) - B(Incr+J)*Array(I,J,K)
   70       Continue
   80    Continue

*  Forward Solution.

         Do 110 J = Nrblk+1, Ntop+Nrblk
            Jmtop = J-Ntop
            X(J,K) = B(Incr+J)/Array(Jmtop,J,K) 
            Do 90 I = J-Ntop+1, Nrblk
               Ii = I + Incr + Ntop
               B(Ii) = B(Ii) - X(J,K)*Array(I,J,K)
   90       Continue
  110    Continue

  120 Continue

*  Forward Modification In Botblk.
   
      Incr = Nbloks*Nrblk

      Do 140 J = 1, Ntop
      Do 140 I = 1, Nbot
         B(Incr+Ntop+I) = B(Incr+Ntop+I) 
     *                      - X(J,Nbloks+1)*Botblk(I,J)
  140 Continue

*  Forward Elimination.

      Do 170 J = Ntop+1, Ntop+Nbot-1
         Jpivot = Ipivot(Incr+J)
         If (Jpivot .ne. Incr+J) Then
            Swap = B(Incr+J)
            B(Incr+J) = B(Jpivot) 
            B(Jpivot) = Swap
         Endif
         Do 160 I = J-Ntop+1, Nbot
            Ii = I + Incr + Ntop
            B(Ii) = B(Ii) - B(Incr+J)*Botblk(I,J)
  160    Continue
  170 Continue

*  Backward Solution In Botblk.

      Do 210 J = Nbot+Ntop, Ntop+1, -1
         X(J,Nbloks+1) = B(Incr+J)/Botblk(J-Ntop,J)
         Do 190 I = 1, J-Ntop-1
            B(Incr+Ntop+I) = B(Incr+Ntop+I)
     *           - X(J,Nbloks+1)*Botblk(I,J)
  190    Continue
  210 Continue

*  Backward Elimination. 

      Do 300 K = Nbloks, 1, -1
         Incr = (K-1)*Nrblk
         Do 240 I = Nrblk, Nrblk-Ntop+1,-1
            Iptp = I+Ntop
            Indx = Iptp-Nrblk
            Dotprd = Zero
            If (Ncblk-Iptp .ge. 1) Then
               Do 220 J = Iptp+1, Ncblk
                  Dotprd = Dotprd + Array(I,J,K)*X(J-Nrblk,K+1)
  220          Continue
            Endif
            X(Indx,K+1) = X(Indx,K+1) - Dotprd 
            Ipvtn = Ipivot(Incr+Iptp)
            If (Incr+Iptp .ne. Ipvtn) Then
*  Convert Ipvtn To An Index Pair.
               Jcol = 1 + (Ipvtn-1)/Nrblk
               Irow = Ipvtn - (Jcol-1)*Nrblk
               Swap = X(Indx,K+1)
               X(Indx,K+1) = X(Irow,Jcol)
               X(Irow,Jcol) = Swap
            Endif
  240    Continue

*   Backward Modification.

         Do 260 J = Nrblk+1, Ncblk
         Do 260 I = 1, Nrblk-Ntop
            Ii = I + Incr + Ntop
            B(Ii) = B(Ii) - X(J,K)*Array(I,J,K)
  260    Continue

*   Backward Solution.

         Do 290 J = Nrblk, Ntop+1, -1
            X(J,K) = B(Incr+J)/Array(J-Ntop,J,K) 
            Do 270 I = 1, J-Ntop-1
               Ii = I + Incr + Ntop
               B(Ii) = B(Ii) - X(J,K)*Array(I,J,K)
  270       Continue
  290    Continue

*  300 Ends The Loop Over K
  300 Continue

*  Backward Elimination In Topblk.


       Do 330 I = Ntop, 1, -1
         Dotprd = X(I,1)
         Do 310 J = I+1, Ntop+Nbot
            Dotprd = Dotprd - Topblk(I,J)*X(J,1)
  310    Continue
         Ipvti = Ipivot(I)
         If (I .ne. Ipvti) Then
*  Convert Ipvti To An Index Pair.
            Jcol = 1 + (Ipvti-1)/Nrblk
            Irow = Ipvti - (Jcol-1)*Nrblk
            Swap = X(I,1)
            X(I,1) = X(Irow,Jcol)
            X(Irow,Jcol) = Swap
         Endif
  330 Continue

      Return
      End

      Subroutine Lufac(N, Ndim, A, Ip, Ier)
      Implicit Double Precision (A-H,O-Z)
      Dimension A(Ndim,N), Ip(N)
      Intrinsic Abs
*  Blas: Daxpy, Dscal, Dswap. Idamax

      Parameter ( Zero = 0.0d+0, One = 1.0d+0 )

*  The Subroutine Lufac Is A Very Simple Code To Compute The 
*  Lu Decomposition (With Partial Pivoting) Of The N By N Matrix A.

*  The Lu Factors Are Overwritten On A.  The Integer Array Ip
*  Reflects The Pairwise Interchanges Performed.  Note That Ip(K)
*  Therefore Does Not Give The Index In The Original Array Of
*  The K-Th Pivot.

*  On Exit, The Error Flag Ier Is Zero When No Zero Pivots Are
*  Encountered.  Otherwise, Ier Is Equal To The Index Of The 
*  Step At Which A Zero Pivot Occurred.

      Ier = 0
      Ip(N) = 0

*  Begin Loop Over Columns 1 Through N-1.  K Is The Current
*  Column Index.

      Do 100 K = 1, N-1

*  Find The Row Index Ipiv Of The Element Of Largest Magnitude In
*  Column K.

         Ipiv = K-1 + Idamax(N-K+1, A(K,K), 1)
         Piv = A(Ipiv,K) 
         If (Piv .eq. Zero) Then
            Ier = K
            Return
         Endif
         Ip(K) = Ipiv

*  Perform Interchanges If Necessary.

         If (Ipiv .ne. K) Then
            Call Dswap(N-K+1, A(Ipiv,K), Ndim, A(K,K), Ndim)
         Endif

*  Save The (Negative) Multipliers In The Subdiagonal Elements Of 
*  Column K.

         Call Dscal(N-K, (-One/Piv), A(K+1,K), 1)

*  Update The Remaining Matrix.  Note That A(I,K) Now Contains
*  The Negative Multipliers.

         Do 50 J = K+1, N
            Call Daxpy(N-K, A(K,J), A(K+1,K), 1, A(K+1,J), 1)
   50    Continue

*  End Of Loop Over Columns.

  100 Continue
      If (A(N,N) .eq. Zero) Ier = N
      Return
      End 



      Subroutine Lusol (N, Ndim, A, Ip, B, X) 
      Implicit Double Precision (A-H, O-Z)
      Dimension A(Ndim,N), Ip(N), B(N), X(N)

*  Blas:  Daxpy, Dcopy

*  The Subroutine Lusol Is A Simple-Minded Routine To Solve A
*  Linear System Whose Lu Factors Have Been Computed By Lufac.
*  On Entry, The Matrix A Should Contain The Lu Factors, And
*  Ip Should Contain The Interchange Array Constructed By Lufac.


*  Copy The Right-Hand Side B Into X.

      Call Dcopy (N, B, 1, X, 1)

*  Forward Solution With L (Unit Lower-Triangular Factor), Which
*  Is Stored In The Strict Lower Triangle Of A.

      Do 20 K = 1, N-1
         Ipiv = Ip(K)
         If (Ipiv .ne. K) Then
            Tem = X(Ipiv)
            X(Ipiv) = X(K)
            X(K) = Tem
         Endif
         Call Daxpy ( N-K, X(K), A(K+1,K), 1, X(K+1), 1 )
   20 Continue

*  Backward Solution With U (Upper-Triangular Factor), Which Is Stored
*  In The Upper Triangle Of A.

      Do 40 Kb = N, 1, -1
         X(Kb) = X(Kb)/A(Kb,Kb)
         Call Daxpy(Kb-1, (-X(Kb)), A(1,Kb), 1, X(1), 1)
   40 Continue

      Return
      End 

      Subroutine Dcopy ( N, X, Incx, Y, Incy )
      Implicit Double Precision (A-H,O-Z)
      Integer            N, Incx, Incy
      Dimension          X( * ), Y( * )

C  Dcopy  Performs The Operation
C
C     Y : =  X
C
C  Nag Fortran 77 Version Of The Blas Routine Dcopy .
C  Nag Fortran 77 O( N ) Basic Linear Algebra Routine.
C
C  -- Written On 26-November-1982.
C     Sven Hammarling, Nag Central Office.

      Integer            I     , Ix    , Iy

      If( N .lt. 1 )Return
      If( ( Incx .eq. Incy ) .and. ( Incy .gt. 0 ) )Then
         Do 10, Iy = 1, 1 + ( N - 1 )*Incy, Incy
            Y( Iy ) = X( Iy )
   10    Continue
      Else
         If( Incx .ge. 0 )Then
            Ix = 1
         Else
            Ix = 1 - ( N - 1 )*Incx
         End If
         If( Incy .gt. 0 )Then
            Do 20, Iy = 1, 1 + ( N - 1 )*Incy, Incy
               Y( Iy ) = X( Ix )
               Ix      = Ix + Incx
   20       Continue
         Else
            Iy = 1 - ( N - 1 )*Incy
            Do 30, I = 1, N
               Y( Iy ) = X( Ix )
               Iy      = Iy + Incy
               Ix      = Ix + Incx
   30       Continue
         End If
      End If

      Return

*     End Of Dcopy .

      End

      Subroutine Daxpy ( N, Alpha, X, Incx, Y, Incy )
      Implicit Double Precision (A-H,O-Z)
      Integer            N, Incx, Incy
      Dimension   X( * ), Y( * )

C  Daxpy  Performs The Operation
C
C     Y : =  Alpha*X + Y
C
C
C  Modified Nag Fortran 77 Version Of The Blas Routine Daxpy .
C
C  -- Written On 3-September-1982.
C     Sven Hammarling, Nag Central Office.

      Integer            I     , Ix    , Iy
      Parameter        ( Zero  = 0.0d+0 )

      If( N     .lt. 1    )Return
      If( Alpha .eq. Zero )Return

      If( ( Incx .eq. Incy ) .and. ( Incx .gt. 0 ) )Then
         If ( Alpha .eq. 1.0d+0 ) Then
            Do 5, Ix = 1, 1 + ( N - 1 )*Incx, Incx
               Y( Ix ) = X( Ix ) + Y( Ix )
 5          Continue
         Else
            Do 10, Ix = 1, 1 + ( N - 1 )*Incx, Incx
               Y( Ix ) = Alpha*X( Ix ) + Y( Ix )
 10         Continue
         Endif
      Else
         If( Incy .ge. 0 )Then
            Iy = 1
         Else
            Iy = 1 - ( N - 1 )*Incy
         End If
         If( Incx .gt. 0 )Then
            Do 20, Ix = 1, 1 + ( N - 1 )*Incx, Incx
               Y( Iy ) = Alpha*X( Ix ) + Y( Iy )
               Iy      = Iy + Incy
   20       Continue
         Else
            Ix = 1 - ( N - 1 )*Incx
            Do 30, I = 1, N
               Y( Iy ) = Alpha*X( Ix ) + Y( Iy )
               Ix      = Ix + Incx
               Iy      = Iy + Incy
   30       Continue
         End If
      End If
      Return

*     End Of Daxpy .

      End

      Double Precision Function Ddot  ( N, X, Incx, Y, Incy )
      Implicit Double Precision (A-H,O-Z)
      Integer           N, Incx, Incy
      Dimension         X( * ), Y( * )

C  Ddot   Returns The Value
C
C     Ddot   = X'Y
C
C
C  Modified Nag Fortran 77 Version Of The Blas Routine Ddot  .
C
C  -- Written On 21-September-1982.
C     Sven Hammarling, Nag Central Office.

      Integer             I     , Ix    , Iy
      Parameter         ( Zero  = 0.0d+0 )

      Sum = Zero
      If( N .ge. 1 )Then
         If( ( Incx .eq. Incy ) .and. ( Incx .gt. 0 ) )Then
            Do 10, Ix = 1, 1 + ( N - 1 )*Incx, Incx
               Sum = Sum + X( Ix )*Y( Ix )
   10       Continue
         Else
            If( Incy .ge. 0 )Then
               Iy = 1
            Else
               Iy = 1 - ( N - 1 )*Incy
            End If
            If( Incx .gt. 0 )Then
               Do 20, Ix = 1, 1 + ( N - 1 )*Incx, Incx
                  Sum = Sum + X( Ix )*Y( Iy )
                  Iy  = Iy  + Incy
   20          Continue
            Else
               Ix = 1 - ( N - 1 )*Incx
               Do 30, I = 1, N
                  Sum = Sum + X( Ix )*Y( Iy )
                  Ix  = Ix  + Incx
                  Iy  = Iy  + Incy
   30          Continue
            End If
         End If
      End If

      Ddot   = Sum
      Return

*     End Of Ddot  .

      End

      Subroutine Dscal ( N, Alpha, X, Incx )
      Implicit Double Precision (A-H,O-Z)
      Integer          N, Incx
      Dimension        X( * )

C  Dscal  Performs The Operation
C
C     X : =  Alpha*X
C
C
C  Modified Nag Fortran 77 Version Of The Blas Routine Dscal .
C
C  -- Written On 26-November-1982.
C     Sven Hammarling, Nag Central Office.

      Integer            Ix
      Parameter        ( One   = 1.0d+0, Zero  = 0.0d+0 )

      If( N .ge. 1 )Then
         If( Alpha .eq. Zero )Then
            Do 10, Ix = 1, 1 + ( N - 1 )*Incx, Incx
               X( Ix ) = Zero
   10       Continue
         Else If( Alpha .eq. ( -One ) )Then
            Do 20, Ix = 1, 1 + ( N - 1 )*Incx, Incx
               X( Ix ) = -X( Ix )
   20       Continue
         Else If( Alpha .ne. One )Then
            Do 30, Ix = 1, 1 + ( N - 1 )*Incx, Incx
               X( Ix ) = Alpha*X( Ix )
   30       Continue
         End If
      End If

      Return

*     End Of Dscal .

      End

      Subroutine Dswap ( N, X, Incx, Y, Incy )
      Implicit Double Precision (A-H,O-Z)
      Integer      N, Incx, Incy
      Dimension    X( * ), Y( * )

C  Dswap  Performs The Operations
C
C     Temp : =  X,   X : =  Y,   Y : =  Temp.
C
C
C  Modified Nag Fortran 77 Version Of The Blas Routine Dswap .
C
C  -- Written On 26-November-1982.
C     Sven Hammarling, Nag Central Office.

      Integer            I     , Ix    , Iy

      If( N .lt. 1 )Return

      If( ( Incx .eq. Incy ) .and. ( Incy .gt. 0 ) )Then
         Do 10, Iy = 1, 1 + ( N - 1 )*Incy, Incy
            Temp    = X( Iy )
            X( Iy ) = Y( Iy )
            Y( Iy ) = Temp
   10    Continue
      Else
         If( Incx .ge. 0 )Then
            Ix = 1
         Else
            Ix = 1 - ( N - 1 )*Incx
         End If
         If( Incy .gt. 0 )Then
            Do 20, Iy = 1, 1 + ( N - 1 )*Incy, Incy
               Temp    = X( Ix )
               X( Ix ) = Y( Iy )
               Y( Iy ) = Temp
               Ix      = Ix + Incx
   20       Continue
         Else
            Iy = 1 - ( N - 1 )*Incy
            Do 30, I = 1, N
               Temp    = X( Ix )
               X( Ix ) = Y( Iy )
               Y( Iy ) = Temp
               Iy      = Iy + Incy
               Ix      = Ix + Incx
   30       Continue
         End If
      End If

      Return

*     End Of Dswap .

      End

      Integer Function Idamax( N, X, Incx )
      Implicit Double Precision (A-H,O-Z)
      Integer         N, Incx
      Dimension       X( * )

C  Idamax Returns The Smallest Value Of I Such That
C
C     Abs( X( I ) ) = Max( Abs( X( J ) ) )
C                      J
C
C  Nag Fortran 77 Version Of The Blas Routine Idamax.
C  Nag Fortran 77 O( N ) Basic Linear Algebra Routine.
C
C  -- Written On 31-May-1983.
C     Sven Hammarling, Nag Central Office.

      Intrinsic           Abs
      Integer             I     , Imax  , Ix

      If( N .lt. 1 )Then
         Idamax = 0
         Return
      End If

      Imax = 1
      If( N .gt. 1 )Then
         Xmax = Abs( X( 1 ) )
         Ix   = 1
         Do 10, I = 2, N
            Ix = Ix + Incx
            If( Xmax .lt. Abs( X( Ix ) ) )Then
               Xmax = Abs( X( Ix ) )
               Imax = I
            End If
   10    Continue
      End If

      Idamax = Imax
      Return

*     End Of Idamax.

      End
      Subroutine Dload ( N, Const, X, Incx )
      Implicit Double Precision (A-H,O-Z)
      Dimension  X(*)
C
C  Dload  Performs The Operation
C
C     X = Const*E,   E' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( N ) Basic Linear Algebra Routine.
C
C  -- Written On 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
      Parameter        ( Zero = 0.0d+0 )

      If( N .lt. 1 )Return

      If( Const .ne. Zero )Then
         Do 10, Ix = 1, 1 + ( N - 1 )*Incx, Incx
            X( Ix ) = Const
   10    Continue
      Else
         Do 20, Ix = 1, 1 + ( N - 1 )*Incx, Incx
            X( Ix ) = Zero
   20    Continue
      End If

      Return

*     End Of Dload .
      End



      Subroutine Maxpy ( Nrow, Ncol, Alpha, Xmat, Nrowy, Ymat )
      Implicit Double Precision (A-H,O-Z)
      Dimension Xmat(Nrow, Ncol), Ymat(Nrowy, Ncol)

*  Subroutine Maxpy Takes As Input The Scalar Alpha And Two Matrices,
*  Xmat And Ymat.  Xmat Has Declared Row Dimension Nrow, And
*  Ymat Has Declared Row Dimension Nrowy, But Both Are
*  Conceptually Nrow By Ncol.  
*  On Output, (New Ymat) Is Alpha*Xmat+ (Old Ymat), By Analogy 
*  With The Vector Blas Routine Saxpy.

      Do 100 J = 1, Ncol
      Do 100 I = 1, Nrow
         Ymat(I,J) = Ymat(I,J) + Alpha*Xmat(I,J)
  100 Continue
      Return
      End


      Subroutine Matcop( Nrow1, Nrow2, Nrow, Ncol, Xmat1, Xmat2 )
      Implicit Double Precision (A-H,O-Z)
      Dimension Xmat1(Nrow1, Ncol), Xmat2(Nrow2, Ncol)

*  Given 2 Matrices Xmat1 And Xmat2, Where Xmat1 Has Declared 
*  Row Dimension Nrow1, Xmat2 Has Declared Row Dimension Nrow2,
*  And Both Have Column Dimension Ncol, The Routine Matcop Copies
*  Rows 1 Through Nrow, And Columns 1 Through Ncol From Xmat1 Into
*  Xmat2.


      If (Nrow .le. 0 .or. Ncol .le. 0) Return
      Do 100 J = 1, Ncol
      Do 100 I = 1, Nrow
           Xmat2(I,J) = Xmat1(I,J)
  100 Continue
      Return
      End

      Subroutine Mtload( Nrow, Ncol, Const, Nrowx, Xmat )
      Implicit Double Precision (A-H,O-Z)
      Dimension Xmat(Nrowx, Ncol)

*  Mtload Sets Elements 1 Through Nrow, 1 Through Ncol, Of The
*  Matrix Xmat (Whose Declared Row Dimension Is Nrowx) To The 
*  Scalar Value Const.  

      If (Nrow .le. 0 .or. Ncol .le. 0)  Return
      Do 100 J = 1, Ncol
      Do 100 I = 1, Nrow
         Xmat(I,J) = Const
  100 Continue
      Return
      End

      Subroutine Mssq  ( Nrow, Ncol, Xmat, Scale, Sumsq )
      Implicit Double Precision (A-H,O-Z)
      Dimension   Xmat(Nrow, *)

*  Given The Nrow By Ncol Matrix Xmat, Mssq Returns Values
*  Scale And Sumsq Such That
*     (Scale**2) * Sumsq = Sum Of Squares Of Xmat(I,J),
*  Where  Scale = Max  Abs(Xmat(I,J)).

*  Mssq Is A Stripped-Down Matrix Version Of The Blas Routine Sssq.

      Intrinsic    Abs
      Parameter   ( One   = 1.0d+0, Zero  = 0.0d+0 )

      Scale = Zero
      Sumsq = One
      If( Nrow .ge. 1 .and. Ncol .ge. 1) Then
         Do 10 I = 1, Nrow
         Do 10 J = 1, Ncol 
            If( Xmat(I,J) .ne. Zero )Then
               Absxij = Abs(Xmat(I,J))
               If( Scale .lt. Absxij ) Then
                  Sumsq = One + Sumsq* (Scale/Absxij)**2
                  Scale = Absxij
               Else
                  Sumsq = Sumsq + (Absxij/Scale)**2
               End If
            End If
   10    Continue
      End If
      Return
      End

      Subroutine Dssq  ( N, X, Incx, Scale, Sumsq )
      Implicit Double Precision (A-H,O-Z)
      Integer            N, Incx
      Dimension   X( * )

*  Given The N-Vector X, Dssq Returns Values Scale And Sumsq Such That
*     (Scale**2) * Sumsq = Sum Of Squares Of X(I),
*  Where  Scale = Max  Abs(X(I)).

      Intrinsic          Abs
      Parameter        ( One   = 1.0d+0, Zero  = 0.0d+0 )

      Scale = Zero
      Sumsq = One
      If( N .ge. 1 )Then
         Do 10, Ix = 1, 1 + ( N - 1 )*Incx, Incx
            If( X( Ix ) .ne. Zero )Then
               Absxi = Abs( X( Ix ) )
               If( Scale .lt. Absxi )Then
                  Sumsq = One   + Sumsq*( Scale/Absxi )**2
                  Scale = Absxi
               Else
                  Sumsq = Sumsq + ( Absxi/Scale )**2
               End If
            End If
   10    Continue
      End If

      Return

*     End Of Dssq.

      End

* File Bvps.f

      Subroutine Bvpsol(Ncomp, Nmsh, Nlbc, Aleft, Aright, 
     *   Nfxpnt, Fixpnt, Ntol, Ltol, Tol, Nmax, Linear, Giveu,  
     *   Givmsh, Xx, Nudim, U, Def, Delu, Rhs, Fval,
     *   Topblk, Botblk, Ajac, Bhold, Chold, Ipvblk, Ipivlu,
     *   Uint, Ftmp, Tmprhs, Dftmp1, Dftmp2, Dgtm, 
     *   Utrial, Rhstri, Xmerit, Xxold, Uold, Tmp, Dsq, Ratdc, 
     *   Etest6, Etest8, Ermx, Ihcomp, Irefin, Def8, Dhold, 
     *   Voldmsh, Fsub, Dfsub, Gsub, Dgsub, Iflbvp, Eps)

      Implicit Double Precision (A-H,O-Z)

      Dimension  Fixpnt(*), Ltol(Ntol), Tol(Ntol)
      Dimension  Xx(*), U(Nudim, *)
      Dimension  Def(Ncomp,*)
      Dimension  Delu(Ncomp, *), Rhs(*), Fval(Ncomp,*)
      Dimension  Topblk(Ncomp,*), Botblk(Ncomp,*)
      Dimension  Ajac(Ncomp, 2*Ncomp, *)
      Dimension  Bhold(Ncomp, Ncomp, *), Chold(Ncomp, Ncomp, *)
      Dimension  Ipivlu(*), Ipvblk(*)
      Dimension  Uint(Ncomp), Ftmp(Ncomp), Dgtm(Ncomp), Tmprhs(*)
      Dimension  Dftmp1(Ncomp, Ncomp), Dftmp2(Ncomp, Ncomp)
      Dimension  Utrial(Ncomp,*), Rhstri(*), Xmerit(Ncomp, *)
      Dimension  Xxold(*), Uold(Ncomp,*), Tmp(Ncomp,12)
      Dimension  Dsq(Ncomp,Ncomp), Ratdc(*)
      Dimension  Etest6(*), Etest8(*), Ermx(*)
      Dimension  Ihcomp(*), Irefin(*)
      Dimension  Def8(Ncomp,*) 
      Dimension  Dhold(3*Ncomp,3*Ncomp), Voldmsh(*)

      Logical Linear, Giveu, Givmsh, Double, Errok

      External Fsub, Dfsub, Gsub, Dgsub

      Common/Mchprs/ Flmax, Epsmch
      Common /Flags/ Ifinal,Iatt,Iback,Iprec
      Common /Convg/ Nits

      Common/Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0

      Intrinsic Max

      Logical Succes, Strctr
      Logical Ludone, Rhsgiv, Maxmsh

      Logical Mchset
      Save Mchset

      Parameter (Zero = 0.0d+0, One = 1.0d+0)
      Parameter (Third = 0.33d+0)
      Parameter (Quan6 = 0.1d+0 )

*  Blas: Dload
*  Double Precision D1mach
      
      Data Mchset/.true./
      Data Fxfct/10.0d+0/
      Data Maxmsh/.false./

*  The Parameter Inumb Is A Counter Used To Limit To Three The Number 
*  Of Mesh Selections Performed For The Final Continuation Problem. 

      Inumb = 0

      If (Mchset) Then
         Flmax = D1mach(2)
         Epsmch = D1mach(3)
         Mchset = .false.
      Endif

*  The Routines Stcon1 And Stcon2 Calculate Integration Constants Stored 
*  In Labeled Common Consts.

      Call Stcon1
      Call Stcon2

*  Set Up Arrays For The Error Tests.


      If (.not. Linear) Then
         Call Dload(Ntol, One, Etest6, 1)
      Else
         Do 10 I = 1, Ntol
            Etest6(I) = One/Max(Quan6, Tol(I)**Third)
   10    Continue
      Endif          

      Nmold = 1
      Strctr = .false.

*
*  If Givmsh Is .true., The Initial Number Of Mesh Points Must Be 
*  Provided By The User In Nmsh, And The Mesh Points Must Be
*  Contained In The Array Xx (Of Dimension Nmsh).
*  Otherwise, Nmsh Is Set To Its Default Value, And A
*  Uniform Initial Mesh Is Created.

      If (.not. Giveu .and. .not. Givmsh) Then
         Nmsh = Nminit
         If (Nmsh .lt. Nfxpnt+2) Nmsh = Nfxpnt + 2
         Call Unimsh(Nmsh, Aleft, Aright, Nfxpnt, Fixpnt, Xx)
      Endif

      If (.not. Giveu) Call Initu(Ncomp, Nmsh, Nudim, U)
      

***** Top Of Logic For 4th Order Solution ****

  400 Continue
      If (Iprint .ge. 0) Write(6,903) Nmsh

*  Set The Def (Deferred Correction) Array To Zero.

      Call Mtload(Ncomp, Nmsh-1, Zero, Ncomp, Def)
      Iorder = 4

*  The Routine Fneval Calls Fsub At The Mesh Xx And The
*  Solution U, And Saves The Values In The Array Fval.

      Call Fneval(Ncomp, Nmsh, Xx, Nudim, U, Fval, Eps)

*  Try To Compute A 4th Order Solution By Solving A System Of Nonlinear
*  Equations.


      If (Linear) Then
         Ludone = .false.

         Call Lineq( Ncomp, Nmsh, Nlbc, 
     *    Ludone, Xx, Nudim, U, Def, 
     *    Delu, Rhs, Fval, 
     *    Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs,
     *    Ajac, Topblk, Botblk, Bhold, Chold, Ipvblk,
     *    Fsub, Dfsub, Gsub, Dgsub, Iflnwt, Eps)

*  Call Fneval To Evaluate The Fval Array At The New Solution U.
*  (Such A Call Is Not Necessary For The Nonlinear Case Because
*  Fval Is Called Within Newteq For The New U.)

         Call Fneval(Ncomp, Nmsh, Xx, Nudim, U, Fval, Eps)

      Else
         Rhsgiv = .false.
         Call Newteq(Ncomp, Nmsh, Nlbc, 
     *        Rhsgiv, Ntol, Ltol, Tol, 
     *        Xx, Nudim, U, Def, 
     *        Delu, Rhs, Fval,
     *        Utrial, Rhstri, 
     *        Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs, Xmerit, 
     *        Ajac, Topblk, Botblk, Bhold, Chold, Ipvblk,
     *        Fsub, Dfsub, Gsub, Dgsub, Itnwt, Iflnwt, Eps)
         If (Iatt .eq. -1) Nits = Max(1,Itnwt)
         If (Iprint .eq. 0) Write(6,999) Itnwt
      Endif

      If (Iflnwt .eq. 0) Then

         Call Dfexcl( Ncomp, Nmsh, Xx, Nudim, U, Def, Linear, Fval, 
     *        Tmp, Fsub, Dfsub, Dsq, Ipivlu, Dhold,
     *        Ntol, Ltol, Tol, Eps )
         
      Else

         Iflbvp = 1
         Return

      Endif

**** Logic For 6th Order ****

      If (Iprint .eq. 1) Write(6,905)

*  Save The 4th Order Solution On This Mesh In Uold. 

      Call Matcop(Nudim, Ncomp, Ncomp, Nmsh, U, Uold)
      Iorder = 6

      If (Linear) Then
        Call Lineq( Ncomp, Nmsh, Nlbc, 
     *    Ludone, Xx, Nudim, U, Def, 
     *    Delu, Rhs, Fval, 
     *    Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs,
     *    Ajac, Topblk, Botblk, Chold, Bhold, Ipvblk,
     *    Fsub, Dfsub, Gsub, Dgsub, Iflnwt, Eps)

      Else
        Call Fixjac(Ncomp, Nmsh, Nlbc, 
     *    Iorder, Ntol, Ltol, Tol, 
     *    Xx, Nudim, U, Def, Def, Delu, 
     *    Rhs, Fval, Utrial, Rhstri, 
     *    Rnsq, Uint, Ftmp, Tmprhs,
     *    Ajac, Topblk, Botblk, Ipvblk,
     *    Fsub, Gsub, Iflnwt, Eps)

*  If The Fixed Jacobian Iterations Fail But Rnsq Is Small,
*  Try A Newton Procedure.  Set Rhsgiv To Indicate That
*  The Right-Hand Side And Fval Have Already Been Evaluated
*  At The Given U.

        If (Iflnwt .eq. -3 .and. Rnsq .lt. Fxfct*Epsmch) Rhsgiv = .true.
        If (Iflnwt .ne. 0) Then 
           Call Newteq(Ncomp, Nmsh, Nlbc,
     *          Rhsgiv, Ntol, Ltol, Tol, 
     *          Xx, Nudim, U, Def, 
     *          Delu, Rhs, Fval,
     *          Utrial, Rhstri, 
     *          Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs, Xmerit, 
     *          Ajac, Topblk, Botblk, Bhold, Chold, Ipvblk,
     *          Fsub, Dfsub, Gsub, Dgsub, Iter, Iflnwt, Eps)
        Endif
      Endif

      If (Iflnwt .ne. 0) Then

         Iflbvp = 1
         Return

      Elseif (Ifinal .eq. 1) Then

         Call Errest (Ncomp, Nmsh, Ntol, Ltol, Tol,
     *        Nudim, U, Uold, Etest6, Errok)
         If (Errok) Then
            Iflbvp = 0
            Return
         Endif

      Endif

***** Logic For Trying To Calculate 8th Order Solution *****

      If (Iprint .eq. 1) Write(6,906)

      Call Matcop(Nudim, Ncomp, Ncomp, Nmsh, U, Uold)

*  For Linear Problems, Calculate The Fval Array For The
*  New Solution U.

      If (Linear) Call Fneval(Ncomp, Nmsh, Xx, Nudim, U, Fval, Eps)

*  Calculate 8th Order Deferred Corrections (The Array Def8).

      Call Df8cal (Ncomp, Nmsh, Xx, Nudim, U, Fval, Def8, Linear,
     *             Tmp, Fsub, Dfsub, Dsq, Ipivlu, Dhold,
     *             Ntol, Ltol, Tol, Eps)

*  For Linear Problems, The Def Array Is The Def8 Array.
*  For Nonlinear Problems, Add The Def8 Array To The 
*  Already-Calculated Def Array.

      If (Linear) Then
         Call Matcop(Ncomp, Ncomp, Ncomp, Nmsh-1, Def8, Def)
      Else
         Call Maxpy(Ncomp, Nmsh-1, One, Def8, Ncomp, Def)
      Endif

      Iorder = 8

      If (Linear) Then
        Call Lineq( Ncomp, Nmsh, Nlbc,  
     *   Ludone, Xx, Nudim, U, Def, 
     *   Delu, Rhs, Fval, 
     *   Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs,
     *   Ajac, Topblk, Botblk, Chold, Bhold, Ipvblk,
     *   Fsub, Dfsub, Gsub, Dgsub, Iflnwt, Eps)
      Else
        Call Fixjac(Ncomp, Nmsh, Nlbc, 
     *    Iorder, Ntol, Ltol, Tol, 
     *    Xx, Nudim, U, Def, Def8, Delu, 
     *    Rhs, Fval, Utrial, Rhstri, 
     *    Rnsq, Uint, Ftmp, Tmprhs,
     *    Ajac, Topblk, Botblk, Ipvblk,
     *    Fsub, Gsub, Iflnwt, Eps)

*  If The Fixed Jacobian Iterations Fail But Rnsq Is Small,
*  Try A Newton Procedure.  Set Rhsgiv To Indicate That
*  The Right-Hand Side And Fval Have Already Been Evaluated
*  At The Given U.

        If (Iflnwt .eq. -3 .and. Rnsq .lt. Fxfct*Epsmch) Rhsgiv = .true.
        If (Iflnwt .ne. 0) Then
           Call Newteq(Ncomp, Nmsh, Nlbc, 
     *          Rhsgiv, Ntol, Ltol, Tol, 
     *          Xx, Nudim, U, Def, 
     *          Delu, Rhs, Fval,
     *          Utrial, Rhstri, 
     *          Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs, Xmerit, 
     *          Ajac, Topblk, Botblk, Bhold, Chold, Ipvblk,
     *          Fsub, Dfsub, Gsub, Dgsub, Iter, Iflnwt, Eps)
        Endif
      Endif
      If (Iflnwt .eq. 0) Then

         Call Conv8( Ncomp, Nmsh, Ntol, Ltol, Tol,
     *              Nfxpnt, Fixpnt, Linear, Nmax,
     *              Xx, Nudim, U, Def8, Uold, 
     *              Ihcomp, Irefin, Ermx,
     *              Etest8, Strctr, Ratdc, Voldmsh,
     *              Double, Nmold, Xxold, Maxmsh, Succes)
         If (Iprec .eq. 2) Then
            If (Iprint .ge. 0) Write(6,1013)
            Iflbvp = 1
            Return
         Endif
      Else 

         Iflbvp = 1
         Return

      Endif
      
      If (Maxmsh) Go To 900

      Iflbvp = 0

      If (Ifinal .eq. 1) Then
        If (Succes) Then
           Return
        Elseif (Iatt .ge. 1) Then
           If (Nmsh .lt. Nmax/2) Then
              Inumb = Inumb+1
           Else
              Iflbvp = 1
              Return
           Endif
        Endif
        If (Iback .ne. 1 .and. Inumb .eq. 3) Then
           Iflbvp = 1
           Return
        Endif
      Elseif (Iatt .eq. 0) Then
        Return
      Endif

      Iatt = Iatt+1

      Goto 400

  900 Continue

* Error Exit---Too Many Mesh Points.

      Iflbvp = 1

      Return

 903  Format(1h ,'The New Mesh Contains',I6,' Points.')
 904  Format(1h ,'Do Not Go On To 6th')
 905  Format(1h ,'Start 6th Order')
 906  Format(1h ,'Start 8th Order')
 999  Format(1h ,'Newton Converged After',I3,' Iterations.')
 1013 Format(/,1x,'** Mesh Cannot Be Defined Within The Bounds Imposed',
     +       ' By The Machine Precision')
      End


      Subroutine Initu(Ncomp, Nmsh, Nudim, U)
      Implicit Double Precision (A-H,O-Z)
      Dimension U(Nudim, *)

      Common/Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0

*  This Routine Must Be Provided To Reset U After Re-Meshing 
*  For Linear Problems Or For Nonlinear Problems
*  When Interpolation Of The Old Solution Is Not Used.
 
*  This Version Sets All Elements Of U To The Constant Uval0.

      If (Iprint .eq. 1) Write(6,99) Uval0
   99 Format(1h ,'Initu, Uval0',1pd15.5)
      Call Mtload(Ncomp, Nmsh, Uval0, Nudim, U)
      Return
      End

      Block Data

*  This Block Data Routine Initializes Nminit (The Initial Number 
*  Of Mesh Points) And Uval0 (The Initial Value For The Trial
*  Solution) To Their Default Values.

      Double Precision Uval0
      Common/Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0

      Data Nminit/7/
      Data Iprint/0/
      Data Itsaim/7/
      Data Maxcon/50/
      Data Uval0/0.0d+0/
      End


      Subroutine Conv8( Ncomp, Nmsh, Ntol, Ltol, Tol,
     *              Nfxpnt, Fixpnt, Linear, Nmax,
     *              Xx, Nudim, U, Def8, Uold, 
     *              Ihcomp, Irefin, Ermx,
     *              Etest8, Strctr, Ratdc, Voldmsh, 
     *              Double, Nmold, Xxold, Maxmsh, Succes)

      Implicit Double Precision (A-H,O-Z)
      Dimension Ltol(Ntol), Tol(Ntol)
      Dimension Fixpnt(*), Etest8(Ntol)
      Dimension Xx(*), U(Nudim,*)
      Dimension Def8(Ncomp,*), Uold(Ncomp,*)
      Dimension Ihcomp(*), Irefin(*)
      Dimension Ermx(*), Xxold(*), Ratdc(*), Voldmsh(*)
      Logical Linear, Strctr, Double, Maxmsh, Succes
***bugfix 2jul01
      Save Nvold

      Common/Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0
      Common /Flags/ Ifinal,Iatt,Iback,Iprec

      Intrinsic Max

      Logical Errok

      Parameter (One = 1.0d+0, Fourth = 0.25d+0, Quan8 = 0.025d+0)

*  Blas: Dload

*  The Newton Iteration Converged For The 8th Order Solution.

      If (Iprint .eq. 1) Write(6,901) 

      Succes = .false.
      Maxmsh = .false.
      If (Ifinal .eq. 1) Then
         If (.not. Linear) Then
            Call Dload(Ntol, One, Etest8, 1)
         Else
            Do 10 I = 1, Ntol
               Etest8(I) = One/Max(Quan8, Tol(I)**Fourth)
 10         Continue
         Endif

*  Check Estimated Error.  For A Nonlinear Problem, All Components
*  Of Etest8 (The Ratios Used In Testing The Error) Are Set To One.  
*  For A Linear Problem, The Components Of Etest8 Are In General
*  Larger Than One.  But If Strctr Is .true. And The Number Of Mesh
*  Points Decreased, We Set The Elements Of Etest8 To One (Which 
*  Makes A Stricter Test).

         If (Linear .and. Strctr .and. Nmsh .lt. Nmold)
     *        Call Dload(Ntol, One, Etest8, 1)
         
         Call Errest (Ncomp, Nmsh, Ntol, Ltol, Tol,
     *        Nudim, U, Uold, Etest8, Errok)
         If (Errok) Then
            Succes = .true.
            Return
         Endif
      Endif

      If (Ifinal .eq. 1 .and. Iatt .ge. 1) Then
         Call Dblmsh (Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh) 
         If (.not. Maxmsh .and. Iprec .ne. 2) Then
***bugfix 2jul01
*           Call Matcop(Nudim, Ncomp, Ncomp, Nmsh, U, Uold)
            Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
            Call Interp(Ncomp, Nmsh, Xx, Nudim, U,
     *           Nmold, Xxold, Uold)
         Endif
         Return
      Endif

*  Perform Selective Mesh Refinement Based On The 8th Order Deferred
*  Corrections.  The Value Of Ipow Indicates That The Error Estimate
*  Is Of Order 6.  Then, For A Nonlinear Problem, Interpolate The
*  Latest Solution Onto The New Mesh.

      Ipow = 6

*  The Array Def8 Will Be Overwritten By Selmsh.

      Call Selmsh(Ncomp, Nmsh, Ntol, Ltol, Tol,
     *     Nfxpnt, Fixpnt, Ipow, Nmax, Nvold,
     *     Xx, Nudim, U, Def8, Irefin, Ihcomp,
     *     Nmold, Xxold, Ermx, Double, Maxmsh, Ratdc, Voldmsh)

      If (.not. Maxmsh .and. Iprec .ne. 2) Then
         If (Linear .and. (Ifinal .eq. 1 .or. Iatt .ne. 0)) Then 
            Call Initu(Ncomp, Nmsh, Nudim, U)
         Else
***bugfix 2jul01
*           Call Matcop(Nudim, Ncomp, Ncomp, Nmsh, U, Uold)
            Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
            Call Interp(Ncomp, Nmsh, Xx, Nudim, U,
     *           Nmold, Xxold, Uold)
         Endif
      Endif

      Return

  901 Format(1h ,'Conv8')
      End

* File Defs.f

      Subroutine Dfexcl (Ncomp, Nmsh, Xx, Nudim, U, Def6, Linear,
     *                 Fval, Tmp, Fsub, Dfsub, Df, Ip, Dhold,
     *                 Ntol, Ltol, Tol, Eps)

      Implicit Double Precision (A-H, O-Z)
      Integer Ncomp, Nmsh
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp, Nmsh)
      Dimension Def6(Ncomp,Nmsh-1), Tmp(Ncomp,8)
      Dimension Df(Ncomp,Ncomp), Ltol(Ntol), Tol(Ntol)
      Dimension Ip(2*Ncomp), Dhold(2*Ncomp,2*Ncomp)
      Dimension St1(200), St2(200), St3(200)
      External Fsub, Dfsub 

      Parameter ( One = 1.0d+0, Two = 2.0d+0 )

      Common /Cons1/ A21,A22,A23,A24,A31,A32,A33,A34,C1,C2,
     +       C16,C26,C123,C223,C14,C24
      Common /Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0
      Common /Flags/ Ifinal,Iatt,Iback,Iprec

      Logical Linear

*  Given The Nmsh Mesh Points Xx, The Estimated Solution
*  U And The Array Fval Of Function Values At (Xx(Im), U(*,Im)),
*  Im = 1,...,Nmsh, Dfexcl Calculates Sixth-Order Implicit 
*  Deferred Correction, Stored In The Array Def6, Indexed
*  Over The Components And Mesh Intervals.
 
*  The Array Tmp Is Workspace For 4 Intermediate Vectors Of
*  Dimension Ncomp.
 
      Do 90 Im = 1, Nmsh-1
         Hmsh = Xx(Im+1) - Xx(Im)

         C16h = C16/Hmsh
         C26h = C26/Hmsh
         Do 10 Ic = 1, Ncomp
            Fvim = Fval(Ic,Im)
            Fvim1 = Fval(Ic,Im+1)
            Uim = U(Ic,Im)
            Uim1 = U(Ic,Im+1)
            Tmp(Ic,3) = C16h*(Uim1-Uim)+C123*Fvim1+C14*Fvim
            Tmp(Ic,4) = C26h*(Uim1-Uim)+C223*Fvim1+C24*Fvim
            St1(Ic) = (Uim+Uim1)/Two
            St2(Ic) = A21*Fvim + A24*Fvim1
            St3(Ic) = A31*Fvim + A34*Fvim1
 10      Continue

         Xxc1 = Xx(Im)+C1*Hmsh
         Xxc2 = Xx(Im)+C2*Hmsh

         Do 60 Nit = 1, 10

         Do 20 Ic = 1,Ncomp
            Tmp3 = Tmp(Ic,3)
            Tmp4 = Tmp(Ic,4)
            Tmp(Ic,1)  = St1(Ic) + Hmsh*(St2(Ic) + A22*Tmp3 + A23*Tmp4)
            Tmp(Ic,2)  = St1(Ic) + Hmsh*(St3(Ic) + A32*Tmp3 + A33*Tmp4)
 20      Continue

         Call Fsub (Ncomp,Xxc1,Tmp(1,1),Tmp(1,5), Eps)
         Call Fsub (Ncomp,Xxc2,Tmp(1,2),Tmp(1,6), Eps)

         Call Dfsub (Ncomp,Xxc1,Tmp(1,1),Df, Eps)
         Do 30 I = 1, Ncomp
            Tmp(I,5) = Tmp(I,5)-Tmp(I,3)
            Tmp(I,6) = Tmp(I,6)-Tmp(I,4)
            Do 25 J = 1,Ncomp
               Dfij = Hmsh*Df(I,J)
               Dhold(I,J) = -A22*Dfij
               Dhold(I,J+Ncomp) = -A23*Dfij
 25         Continue
 30      Continue

         Call Dfsub(Ncomp,Xxc2,Tmp(1,2),Df, Eps)
         Do 35 I = 1, Ncomp
            Do 32 J = 1, Ncomp
               Dfij = Hmsh*Df(I,J)
               Dhold(I+Ncomp,J) = -A32*Dfij
               Dhold(I+Ncomp,J+Ncomp) = -A33*Dfij
 32         Continue
 35      Continue

         Do 40 I = 1,Ncomp
           Dhold(I,I) = Dhold(I,I) + One
           Dhold(I+Ncomp,I+Ncomp) = Dhold(I+Ncomp,I+Ncomp) + One
 40      Continue   

         Call Lufac(2*Ncomp,2*Ncomp,Dhold,Ip,Ier)
         Call Lusol(2*Ncomp,2*Ncomp,Dhold,Ip,Tmp(1,5),Tmp(1,7))

         Do 45 I = 1,Ncomp
            Tmp(I,3) = Tmp(I,3) + Tmp(I,7)
            Tmp(I,4) = Tmp(I,4) + Tmp(I,8)
 45      Continue

         If (Linear) Goto 70

         Jc = 0
         Do 50 I = 1, Ntol 
           Ii = Ltol(I)
           Er = Tol(I)
           If (Abs(Tmp(Ii,7)) .gt. Er*Max(One,Abs(Tmp(Ii,3)))  .or. 
     *         Abs(Tmp(Ii,8)) .gt. Er*Max(One,Abs(Tmp(Ii,4)))) Jc = 1
 50      Continue

         If (Jc .eq. 0) Goto 70

 60      Continue

         If (Iprec .eq. 0) Then
            If (Iprint .ge. 0) Write(6,900) Eps
            Iprec = 1
         Endif

 70      Continue

         Do 80 Ic = 1, Ncomp
            Def6(Ic,Im) = (Hmsh/12.d+0)*(Fval(Ic,Im)+
     *                      5.d+0*(Tmp(Ic,3)+Tmp(Ic,4))+Fval(Ic,Im+1))-
     *                      U(Ic,Im+1)+U(Ic,Im)
 80      Continue
 
 90   Continue

      Return

 900  Format(/,'** Warning - Possibly Approaching Machine Precision ',
     *         'Beyond Epsilon  = ',D10.3,/) 

      End 


      Subroutine Df8cal (Ncomp, Nmsh, Xx, Nudim, U, Fval, Def8, Linear,
     *                   Tmp, Fsub, Dfsub, Df, Ip, Dhold, Ntol,
     *                   Ltol, Tol, Eps) 

*   Given The Mesh Points Xx, The Solution U, And The Function 
*   Values Fval, Df8cal Computes Eighth-Order Deferred Corrections, 
*   Which Are Stored In Def8.
*   The Array Tmp Is Workspace For 9 Intermediate Vectors.

      Implicit Double Precision (A-H, O-Z)
      Integer Ncomp, Nmsh
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp,Nmsh)
      Dimension Def8(Ncomp, Nmsh-1), Tmp(Ncomp,12)
      Dimension Df(Ncomp,Ncomp), Ltol(Ntol), Tol(Ntol)
      Dimension Ip(3*Ncomp), Dhold(3*Ncomp,3*Ncomp)
      Dimension St1(200), St2(200), St3(200), St4(200)
      External Fsub, Dfsub 

      Parameter ( One = 1.0d+0, Two = 2.0d+0 )

      Common/Cons2/ A21,A22,A23,A24,A25,A31,A32,A34,A35,A41,A42,
     +      A43,A44,A45,B1,B2,B3,C1,C2,C3,C16,C26,C36,C123,C223,
     +      C323,C14,C24,C34
      Common /Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0
      Common /Flags/ Ifinal,Iatt,Iback,Iprec

      Logical Linear

      Do 110 Im = 1, Nmsh-1
         Hmsh = Xx(Im+1) - Xx(Im)

         C16h = C16/Hmsh
         C26h = C26/Hmsh
         C36h = C36/Hmsh

         Do 10 Ic = 1, Ncomp
             Fvim = Fval(Ic,Im)
             Fvim1 = Fval(Ic,Im+1)
             Uim = U(Ic,Im)
             Uim1 = U(Ic,Im+1)
             Tmp(Ic,4) = C16h*(Uim1-Uim)+C123*Fvim1+C14*Fvim
             Tmp(Ic,5) = C26h*(Uim1-Uim)+C223*Fvim1+C24*Fvim
             Tmp(Ic,6) = C36h*(Uim1-Uim)+C323*Fvim1+C34*Fvim
             St1(Ic) = (Uim+Uim1)/Two
             St2(Ic) = A21*Fvim + A25*Fvim1
             St3(Ic) = A31*Fvim + A35*Fvim1
             St4(Ic) = A41*Fvim + A45*Fvim1
 10      Continue

         Xxc1 = Xx(Im)+C1*Hmsh
         Xxc2 = Xx(Im)+C2*Hmsh
         Xxc3 = Xx(Im)+C3*Hmsh

         Do 80 Nit = 1, 10

            Do 20 Ic = 1, Ncomp   
              Tmp4 = Tmp(Ic,4)
              Tmp5 = Tmp(Ic,5)
              Tmp6 = Tmp(Ic,6)
              Tmp(Ic,1) = St1(Ic) + Hmsh*(St2(Ic) + A22*Tmp4 + A23*Tmp5 
     +             + A24*Tmp6)
              Tmp(Ic,2) = St1(Ic) + Hmsh*(St3(Ic) + A32*Tmp4 + A34*Tmp6)
              Tmp(Ic,3) = St1(Ic) + Hmsh*(St4(Ic) + A42*Tmp4 + A43*Tmp5 
     +             + A44*Tmp6)
 20         Continue

            Call Fsub (Ncomp,Xxc1,Tmp(1,1),Tmp(1,7), Eps)
            Call Fsub (Ncomp,Xxc2,Tmp(1,2),Tmp(1,8), Eps)
            Call Fsub (Ncomp,Xxc3,Tmp(1,3),Tmp(1,9), Eps)

            Call Dfsub(Ncomp,Xxc1,Tmp(1,1),Df, Eps)
            Do 30 I = 1, Ncomp
               Tmp(I,7) = Tmp(I,7)-Tmp(I,4)
               Tmp(I,8) = Tmp(I,8)-Tmp(I,5)
               Tmp(I,9) = Tmp(I,9)-Tmp(I,6)
               Do 25 J = 1,Ncomp
                  Dfij = Hmsh*Df(I,J)
                  Dhold(I,J) = -A22*Dfij
                  Dhold(I,J+Ncomp) = -A23*Dfij
                  Dhold(I,J+2*Ncomp) = -A24*Dfij
 25            Continue
 30         Continue

            Call Dfsub(Ncomp,Xxc2,Tmp(1,2),Df, Eps)          
            Do 40 I = 1, Ncomp
                Do 35 J = 1, Ncomp
                   Dfij = Hmsh*Df(I,J)
                   Dhold(I+Ncomp,J) = -A32*Dfij
                   Dhold(I+Ncomp,J+Ncomp) = 0.d+0
                   Dhold(I+Ncomp,J+2*Ncomp) = -A34*Dfij
 35             Continue
 40         Continue     

            Call Dfsub(Ncomp,Xxc3,Tmp(1,3),Df, Eps)
            Do 50 I = 1, Ncomp
                 Do 45 J =  1, Ncomp
                    Dfij = Hmsh*Df(I,J)
                    Dhold(I+2*Ncomp,J) = -A42*Dfij
                    Dhold(I+2*Ncomp,J+Ncomp) = -A43*Dfij
                    Dhold(I+2*Ncomp,J+2*Ncomp) = -A44*Dfij
 45              Continue
 50         Continue     

            Do 60 I = 1, Ncomp
               Dhold(I,I) = Dhold(I,I) + One
               Dhold(I+Ncomp,I+Ncomp) = Dhold(I+Ncomp,I+Ncomp) + One
               Dhold(I+2*Ncomp,I+2*Ncomp) = 
     *              Dhold(I+2*Ncomp,I+2*Ncomp) + One
 60         Continue

            Call Lufac(3*Ncomp,3*Ncomp,Dhold,Ip,Ier) 
            Call Lusol(3*Ncomp,3*Ncomp,Dhold,Ip,Tmp(1,7),Tmp(1,10))

            Do 65 I = 1,Ncomp
               Tmp(I,4) = Tmp(I,4) + Tmp(I,10)
               Tmp(I,5) = Tmp(I,5) + Tmp(I,11)
               Tmp(I,6) = Tmp(I,6) + Tmp(I,12)
 65         Continue

            If (Linear) Goto 90

            Jc = 0
            Do 70 I = 1, Ntol 
               Ii = Ltol(I)
               Er = Tol(I)
               If (Abs(Tmp(Ii,10)) .gt. Er*Max(One,Abs(Tmp(Ii,4))) .or. 
     *           Abs(Tmp(Ii,11)) .gt. Er*Max(One,Abs(Tmp(Ii,5))) .or. 
     *           Abs(Tmp(Ii,12)) .gt. Er*Max(One,Abs(Tmp(Ii,6)))) Jc = 1
 70         Continue

            If (Jc .eq. 0) Goto 90

 80      Continue

         If (Iprec .eq. 0) Then
            If (Iprint .ge. 0) Write(6,900) Eps
            Iprec = 1
         Endif

 90      Continue

         Do 100 Ic = 1, Ncomp
           Def8(Ic,Im) = Hmsh*(B1*(Fval(Ic,Im)+Fval(Ic,Im+1))+
     *                   B2*(Tmp(Ic,4)+Tmp(Ic,6))+B3*Tmp(Ic,5))-
     *                   U(Ic,Im+1)+U(Ic,Im)
 100     Continue

 110  Continue  

      Return

 900  Format(/,'** Warning - Possibly Approaching Machine Precision ',
     *         'Beyond Epsilon  = ',D10.3,/) 

      End 


      Subroutine Stcon1
      Implicit Double Precision(A-H,O-Z)
      Common/Cons1/A21,A22,A23,A24,A31,A32,A33,A34,C1,C2,
     +      C16,C26,C123,C223,C14,C24

      Parameter ( One = 1.0d+0, Two = 2.0d+0,  Three = 3.0d+0 )
      Parameter ( Four = 4.0d+0, Five = 5.0d+0,  Six = 6.0d+0 )

      Rt5 = Sqrt(5.0d+0)

      A21 = (Six + Rt5)/120.0d0
      A22 = -Rt5/120.0d0
      A23 = (-13.d0*Rt5)/120.0d0
      A24 = (-Six + Rt5)/120.d0

      A31 = (Six-Rt5)/120.0d0
      A32 = (13.0d0*Rt5)/120.0d0
      A33 = Rt5 / 120.0d0
      A34 = (-Six - Rt5)/120.d0

      C1 = (Five - Rt5)/10.0d0
      C2 = (Five + Rt5)/10.0d0

      C12 = C1*C1
      C22 = C2*C2

      C16 = Six*(C1 - C12)
      C26 = Six*(C2 - C22)

      C123 = Three*C12 - Two*C1
      C223 = Three*C22 - Two*C2

      C14 = One - Four*C1 + Three*C12
      C24 = One - Four*C2 + Three*C22

      Return
      End 

      Subroutine Stcon2
      Implicit Double Precision(A-H,O-Z)
      Common/Cons2/ A21,A22,A23,A24,A25,A31,A32,A34,A35,A41,A42,
     +      A43,A44,A45,B1,B2,B3,C1,C2,C3,C16,C26,C36,C123,C223,
     +      C323,C14,C24,C34

      Parameter ( One = 1.0d+0, Two = 2.0d+0, Three = 3.0d+0 )
      Parameter ( Four = 4.0d+0, Six = 6.0d+0 )

      Rt21 = Sqrt(21.0d+0)

      A21 = One/28.d0 + Three*Rt21/1960.d0
      A22 = -Rt21/280.d0
      A23 = -32.d0*Rt21/735.d0
      A24 = -23.d0*Rt21/840.d0
      A25 = -One/28.d0 + Three*Rt21/1960.d0

      A31 = One/64.d0
      A32 = 7.d0*Rt21/192.d0
      A34 = -7.d0*Rt21/192.d0
      A35 = -One/64.d0

      A41 = One/28.d0 - Three*Rt21/1960.d0
      A42 = 23.d0*Rt21/840.d0
      A43 = 32.d0*Rt21/735.d0
      A44 = Rt21/280.d0
      A45 = -(One/28.d0) - Three*Rt21/1960.d0

      B1 = One/20.0d0
      B2 = 49.0d0/180.d0
      B3 = 16.0d0/45.d0

      C1 = One/Two - Rt21/14.d0
      C2 = One/Two
      C3 = One/Two + Rt21/14.d0

      C12 = C1*C1
      C22 = C2*C2
      C32 = C3*C3

      C16 = Six*(C1 - C12)
      C26 = Six*(C2- C22)
      C36 = Six*(C3 - C32)

      C123 = Three*C12 - Two*C1
      C223 = Three*C22 - Two*C2
      C323 = Three*C32 - Two*C3

      C14 = One - Four*C1 + Three*C12
      C24 = One - Four*C2 + Three*C22
      C34 = One - Four*C3 + Three*C32

      Return
      End 

* File Eqs.f

      Subroutine Fixjac(Ncomp, Nmsh, Nlbc,  
     *    Iorder, Ntol, Ltol, Tol, 
     *    Xx, Nudim, U, Defcor, Defnew, Delu, 
     *    Rhs, Fval, Utrial, Rhstri, 
     *    Rnsq, Uint, Ftmp, Tmprhs,
     *    Ajac, Topblk, Botblk, Ipivot,
     *    Fsub, Gsub, Iflag, Eps)

* Fixed Jacobian Iterations.

      Implicit Double Precision (A-H,O-Z)

      Dimension  Ltol(Ntol), Tol(Ntol)
      Dimension  Xx(Nmsh), U(Nudim,Nmsh), Defcor(Ncomp,Nmsh-1) 
      Dimension  Defnew(Ncomp,Nmsh-1), Delu(Ncomp,Nmsh)
      Dimension  Rhs(Ncomp*Nmsh), Fval(Ncomp,Nmsh)
      Dimension  Utrial(Ncomp,Nmsh), Rhstri(Ncomp*Nmsh)
      Dimension  Uint(Ncomp), Ftmp(Ncomp), Tmprhs(Ncomp*Nmsh)
      Dimension  Ajac(Ncomp, 2*Ncomp, Nmsh-1)
      Dimension  Topblk(Ncomp,Ncomp), Botblk(Ncomp,Ncomp)
      Dimension  Ipivot(Ncomp*Nmsh)
      Logical    Better

      External   Fsub, Gsub

      Common/Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0
      Common/Mchprs/ Flmax, Epsmch

      Intrinsic  Abs, Max 

*  Blas: Dcopy, Dssq

      Parameter ( One    = 1.0d+00 )
      Parameter ( Xlarge = 1.0d+6, Huge = 1.0d+30, Lmtfrz = 8)
      Parameter ( Rngrow = 16.0d+0 )
      Parameter ( Tolfct = 0.1d+0 )

*  The Iteration Scheme Uses A Fixed Jacobian Matrix To Solve For 
*  Correction Vectors, Once There Has Been Convergence Of The Newton
*  Iterations On This Mesh.   It Is Assumed That The Lu
*  Factors Of The Jacobian Have Been Left Unaltered Since
*  Their Calculation.

      If (Iprint .eq. 1) Write(6,901)
      Ninter = Nmsh - 1
      Rnold = Flmax

*  Evaluate The Right-Hand Side Rhstri At The Initial Solution U By
*  Adding The New Deferred Corrections To The Already-Calculated
*  Rhs Vector.

      Call Dcopy(Nlbc, Rhs, 1, Rhstri, 1)
      Ind = Nlbc
      Do 10 Im = 1, Ninter
      Do 10 Ic = 1, Ncomp
         Ind = Ind + 1
         Rhstri(Ind) = Rhs(Ind) + Defnew(Ic, Im)
   10 Continue
      Ind = Ninter*Nmsh + Nlbc + 1
      Call Dcopy(Ncomp-Nlbc, Rhs, 1, Rhstri, 1)

      Call Dssq  ( Nmsh*Ncomp, Rhstri, 1, Scale, Sumsq )
      Rnsq = (Scale**2)*Sumsq

      Iter = 0
      
*  If The Initial Right-Hand Side Is Too Large, Do Not Even Attempt To 
*  Solve The Nonlinear Equations.

      If (Rnsq .gt. Huge .or. 
     *      (Iorder .eq.  8 .and. Rnsq .gt. Xlarge)) Then
         If (Iprint .eq. 1) Write(6,902) Rnsq
         Iflag = -2
         Return
      End If
      Call Dcopy(Ncomp*Nmsh, Rhstri, 1, Rhs, 1)

*  Statement 100 Is The Top Of The Iteration Loop.

  100 Continue

*  If Rnsq Is Sufficiently Small, Terminate Immediately.

      If (Iprint .eq. 1) Write(6,903) Iter, Rnsq
      If (Rnsq .le. Epsmch) Then
         Iflag = 0
         Return
      Endif

      Iter = Iter + 1

*  Solve For The Step Delu By Solving A System Involving The Fixed 
*  Jacobian (Whose Lu Factors Are Saved).  Copy The Rhs Array Into 
*  Tmprhs, Which Is Then Overwritten By Blkslv.

      Call Dcopy(Ncomp*Nmsh, Rhs, 1, Tmprhs, 1) 

      Call Blkslv (Nlbc, Ncomp-Nlbc, Ncomp, Ninter, 
     *       Topblk, Ajac, Botblk, Ipivot, Tmprhs, Delu)


*  Compute The Trial Point Utrial By Adding Delu To U.

      Call Matcop( Nudim, Ncomp, Ncomp, Nmsh, U, Utrial )
      Call Maxpy( Ncomp, Nmsh, One, Delu, Ncomp, Utrial )

*  Compute The Right-Hand Side Vector Rhstri And Its Squared
*  Two-Norm At The Trial Point.

      Rnold = Rnsq
      Call Fneval(Ncomp, Nmsh, Xx, Ncomp, Utrial, Fval, Eps)
      Call Rhscal (Ncomp, Nmsh, Nlbc, Xx, Ncomp, Utrial, Defcor,
     *   Fsub, Gsub, Rhstri, Rnsq, Fval, Ftmp, Uint, Eps) 

*  If Rnsq Strictly Decreased, Update The Solution Vector U 
*  And The Right-Hand Side Rhs.

      Better = .false.
      If (Rnsq .lt. Rnold) Then
         Better = .true.
         Call Matcop( Ncomp, Nudim, Ncomp, Nmsh, Utrial, U )
         Call Dcopy( Ncomp*Nmsh, Rhstri, 1, Rhs, 1 )
      Endif

*  Stop The Fixed Jacobian Iterations If There Have Been Too
*  Many Iterations, Or If Rnsq Has Not Decreased By A Factor
*  Of At Least Rngrow.

      If (Iter .ge. Lmtfrz .or. Rnsq .gt. (Rnold/Rngrow)) Then
         If (Better) Then

*  Setting Iflag To -3 Signals That, Although The Fixed Jacobian
*  Iterations Did Not Succeed, The Current Point Was An Improvement
*  On The Previous One.  Hence, If We Switch To A Newton Procedure,
*  The Right-Hand Side Does Not Need To Be Recalculated.

            Iflag = -3
         Else
            Iflag = -2
         Endif
         If (Iprint .eq. 1) Write(6,904) Iflag
         Return
      Endif

*  Test For Convergence Using The Ratio Abs((Change In U)/Max(U,1)).

      Do 150 Im = 1, Nmsh
      Do 150 It = 1, Ntol
         Itol = Ltol(It)
         Er = Abs(Delu(Itol,Im))/Max(Abs(U(Itol,Im)), One)
         If (Er .gt. Tolfct*Tol(It)) Go To 100
  150 Continue

*  To Exit From The Loop Here, The Convergence Tests Have
*  Been Passed.

      If (Iprint .eq. 1) Write(6,905) Iter, Rnsq

      Iflag = 0
      Return
  901 Format(1h ,'Fixed Jacobian Iterations')
  902 Format(1h ,'Large Residual, Rnsq  = ',1pe12.4)
  903 Format(1h ,'Iter, Rnsq',I5,1pe11.3)
  904 Format(1h ,'Failure Of Fixed Jacobian, Iflag  = ',I5)
  905 Format(1h ,'Fixed Jacobian Convergence',I5,1pe11.3)
      End 



      Subroutine Lineq( Ncomp, Nmsh, Nlbc, 
     *    Ludone, Xx, Nudim, U, Defcor, 
     *    Delu, Rhs, Fval, 
     *    Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs,
     *    Ajac, Topblk, Botblk, Bhold, Chold, Ipivot,
     *    Fsub, Dfsub, Gsub, Dgsub, Iflag, Eps)

      Implicit Double Precision (A-H,O-Z)

      Dimension  Xx(Nmsh), U(Nudim, Nmsh), Defcor(Ncomp, Nmsh-1)
      Dimension  Delu(Ncomp, Nmsh), Rhs(Ncomp*Nmsh)
      Dimension  Fval(Ncomp,Nmsh), Uint(Ncomp), Ftmp(Ncomp)
      Dimension  Topblk(Ncomp,*), Botblk(Ncomp,*)
      Dimension  Dftmp1(Ncomp, Ncomp), Dftmp2(Ncomp, Ncomp)
      Dimension  Dgtm(Ncomp), Tmprhs(Ncomp*Nmsh)
      Dimension  Ajac(Ncomp, 2*Ncomp, Nmsh-1),
     *               Bhold(Ncomp, Ncomp, Nmsh-1),
     *               Chold(Ncomp, Ncomp, Nmsh-1)
      Dimension  Ipivot(Ncomp*Nmsh)

      Logical    Ludone
      External   Fsub, Dfsub, Gsub, Dgsub

      Common/Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0

*  Blas: Dcopy, Dload

      Parameter  ( One = 1.0d+0, Zero = 0.0d+0 )

*  The Routine Lineq Calculates The Newton Step For A Linear
*  Problem.  The Newton Step Is Exact Unless The Jacobian
*  Matrix Is Singular.

      Ninter = Nmsh - 1

      If (.not. Ludone) Then

*  Compute The Right-Hand Side Vector Rhs.
 
         Call Lnrhs (Ncomp, Nmsh, Nlbc, Xx, Nudim, U, 
     *          Fsub, Gsub, Rhs, Rnsq, Fval, Ftmp, Uint, Eps) 


*  If The Jacobian For This Mesh Has Not Previously Been
*  Calculated And Factorized Successfully, Call Jaccal.
*  The Block-Structured Jacobian Matrix Is Stored In Three 
*  Matrices (Topblk, Ajac, And Botblk).
*  The Matrices Bhold And Chold Are Also Calculated In Jaccal,
*  And Are Saved For Later Use In Outer Routines.
         
         Call Jaccal (Ncomp, Nmsh, Nlbc, 
     *      Xx, Nudim, U, Fval, Dgtm, Dftmp1, Dftmp2, Uint,
     *      Ajac, Topblk, Botblk, Bhold, Chold,
     *      Dfsub, Dgsub, Eps)

*  Call Blkdcm To Calculate The Lu Factors Of The Jacobian.
*  The Factors Are Overwritten On The Matrices Topblk, Ajac And Botblk.
*  Interchanges Are Represented In The Integer Array Ipivot. 

         Call Blkdcm (Nlbc, Ncomp-Nlbc, Ncomp, Ninter, 
     *      Topblk, Ajac, Botblk, Ipivot, Iflfac)
         
*  If The Jacobian Is Considered Singular, Iflfac Will Be -1 On 
*  Exit From Blkdcm.
         
         If (Iflfac .eq. (-1)) Then
            If (Iprint .eq. 1) Write(6,901)
            Iflag = -1 
            Ludone = .false.
            Return
         Endif
         Ludone = .true.

*  Copy The Rhs Into The Temporary Vector Tmprhs, Which Will Be
*  Overwritten By Blkslv.

         Call Dcopy(Ncomp*Nmsh, Rhs, 1, Tmprhs, 1)

      Else

*  The Right-Hand Side Is The Deferred Correction Array,
*  Padded With Zeros At The Boundary Conditions.

         Call Dload(Nlbc, Zero, Tmprhs(1), 1)
         Do 100 Im = 1, Ninter
            Loc = (Im-1)*Ncomp + Nlbc + 1
            Call Dcopy(Ncomp, Defcor(1,Im), 1, Tmprhs(Loc), 1)
  100    Continue
         Nrhs = Ninter*Ncomp + Nlbc + 1
         Call Dload(Ncomp-Nlbc, Zero, Tmprhs(Nrhs), 1)

      Endif

*  Solve For The Newton Step Delu. 

      Call Blkslv (Nlbc, Ncomp-Nlbc, Ncomp, Ninter,
     *       Topblk, Ajac, Botblk, Ipivot, Tmprhs, Delu)


*  Since The Problem Is Linear, The Newton Step  Is Exact.  The 
*  New U Array Is Obtained By Adding Delu To U.

      Call Maxpy ( Ncomp, Nmsh, One, Delu, Nudim, U )

      Iflag = 0
      Return

  901 Format(1h ,'Singular Matrix')
      End


      Subroutine Newteq(Ncomp, Nmsh, Nlbc,  
     *    Rhsgiv, Ntol, Ltol, Tol, 
     *    Xx, Nudim, U, Defcor, 
     *    Delu, Rhs, Fval,
     *    Utrial, Rhstri, 
     *    Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs, Xmerit, 
     *    Ajac, Topblk, Botblk, Bhold, Chold, Ipivot,
     *    Fsub, Dfsub, Gsub, Dgsub, Iter, Iflag, Eps)

      Implicit Double Precision (A-H,O-Z)

      Dimension  Ltol(*), Tol(*), Xx(*)
      Dimension  Fval(Ncomp,*)
      Dimension  U(Nudim, *), Delu(Ncomp, *), Utrial(Ncomp,Nmsh)
      Dimension  Rhs(Ncomp*Nmsh),  Defcor(Ncomp,*)
      Dimension  Ajac(Ncomp, 2*Ncomp, *)
      Dimension  Topblk(Ncomp,*), Botblk(Ncomp,*)
      Dimension  Ftmp(*), Uint(*), Dgtm(Ncomp)
      Dimension  Dftmp1(Ncomp, Ncomp), Dftmp2(Ncomp, Ncomp)
      Dimension  Bhold(Ncomp, Ncomp, *), Chold(Ncomp, Ncomp, *)
      Dimension  Ipivot(*)
      Dimension  Rhstri(Ncomp*Nmsh)
      Dimension  Tmprhs(Ncomp*Nmsh), Xmerit(Ncomp, Nmsh)

      Logical Rhsgiv

      External   Fsub, Dfsub, Gsub, Dgsub

      Parameter  ( Zero   = 0.0d+0, One    = 1.0d+0 )
      Parameter ( Two = 2.0d+0 )

      Common/Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0
      Common/Mchprs/ Flmax, Epsmch

      Logical Gtpdeb, Imprvd, Braktd, Crampd, Extrap, Vset, Wset
      Save  Gtpdeb, Mfsrch, Epsaf, Epsag, Rmu, Tolabs, Alfmax
      Save  Tolrel, Toltny
      
      Intrinsic  Abs, Max

*  Blas: Dcopy
  
      Logical Frscal
      Save Frscal

      Parameter (Cnvfct = 0.1d+0 )

      Data  Alfsml/1.0d-4/,  Alfmax/1.1d+0/
      Data  Imerit/1/, Lmtnwt/39/
      Data  Shrfct/100.0d+0/, Stpfct/2.0d+0/
      Data  Gtpdeb/.false./, Mfsrch/5/
      Data  Rmu/1.0d-6/
      Data  Frscal/.true./

*  The Routine Newteq Performs Newton Iterations With A Line
*  Search, To Solve The Nonlinear Equations.


*  Set Up Constants If This Is The First Call To Newteq.

      If (Frscal) Then
         Frscal = .false.
         Epsaf = Epsmch
         Epsag = Epsmch
         Tolabs = Epsmch
         Tolrel = Epsmch
         Toltny = Epsmch
      Endif
      Ninter = Nmsh - 1

      If (Iprint .eq. 1) Write(6,901) 
  
*  A Newton Method With Line Search And Watchdog Safeguarding
*  Is Used To Try To Solve The Nonlinear Equations.

*  Initialize Iter (The Counter Of Newton Iterations) And Alfold
*  (The Step Taken At The Previous Iteration).

      Iter = -1
      Alfold = One
      Alfa = Zero
      Rnbest = Flmax

      If (.not. Rhsgiv) Then

*  If Necessary, Evaluate The Right-Hand Side At The Initial U. 

         Call Rhscal (Ncomp, Nmsh, Nlbc, Xx, Nudim, U, Defcor,
     *      Fsub, Gsub, Rhs, Rnsq, Fval, Ftmp, Uint, Eps) 
      Endif

*  At Any Given Newton Iteration, Rnprev Is The Value Of Rnsq At 
*  The Immediately Preceding Newton Iteration.

      Rnprev = Flmax

      If (Iprint .eq. 1) Write(6,902)

*  Initialize Counter Of Watchdog Iterations.

      Itwtch = 0

*  Statement 100 Is The Top Of The Newton Iteration Loop.

  100 Continue

      Iter = Iter + 1

      If (Iprint .eq. 1) Write(6,910) Iter

*  If There Have Been Too Many Newton Iterations, Terminate.

      If (Iter .ge. Lmtnwt) Then
         If (Iprint .eq. 1) Write(6,903)
         Iflag = -2
         Return
      Endif

*  The Vector Rhs Is The Right-Hand Side At The Current Iterate,
*  And Rnsq Is Its Squared Two-Norm.
*  Perform Watchdog Tests, Using The Unscaled Merit Function (Rnsq)
*  As The Watchdog Function.  The Routine Wtchdg Updates Rnbest
*  And Itwtch.  If Iflwat Is Not Zero, This Sequence Of Newton
*  Iterations Is Terminated.
                                                                       
      Call Wtchdg ( Iter, Rnsq, Rnbest, Rnprev, Itwtch, 
     *                Alfold, Iflwat )
 
      If (Iflwat .ne. 0) Then
         If (Iprint .eq. 1) Write(6,904) Iter
         Iflag = -3
         Return
      Endif

*  Watchdog Tests Are Passed.  Proceed With The Newton Iteration.
*  Call Jaccal To Evaluate The Block-Structured Jacobian Matrix, 
*  Which Is Stored In Three Matrices (Topblk, Ajac, And Botblk).
*  The Matrices Bhold And Chold Are Saved For Use In Later 
*  Calculations In The Outer Routine.    
      

*  If Rnsq Is Sufficiently Small, Terminate Immediately.
*  Note That The Stored Jacobian Does Not Correspond Exactly
*  To The Final Point.

      If (Rnsq .le. Epsmch) Then
         If (Iprint .eq. 1)  Write(6,906) Iter, Rnsq
         Iflag = 0
         Return
      Endif

      Call Jaccal (Ncomp, Nmsh, Nlbc, 
     *    Xx, Nudim, U, Fval, Dgtm, Dftmp1, Dftmp2, Uint,
     *    Ajac, Topblk, Botblk, Bhold, Chold,
     *    Dfsub, Dgsub, Eps)

*  Blkdcm Is Called To Calculate The Lu Factors Of The Jacobian, 
*  Which Are Overwritten On Topblk, Ajac And Botblk.
*  Interchanges Are Represented In The Integer Array Ipivot.  


      Call Blkdcm (Nlbc, Ncomp-Nlbc, Ncomp, Ninter, 
     *    Topblk, Ajac, Botblk, Ipivot, Iflfac)


*  If The Jacobian Is Found To Be Singular, Iflfac Will Be -1 On 
*  Exit From Blkdcm.

      If (Iflfac .eq. (-1)) Then
         If (Iprint .eq. 1) Write(6,905) Iter
         Iflag = -1 
         Return
      Endif

*   Solve For The Newton Step Delu.  Copy The Rhs Array Into Tmprhs,
*   Which Is Then Overwritten By Blkslv.

      Call Dcopy(Ncomp*Nmsh, Rhs, 1, Tmprhs, 1)  
      Call Blkslv (Nlbc, Ncomp-Nlbc, Ncomp, Ninter, 
     *       Topblk, Ajac, Botblk, Ipivot, Tmprhs, Delu)


*  If Imerit = 1, The Line Search Is Based On The Scaled Merit Function,
*  The Squared Two-Norm Of The Solution Xmerit Of The Linear System
*      (Jacobian)*Xmerit = Rhs, 
*  Where (Jacobian) Is The Jacobian At The Current Newton Iterate. 
*  Thus The Initial Value Of The Scaled Merit Function Is Simply
*  The Squared Two-Norm Of The Newton Step Delu Itself. 

      If (Imerit .eq. 1) Then
         Call Mssq( Ncomp, Nmsh, Delu, Xmscal, Xmsq )
         Fmtry = (Xmscal**2)*Xmsq
      Else
*  The Unscaled Merit Function Is Simply The Squared Two-Norm Of Rhs. 
         Fmtry = Rnsq
      End If

C  Fa And Oldg Represent The Merit Function And Its Gradient
C  At The Initial Point Of The Line Search.

      Fa = Fmtry
      Oldg = -Two*Fa
      Alfa = Zero
      If (Iprint .eq. 1) Write(6,908) Alfa, Fmtry, Rnsq

*  On The First Newton Iteration, The Initial Trial Step Is Unity.   
*  On Subsequent Iterations, The Initial Step Is Not Allowed To 
*  Be More Than The Factor Stpfct Larger Than The Final Step At 
*  The Immediately Preceding Iteration.
     
      Alfa = One
      If(Stpfct*Alfold .lt. One) Alfa = Stpfct*Alfold

      If (Alfa .lt. Alfsml) Alfa = Alfsml

      Fmold = Fa
      Inform = -1

*  Statement 150 Is The Top Of The Inner Line Search Iteration.
*  The Line Search Routine Getptq Has Been Altered So That It
*  Terminates With An Indication Of Success As Soon As A
*  Strictly Lower Value Of The Merit Function Is Found.  Note That
*  This Is A Much Less Strict Requirement Than The Usual Sufficient
*  Decrease Conditions.

      
  150 Continue
      Iwr = 6
      Call Getptq (Gtpdeb, Mfsrch, Iwr, Alfmax, Alfsml, Alfuzz,
     *      Epsaf, Epsag, 
     *      Fmtry, Fmold, Oldg, Rmu, Tolabs, Tolrel, Toltny,
     *      Imprvd, Inform, Nfsrch, Alfa, Alfbst, Fbest, 
     *      Braktd, Crampd, Extrap, Vset, Wset, Nsamea, Nsameb,
     *      Alin, Blin, Fa, Factor, Fv, Fw, Xtry, Xv, Xw)

*  Inform = 1, 2 Or 3 Indicates Success In Finding An Acceptable Point.
*  Inform = 4 Means Alfmax Is Too Small (This Should Never Happen Here,
*  Since Alfmax Is Set Always To 1.1).
*  Inform = 5 Means That A Decrease Was Not Achieved For Any Step
*  Greater Than Alfsml.
*  Inform = 6 Means A Better Point Could Not Be Found (The Minimum
*  Probably Lies Too Close To Alfa = 0).
*  Inform = 7 Means That The Gradient At Alfa = 0 (Oldg) Is Positive 
*  (This Cannot Happen Here, Since Oldg = -Two*Fa, And Fa Is A Non-Negative
*  Number)

      If (Inform .eq. 5) Then
          Iflag = -5
          Return
      Elseif (Inform .eq. 4 .or. Inform .eq. 7) Then
         Iflag = -4
         Return
      Elseif (Inform .eq. 0) Then

*  Inform = 0 Means That A New Function Value Should Be Obtained
*  With The Step Alfa.
*  We May Override Alfa From Getptq By Requiring That The Step Is Not 
*  Allowed To Decrease By More Than A Factor Of Shrfct During 
*  A Line Search Iteration.
*
         If (Alfa .lt. Alfold/Shrfct) Alfa = Alfold/Shrfct
         Alfold = Alfa

*  Define The Next Iterate Utrial = U + Alfa*Delu. 
*  Call Fneval And Rhscal To Evaluate The Right-Hand Side 
*  Rhstri At Utrial. 
*  The Vector Rhstri Is Stored Separately, And Rhs Is Overwritten
*  Only When An Improved Point Is Found.

         Call Matcop ( Nudim, Ncomp, Ncomp, Nmsh, U, Utrial)
         Call Maxpy ( Ncomp, Nmsh, Alfa, Delu, Ncomp, Utrial )
         Call Fneval(Ncomp, Nmsh, Xx, Ncomp, Utrial, Fval, Eps)
         Call Rhscal (Ncomp, Nmsh, Nlbc, Xx, Ncomp, Utrial, Defcor,
     *      Fsub, Gsub, Rhstri, Rnsqtr, Fval, Ftmp, Uint, Eps) 

         Fmold = Fmtry
         If (Imerit .eq. 1) Then

*  Solve A Linear System To Obtain The 2-D Array Xmerit Whose Squared 
*  Norm Is The Scaled Merit Function.   The Lu Factors Of The Jacobian 
*  Have Already Been Calculated By Blkdcm.  
*  Copy Rhstri Into Tmprhs, Which Is Overwritten By Blkslv. 

            Call Dcopy(Ncomp*Nmsh, Rhstri, 1, Tmprhs, 1) 
            Call Blkslv (Nlbc, Ncomp-Nlbc, Ncomp, Ninter,
     *        Topblk, Ajac, Botblk, Ipivot,
     *        Tmprhs, Xmerit)
            Call Mssq( Ncomp, Nmsh, Xmerit, Xscale, Xsolsq )
            Fmtry = (Xscale**2)*Xsolsq
         Else

*  The Unscaled Merit Function Is The Squared Two-Norm Of The Right-Hand
*  Side.
            Fmtry = Rnsqtr
         End If
         If (Iprint .eq. 1) Write(6,908) Alfa, Fmtry, Rnsqtr
         Go To 150
      Endif

*  To Reach Here, Inform Must Be 1, 2, 3, Or 6, And The Line Search 
*  Has Found A Strictly Lower Value Of The Merit Function.
*  Store The New Newton Iterate In U, And The Corresponding Rhs
*  Vector In Rhs.
      
      Rnprev = Rnsq
      Rnsq = Rnsqtr
      Call Matcop (Ncomp, Nudim, Ncomp, Nmsh, Utrial, U)
      Call Dcopy(Ncomp*Nmsh, Rhstri, 1, Rhs, 1)
      If (Iprint .eq. 1) Write(6,909) Iter, Alfa, Fmtry, Rnsq

*  Now Test For Convergence Using The Ratio Of The Newton Step 
*  For Each Component With Max(1, Abs(Current Solution Estimate)).
*  If The Test Fails For Any Element Of U, Branch Back To The
*  Top Of The Newton Iteration.

      Do 160 Im = 1, Nmsh
      Do 160 It = 1, Ntol
         Icmp = Ltol(It)
         Er = Abs(Delu(Icmp,Im))/Max(Abs(U(Icmp,Im)), One)
         If (Er .gt. Cnvfct*Tol(It)) Go To 100
  160 Continue
  
      If (Iprint .eq. 1) Write(6, 906) Iter+1, Rnsq

      Iflag = 0

*  To Fall Through The Above Loop, The Termination Test For A 
*  Sufficiently Small Delu Is Satisfied. 
*  Note That The Stored Jacobian And Its Factorization Do Not
*  Correspond To The Final Solution. 

      Return

  901 Format(1h ,'Start Newton Iterations')
  902 Format(1h ,' Iter',
     *              7x,'Alfa',6x,'Merit',7x,'Rnsq')
  903 Format(1h ,'Too Many Newton Iterations')
  904 Format(1h ,'Watchdog Tests Fail, Iter  = ', I5)
  905 Format(1h ,'Singular Jacobian, Iter = ',I5)
  906 Format(1h ,'Convergence, Iter  = ',I5,4x,'Rnsq  = ',1pe12.3,/)
  908 Format(1h ,'Alfa, Merit, Rnsq',3(1pe11.3))
  909 Format(1h ,I5,3(1pe11.3))
  910 Format(1h ,'Newton Iteration',I5)
      End 


      Subroutine Wtchdg ( Iter, Wmerit, Wmbest, Wmprev, 
     *      Itwtch, Alfold, Iflag )

*  Logic For Watchdog Tests.

      Implicit Double Precision (A-H,O-Z)
      Parameter ( Itonew = 5, Itwtmx = 8, Grfct = 100.0d+0 )
      Parameter ( Half = 0.5d+0 )

*  Perform Watchdog Tests In Two Forms: 
*  (1) To Determine Whether A Sufficient Decrease In The 
*  Watchdog Merit Function Has Occurred Within The Most Recent
*  Sequence Of Itwtmx Iterations;
*  (2) To Determine Whether The Watchdog Merit Function Has Increased 
*  Too Much In A Single Iteration After Itonew Newton Iterations 
*  Have Been Performed.  This Allows The Merit Function To Increase
*  Wildly Only During The First Itonew Iterations.

*  Wmbest Is The Smallest Watchdog Merit Function Achieved In This 
*  Sequence Of Newton Iterations.
*  Wmprev Is The Watchdog Merit Function From The Immediately 
*  Preceding Newton Iteration.

*  Itwtch Counts The Number Of Iterations Without An Improvement
*  In The Unscaled Merit Function.

      Iflag = 0      
      If (Wmerit .le. Wmbest) Then 

*  The Current Watchdog Merit Function Is The Best.

         Wmbest = Wmerit
         Itwtch = 0 
         Return
      Endif

*  The Current Merit Function Is Not The Best.

      Itwtch = Itwtch + 1

*  Do Not Apply Watchdog Tests If (1) The Previous Step Alfold 
*  Exceeds 1/2, Or (2) The Watchdog Merit Function Decreased In 
*  The Immediately Preceding Iteration And Itwtch Does Not
*  Exceed Twice Its Maximum.

      If (Alfold .ge. Half) Return
      If (Wmerit .le. Wmprev .and. Itwtch .le. 2*Itwtmx) Return


*  If More Than Itwtmx Iterations Have Occurred Without 
*  An Overall Improvement In The Watchdog Merit Function,
*  Signal For Termination.

      If (Itwtch .ge. Itwtmx) Then
         Iflag = -1

*  If A Too-Large Increase In The Watchdog Merit Function 
*  Compared To The Best Value Occurred, And Iter .ge. Itonew,
*  Signal For Termination.

      Elseif (Iter .ge. Itonew .and. 
     *          Wmerit .gt. Grfct*Wmbest) Then
          Iflag = -1
      Endif
      Return
      End
* File Evals.f

      Subroutine Fneval(Ncomp, Nmsh, Xx, Nudim, U, Fval, Eps)
      Implicit Double Precision (A-H,O-Z)
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp,Nmsh)

*  Fneval Evaluates The Function Values (From Fsub) For
*  A Given Mesh Xx And Array U, And Stores The Values
*  In The Array Fval.

      Call Fsub (Ncomp, Xx(1), U(1,1), Fval(1,1), Eps)
      Do 50 Im = 1, Nmsh-1
         Call Fsub (Ncomp, Xx(Im+1), U(1,Im+1), Fval(1,Im+1), Eps)
   50 Continue
      Return
      End


      Subroutine Jaccal (Ncomp, Nmsh, Nlbc, 
     *   Xx, Nudim, U, Fval,
     *   Dgtm, Dftm1, Dftm2, Uint,
     *   Ajac, Topblk, Botblk, Bhold, Chold,
     *   Dfsub, Dgsub, Eps)

      Implicit Double Precision (A-H,O-Z)
                                                                  
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp,Nmsh)
      Dimension Dgtm(Ncomp)
      Dimension Dftm1(Ncomp, Ncomp), Dftm2(Ncomp, Ncomp), 
     *             Uint(Ncomp)
      Dimension Ajac(Ncomp, 2*Ncomp, Nmsh-1)
      Dimension Topblk(Ncomp, Ncomp), Botblk(Ncomp,Ncomp)
      Dimension Bhold(Ncomp, Ncomp, Nmsh-1), 
     *             Chold(Ncomp, Ncomp, Nmsh-1)

      External  Dfsub, Dgsub 

*  Blas: Dcopy, Ddot

      Parameter ( Half = 0.5d+0, Eighth = 0.125d+0 )
      Parameter ( Six = 6.0d+0 )
      Parameter ( One = 1.0d+0, Three = 3.0d+0, Twelve = 12.0d+0 )



      Ninter = Nmsh - 1
      
      Do 110 I = 1, Nlbc
         Call Dgsub (I, Ncomp, U(1,1), Dgtm, Eps)
         Call Dcopy(Ncomp, Dgtm(1), 1, Topblk(I,1), Ncomp)
  110 Continue

      Call Dfsub (Ncomp, Xx(1), U(1,1), Dftm1(1,1), Eps)

*  On Entry To Jaccal, The Array Fval Contains The Function Values
*  At (Xx(Im), U(Ic,Im)), Ic = 1,...,Ncomp And Im = 1,...,Nmsh,
*  Calculated By A Preceding Call Of Rhscal With The Same Xx And U
*  Arrays.

      Do 200 Im = 1, Ninter

         Hmsh = Xx(Im+1) - Xx(Im)

         Do 120 Ic = 1, Ncomp
            Uint(Ic) = Half*(U(Ic,Im) + U(Ic,Im+1))  
     *         - Eighth*Hmsh*(Fval(Ic,Im+1) - Fval(Ic,Im))
  120    Continue
         Xhalf = Half*(Xx(Im+1) + Xx(Im))
         Call Dfsub (Ncomp, Xhalf, Uint, Dftm2(1,1), Eps)
         Do 140 Ic = 1, Ncomp
            Do 130 Jc = 1, Ncomp
               Dsq = Ddot(Ncomp, Dftm2(Ic,1), Ncomp,
     *                       Dftm1(1,Jc), 1)
               Ajac(Ic,Jc,Im) = -Hmsh*(Dftm1(Ic,Jc)/Six
     *             + Dftm2(Ic,Jc)/Three + Hmsh*Dsq/Twelve)
  130       Continue
            Ajac(Ic,Ic,Im) = Ajac(Ic,Ic,Im) - One
  140    Continue

         Call Dfsub (Ncomp, Xx(Im+1), U(1,Im+1), Dftm1(1,1), Eps)
         Do 170 Ic = 1, Ncomp
            Do 160 Jc = 1, Ncomp
               Dsq = Ddot(Ncomp, Dftm2(Ic,1), Ncomp,
     *                         Dftm1(1,Jc), 1)
               Ajac(Ic,Jc+Ncomp,Im) = -Hmsh*(Dftm1(Ic,Jc)/Six
     *               + Dftm2(Ic,Jc)/Three - Hmsh*Dsq/Twelve)
  160       Continue
            Call Dcopy(Ncomp, Ajac(Ic,Ncomp+1,Im), Ncomp,
     *                   Chold(Ic,1,Im), Ncomp)
            Call Dcopy(Ncomp, Dftm1(Ic,1), Ncomp,
     *                   Bhold(Ic,1,Im), Ncomp)
            Ajac(Ic,Ic+Ncomp,Im) = Ajac(Ic,Ic+Ncomp,Im) + One
            Chold(Ic,Ic,Im) = Ajac(Ic,Ic+Ncomp,Im)
  170    Continue


  200 Continue
      Do 220 I = Nlbc+1, Ncomp
         Call Dgsub (I, Ncomp, U(1, Nmsh), Dgtm, Eps) 
         Call Dcopy(Ncomp, Dgtm(1), 1, Botblk(I-Nlbc,1), Ncomp)
  220 Continue

      Return

      End 



      Subroutine Lnrhs (Ncomp, Nmsh, Nlbc,
     *   Xx, Nudim, U, Fsub, Gsub, 
     *   Rhs, Rnsq, Fval, Ftmp, Uint, Eps) 

       Implicit Double Precision(A-H,O-Z)

*  This Subroutine Is Designed To Calculate The Right-Hand
*  Side For Linear Problems.

      Dimension Xx(*), U(Nudim,*)
      Dimension Rhs(*), Fval(Ncomp,*), Ftmp(*), Uint(*)
      External Fsub, Gsub 

      Parameter ( Zero = 0.0d+0, Half = 0.5d+0, Eighth = 0.125d+0 )
      Parameter ( Four = 4.0d+0, Six = 6.0d+0 )

      Intrinsic Abs

*  Blas: Dssq

      Ninter = Nmsh - 1
      Rnsq = Zero

*  First, Process The Left-Hand Boundary Conditions.

      Do 20 I = 1, Nlbc
         Call Gsub (I, Ncomp, U(1,1), Wg, Eps)
         Rhs(I) = -Wg
   20 Continue

*  Next, Process The Interior Mesh Points.  The Fval Array
*  Contains The Function Values From Fsub At Xx And U.

      Do 50 Im = 1, Ninter
         Hmsh = Xx(Im+1) - Xx(Im)
         Do 30 Ic = 1, Ncomp
            Uint(Ic) = Half*(U(Ic,Im) + U(Ic,Im+1))  
     *         - Eighth*Hmsh*(Fval(Ic,Im+1) - Fval(Ic,Im))
   30    Continue
         Xhalf = Half*(Xx(Im) + Xx(Im+1)) 
         Call Fsub (Ncomp, Xhalf, Uint, Ftmp, Eps)
         Loc = (Im-1)*Ncomp + Nlbc
         Do 40 Ic = 1, Ncomp
            Rhs(Loc+Ic) = -U(Ic,Im+1) + U(Ic,Im)
     *        + Hmsh*
     *            (Fval(Ic,Im) + Fval(Ic,Im+1) + Four*Ftmp(Ic))/Six
   40    Continue
   50 Continue

      Nrhs = Ninter*Ncomp
      Do 60 Ii = Nlbc+1, Ncomp
         Call Gsub (Ii, Ncomp, U(1,Nmsh), Wg, Eps) 
         Rhs(Nrhs+Ii) = -Wg
   60 Continue

      Call Dssq  ( Nmsh*Ncomp, Rhs, 1, Scale, Sumsq )
      Rnsq = (Scale**2)*Sumsq

      Return
      End


      Subroutine Rhscal (Ncomp, Nmsh, Nlbc,
     *   Xx, Nudim, U, Defcor,
     *   Fsub, Gsub, 
     *   Rhs, Rnsq, Fval, Ftmp, Uint, Eps) 

       Implicit Double Precision(A-H,O-Z)

*  This Subroutine Constructs The (Ncomp*Nmsh)-Dimensional
*  Vector Rhs, Which Is The Right-Hand Side Of The Newton Equations.
*  The Ncomp By Nmsh Array Fval Is Assumed To Have Been Calculated
*  Elsewhere By Routine Fneval.

      Dimension  Xx(Nmsh), U(Nudim,Nmsh), Defcor(Ncomp,Nmsh-1) 
      Dimension  Rhs(Ncomp*Nmsh), Fval(Ncomp,Nmsh)
      Dimension  Ftmp(Ncomp), Uint(Ncomp)
      External   Fsub, Gsub 

      Intrinsic Abs

*  Blas: Dssq

      Parameter ( Zero = 0.0d+0, Half = 0.5d+0, Eighth = 0.125d+0 )
      Parameter ( Four = 4.0d+0, Six = 6.0d+0 )

*  Ninter Is The Number Of Intervals In The Mesh (One Less Than The
*  Number Of Mesh Points)

      Ninter = Nmsh - 1
      Rnsq = Zero

*  First, Process The Left-Hand Boundary Conditions.

      Do 20 I = 1, Nlbc
         Call Gsub (I, Ncomp, U(1,1), Wg, Eps)
         Rhs(I) = -Wg
   20 Continue

*  Next, Process The Interior Mesh Points.  The Fval Array
*  Contains The Function Values From Fsub At Xx And U.

      Do 50 Im = 1, Ninter
         Hmsh = Xx(Im+1) - Xx(Im)
         Do 30 Ic = 1, Ncomp
            Uint(Ic) = Half*(U(Ic,Im) + U(Ic,Im+1))  
     *         - Eighth*Hmsh*(Fval(Ic,Im+1) - Fval(Ic,Im))
   30    Continue
         Xhalf = Half*(Xx(Im) + Xx(Im+1)) 
         Call Fsub (Ncomp, Xhalf, Uint, Ftmp, Eps)
         Loc = (Im-1)*Ncomp + Nlbc
         Do 40 Ic = 1, Ncomp
            Rhs(Loc+Ic) = -U(Ic,Im+1) + U(Ic,Im) + Defcor(Ic,Im)
     *        + Hmsh*
     *            (Fval(Ic,Im) + Fval(Ic,Im+1) + Four*Ftmp(Ic))/Six

   40    Continue
   50 Continue

      Nrhs = Ninter*Ncomp
      Do 60 Ii = Nlbc+1, Ncomp
         Call Gsub (Ii, Ncomp, U(1,Nmsh), Wg, Eps) 
         Rhs(Nrhs+Ii) = -Wg
   60 Continue

      Call Dssq  ( Nmsh*Ncomp, Rhs, 1, Scale, Sumsq )
      Rnsq = (Scale**2)*Sumsq

      Return

      End
* File Msh.f

      Subroutine Dblmsh (Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh) 
      Implicit Double Precision (A-H,O-Z)
      Dimension Xx(*), Xxold(*)
      Logical Maxmsh

      Common/Algprs/ Nminit, Iprint, Maxcon, Itsaim, Uval0
      Common /Flags/ Ifinal,Iatt,Iback,Iprec

*  Blas: Dcopy

      Parameter (Half = 0.5d+0)

*  This Routine Is Used To Double The Mesh, I.e., Produce A Mesh
*  With Twice As Many Intervals In Which Each New Interval Is
*  Half The Corresponding Old Interval.

*  On Entry To Dblmsh, The Integer Nmsh And The Array Xx
*  Specify A Set Of Mesh Points Xx(1),..., Xx(Nmsh) (Assumed 
*  To Be In Ascending Order).

*  If The Number Of Mesh Points In The Doubled Mesh Would
*  Exceed The Maximum Allowed Number Nmax, The Flag Maxmsh Is
*  Set To True, And We Exit Without Changing Any Other Parameters.

*  Otherwise, Nmold Is Set To The Old Number Of Mesh Points,
*  Xxold Is Set To The Old Mesh, Nmsh Is The New Number Of Mesh
*  Points, And Xx Contains The New Mesh Points.

      Nmold = Nmsh
      Iprec = Min(Iprec,1)
      Call Dcopy(Nmold, Xx, 1, Xxold, 1)

      Ninnew = 2*(Nmsh-1)
      Nmnew = Ninnew + 1
      If(Nmnew .ge. Nmax) Then
         If (Iprint .ge. 0)  Write(6,901) Nmnew
         Maxmsh = .true.
         Return
      Endif
      Maxmsh = .false.

*  Loop Backwards Through The Old Mesh Points To Create The New Ones.

      Xx(Nmnew) = Xx(Nmsh)
      Do 100 I = Ninnew, 4, -2
         Id2 = I/2
         Xx(I) = Half*(Xx(I+1) + Xx(Id2))
         Xx(I-1) = Xx(Id2)
         If (Xx(I) .eq. Xx(I+1) .or. Xx(I) .eq. Xx(I-1)) Then
            Iprec = 2
            Return
         Endif
  100 Continue

*  Calculate The New Xx(2). Xx(1) Remains Unchanged.

      Xx(2) = Half*(Xx(3) + Xx(1))
      If (Xx(2) .eq. Xx(3) .or. Xx(2) .eq. Xx(1)) Then
         Iprec = 2
         Return
      Endif
      Nmsh = Nmnew
      Return
  901 Format (1h , ' Dblmsh.  Maximum Mesh Exceeded, Nmnew  = ', I8)
      End 


      Subroutine Selmsh(Ncomp, Nmsh, Ntol, Ltol, Tol,
     *     Nfxpnt, Fixpnt, Ipow, Nmax, Nvold,
     *     Xx, Nudim, U, Ermeas, Irefin, Ihcomp,
     *     Nmold, Xxold, Ermx, Double, Maxmsh, Phiold, Voldmsh)

      Implicit Double Precision (A-H,O-Z)

      Dimension  Ltol(Ntol), Tol(Ntol), Fixpnt(*)
      Dimension  Xx(*), U(Nudim, *), Ermeas(Ncomp,*)
      Dimension  Irefin(Nmsh-1), Ihcomp(Nmsh-1)
      Dimension  Xxold(*), Ermx(*), Phiold(*), Voldmsh(*)
      Dimension  Npr(2), Pmax(2), Hord(2)
      Logical    Double, Maxmsh
      
      Intrinsic Abs, Max, Int

*  Blas: Dcopy

      Common /Flags/ Ifinal,Iatt,Iback,Iprec
      Common /Mshvar/ Hsml,Npr,Pmax,Hord

      Parameter  ( Zero = 0.0d+0, One = 1.0d+0, Onep1 = 1.1d+0 ) 
      Parameter  ( Erdcid = 5.0d+0 )
      Parameter  ( Phitst = 0.1d+0 )
  
*  The Routine Selmsh Performs Selective Mesh Refinement, Depending
*  On The Error Measure Ermeas.
   
      Maxmsh = .false.

      Frcpow = One/Ipow
      Double = .false.
      Nmold = Nmsh
      Ninter = Nmsh - 1
      Iprec = Min(Iprec,1)

*  Copy The Current Mesh Into The Xxold Array.

      If (Iatt .eq. -1) Nvold = Nmsh
      If (Iatt .eq. 0) Call Dcopy(Nvold, Xxold, 1, Voldmsh, 1)

      Call Dcopy(Nmold, Xx, 1, Xxold, 1)
      Ithres = 0
      Thres = One

*  On Input, The Array Ermeas Represents Some Error Measure Defined
*  Over The Components And Mesh Intervals (Not Mesh Points).  
*  It Is Normalized In The Following Loop With Respect To The
*  Tolerance Array And The Current Solution.
*  The Value Errmax Gives The Maximum Normalized Error.

      Errmax = Zero
      Do 120 Im = 1, Ninter
         Ermx(Im) = Zero
         Do 110 It = 1, Ntol
            Jcomp = Ltol(It)
            Umin = Min(Abs(U(Jcomp,Im)),Abs(U(Jcomp,Im+1)))
            Denom = Tol(It)*Max(One, Umin)
            Ems = Ermeas(Jcomp,Im)
            Ermeas(Jcomp,Im) = Abs(Ems)/Denom
            Err = Ermeas(Jcomp, Im)
            If (Err .ge. Ermx(Im)) Then
                Ermx(Im) = Err
                Ihcomp(Im) = Jcomp
            Endif
  110    Continue
         Errmax = Max(Ermx(Im), Errmax)
  120 Continue

  200 Continue

*  For Each Interval Im,  The Integer Irefin(Im) Is Calculated 
*  Based On An Equidistribution Procedure Involving The
*  Threshold Value Thres.  If Irefin(Im) > 1, We Add 
*  Points To Interval Im.  If Irefin(Im) = 1, We Check Whether
*  Point Im Can Be Removed From The Mesh.

*  Nmest Is A Lower Bound On The Number Of Points In The New Mesh.
*  We Do Not Know In Advance How Many Points Will Be Removed,
*  So Nmest Is Computed By Assuming All Eligible Points Are Removed.

      Philrg = 0.d0
      Esum = 0.d0
      Nmest = Nmsh
      Mshchng = 0
      Do 220 Im = 1, Ninter
         Errim = Ermx(Im)**Frcpow
         If (Ermx(Im) .ge. Thres) Then
            Irefin(Im) = Int(Errim) + 1
            Nmest = Nmest + Irefin(Im) - 1
            Mshchng = 1
         Else
            Irefin(Im) = 1
            Nmest = Nmest - 1
         Endif
         Him = Xxold(Im+1)-Xxold(Im)
         Phiim = Errim/Him
         Esum = Esum+Errim
         If (Iatt .eq. -1) Then
            Phiold(Im) = Phiim
         Elseif (Iatt .eq. 0 .and. Phiim .gt. Philrg) Then
           Philrg = Phiim
           Hordlrg = Errim
           Imreg = Im
         Endif
  220 Continue
C
      If (Iatt .eq. 0) Then
        Hord(2) = Hordlrg
        Pmax(2) = Philrg
        Xloc1 = Xxold(Imreg)
        Xloc2 = Xxold(Imreg+1)
        Ichkpt = Nvold/2
        If (Xloc1 .lt. Voldmsh(Ichkpt)) Ichkpt = 1
        Xb = Voldmsh(Ichkpt)
        Do 225 J = Ichkpt, Nvold-1
          Xa = Xb
          Xb = Voldmsh(J+1)
          If (Xloc1 .ge. Xa .and. Xloc1 .lt. Xb) Then
            If (Xloc2-Xb .lt. (Xloc2-Xloc1)/2.d0) Then
              Pmax(1) = Phiold(J)
              Hord(1) = (Xb-Xa)*Pmax(1)
            Else
              Pmax(1) = Phiold(J+1)
              Hord(1) = (Voldmsh(J+2)-Xb)*Pmax(1)
            Endif
            Goto 228
          Endif
 225    Continue
 228    Hsml = (Xloc2-Xloc1)/Dble(Irefin(Imreg))
        Npr(2) = Int(Esum)
      Endif

      If (Iatt .eq. -1) Npr(1) = Int(Esum)
      If (Ifinal .eq. 1 .and. Mshchng .eq. 0) Iatt = 0
      If (Nmest .gt. Nmax) Go To 360

*  It Appears That We Can Perform The Desired Selective Mesh
*  Refinement.

*  Now Begin Running Through The Mesh, Adding And Possibly Deleting
*  Points As Indicated By The Irefin Array.

*  The Integer New Is A Count Of The Number Of Intervals In
*  The Tentative Mesh Being Generated By The Refinement Strategy.

      New = 1

*  The First Interval Is Treated As A Special Case, Since Xx(1)
*  Always Remains In The Mesh, And Cannot Be A Fixed Point.

      Rlen = Xxold(2) - Xx(1)
      Slen = Rlen
      If (Irefin(1) .gt. 1) Then
         Dx = Rlen/Irefin(1)
         Do 230 J = 2, Irefin(1) 
            New = New + 1
            Xx(New) = Xx(1) + (J-1)*Dx
            If (Xx(New) .eq. Xx(New-1)) Then
               Iprec = 2
               Return
            Endif
  230    Continue
      Endif

*  The Fixed Points Specified By The Fixpnt Array Cannot Be
*  Removed From The Mesh.  The Value Fxnext Indicates The 'Next'
*  Fixed Point.  When No Further Fixed Points Remain To Be Processed
*  (Or If Nfxpnt = 0), Fxnext Is Set To A Value Strictly Larger Than
*  The Last Mesh Point, So That No Mesh Point Can Equal Fxnext.  
*  This Way We Need To Compare Only The Values Of Xxold(I)
*  And Fxnext.

      Ifxcnt = 1
      If (Nfxpnt .eq. 0) Then
         Fxnext = Onep1*Abs(Xxold(Nmsh))
      Else
         Fxnext = Fixpnt(Ifxcnt)
      Endif

*  Jtkout Is A Counter Of The Number Of Consecutive Points That 
*  Have Been Removed From The Mesh.

      Jtkout = 0
      Do 330 Im = 2, Ninter
         Rlold = Rlen
         Rlen = Xxold(Im+1) - Xxold(Im)

*  If Xxold(Im) Is The Next Fixed Point, It Cannot Be Removed
*  And So We Don'T Test Its Error Estimates.

         If(Xxold(Im) .eq. Fxnext)  Then

            Ifxcnt = Ifxcnt + 1
            If(Ifxcnt .gt. Nfxpnt) Then
               Fxnext = Onep1*Abs(Xxold(Nmsh))
            Else
               Fxnext = Fixpnt(Ifxcnt)
            Endif

         Elseif (Irefin(Im) .eq. 1) Then

*  If Xxold(Im) Is Not A Fixed Point And Irefin(Im) = 1, We May Wish 
*  To Remove Point Im From The Mesh.

*  If We Are Considering Removing Points And Jtkout = 0, This
*  Is The First Point In A Possible Consecutive Set To Be Removed,
*  And We Initialize Phihat, Which Represents A Maximum Of
*  Certain Estimates.
*  If Jtkout Is Not Zero, Previous Points Contiguous To This
*  Point Have Been Removed, And Phihat Retains Its Previous Value.

            Slen = Slen + Rlen

            If (Jtkout .eq. 0) Then
                Ind1 = Ihcomp(Im-1)
                Phihat = Ermeas(Ind1,Im-1)/(Rlold**Ipow)
            Endif
            Phihat = Max(Phihat,
     *                 Ermeas(Ihcomp(Im),Im)/(Rlen**Ipow))
            Val1 = Phihat*(Slen**Ipow)
            If (Val1 .le. Phitst 
     *             .and. Jtkout .lt. 4) Then

*  Increment The Counter Of Removed Points.
*  'Remove' The Mesh Point Xxold(Im) By Not Including It.

               Jtkout = Jtkout+1
               Go To 330
            Endif
*        End Of Logic For Irefin(Im) = 1.
         Endif

         Jtkout = 0 
         New = New + 1
         Xx(New) = Xxold(Im)
         If (Xx(New) .eq. Xx(New-1)) Then
            Iprec = 2
            Return
         Endif
         If (Irefin(Im) .gt. 1) Then
            Dx = Rlen/Irefin(Im)
            Do 300 J = 2, Irefin(Im)
               New = New + 1
               Xx(New) = Xxold(Im) + (J-1)*Dx
               If (Xx(New) .eq. Xx(New-1)) Then
                  Iprec = 2
                  Return
               Endif
 300        Continue
         Endif
         Slen = Rlen
         
*  If The New Mesh Contains Too Many Points, Branch Out Of The
*  Loop To Try Alternative Strategies.

         If (New .gt. Nmax-1) Goto 360

  330 Continue

*  To End Up Here, We Have Processed The Entire Interval,
*  And Have Neither Exceeded The Specified Maximum Nor
*  Exceeded Three Times The Number Of Intervals In The Old
*  Mesh.  The Last Mesh Point Remains Unchanged.
     
      New = New + 1
      Xx(New) = Xxold(Nmsh)
      If (Xx(New) .eq. Xx(New-1)) Then
         Iprec = 2
         Return
      Endif
      Nmsh = New
      Maxmsh = .false.
      Return

  360 Continue

*  To Reach Here, The Number Of Mesh Points Created At Some Stage
*  Of The Refinement Process Was Larger Than The Maximum Permitted 
*  Value Nmax.  

*  Check Whether The Mesh Can Safely Be Doubled.
  
      If ((2*Nmsh-1) .lt. Nmax) Then

*  Double The Mesh.
         Call Dblmsh (Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh) 
         Double = .true.

*  If The Number Of Intervals Is Too Large And The Mesh Cannot Be 
*  Doubled, Increase The Threshold Thres By A Factor Of Erdcid And
*  Try The Selective Refinement Again.  
*  If This Happens Three Times Without Success Or If Thres Exceeds
*  Or Is Equal To Errmax, Stop.  (In This Case, We Know Already 
*  That Doubling The Mesh Produces Too Many Points.)

      Elseif (Thres .lt. Errmax .and. Ithres .lt. 3) Then
         Ithres = Ithres + 1
         Thres = Erdcid*Thres
         If(Thres .gt. Errmax) Thres = Errmax
         Call Dcopy(Nmsh, Xxold, 1, Xx, 1)
         Go To 200
      Else
         Nmsh = 2*Nmsh - 1
         Maxmsh = .true.
      Endif
      Return

  904 Format(1h ,'Nmest, Irefin',(10i5))
  910 Format(1h ,'Ihcomp',(10i5))
      End 


      Subroutine Unimsh(Nmsh, Aleft, Aright, Nfxpnt, Fixpnt, Xx)
      Implicit Double Precision (A-H,O-Z)
      Integer  Nmsh, Nfxpnt
      Dimension Fixpnt(*), Xx(Nmsh)

      Intrinsic Max

*  Given A Left Endpoint Aleft, A Right Endpoint Aright,
*  A Set Of Nfxpnt Fixed Points Fixpnt(I), I = 1,...,Nfxpnt,
*  (Where Fixpnt(I) Is Different From Aleft And Aright For All I),
*  And An Initial Target Number Nmsh Of Mesh Points,
*  The Subroutine Unimsh Generates A Piecewise Uniform Mesh 
*  Beginning At Aleft, Ending At Aright, And With Equally 
*  Spaced Points Between Aleft And Fixpnt(1), Then Between 
*  Fixpnt(1) And Fixpnt(2), ..., And Finally Between 
*  Fixpnt(Nfxpnt) And Aright.  The Final Number Of Intervals
*  Is The Maximum Of Nfxpnt+2 And The Initial Value Of Nmsh.

*  In The Simplest Case When Nfxpnt = 0, Unimsh Generates A
*  Uniform Mesh With Nmsh Intervals In The Closed Interval
*  (Aleft, Aright).

*  On Exit, The Integer Nmsh Contains The Number Of Mesh Points
*  (Which Is The Maximum Of The Initial Nmsh And Nfxpnt).
*  The Array Xx (Of Dimension Nmsh) Contains The Mesh Points.

      If (Nfxpnt .eq. 0) Then

*  If There Are No Interior Fixed Points, The Spacing Is Uniform 
*  Throughout The Interval.  Calculate The Spacing Dx
*  And Set Up The Xx Array.

        Ninter = Nmsh - 1

         Dx = (Aright - Aleft)/Ninter
         Do 10 I = 1, Ninter
            Xx(I) = Aleft + (I-1)*Dx
   10    Continue
         Xx(Nmsh) = Aright
         Return
      Endif

*  We Know That There Is At Least One Fixed Point Strictly Between
*  The Endpoints.

      If (Nmsh .lt. Nfxpnt+2)  Nmsh = Nfxpnt + 2
      Ninter = Nmsh - 1
      Xx(1) = Aleft
      Ileft = 1
      Xleft = Aleft
      Totint = Aright - Aleft
      Ndif = Ninter - Nfxpnt
      Do 50 J = 1, Nfxpnt + 1

*  Deal In Turn With The Subintervals Defined By The Interval
*  Boundaries And The Fixed  Points.

         If (J .lt. Nfxpnt+1) Then

*  The J-Th Fixed Point Is Xright.  Calculate Where It Should
*  Fall In The Mesh.

            Xright = Fixpnt(J)
            Nmin = Int(Ninter*(Xright-Aleft)/Totint + 1.5d+0)
            If (Nmin .gt. Ndif+J) Nmin = Ndif + J
            Iright = Max(Ileft+1, Nmin)
         Else
            Xright = Aright
            Iright = Nmsh
         Endif

*  Npt Is The Number Of Equally Spaced Points That Should
*  Lie Strictly Between The (J-1)-Th And J-Th Fixed Points.

         Xx(Iright) = Xright
         Npt = Iright - Ileft - 1
         Dx = (Xright - Xleft)/(Npt + 1)
         Do 30 I = 1, Npt
            Xx(Ileft+I) = Xleft + I*Dx
   30    Continue
         Ileft = Iright
         Xleft = Xright
   50 Continue

      Return
C
      End 


* File Rest.f

      Subroutine Errest (Ncomp, Nmsh, Ntol, Ltol, Tol,
     *   Nudim, U, Uold, Etest, Errok)

      Implicit Double Precision (A-H,O-Z)
      Dimension Ltol(Ntol), Tol(Ntol), U(Nudim,Nmsh), 
     *              Uold(Ncomp,Nmsh), Etest(Ntol)
      Logical Errok

      Intrinsic Abs, Max

      Parameter ( One = 1.0d+0 )

 
*  Given Current And Previous Solutions U And Uold On The Same
*  Mesh, Errest Calculates An Error Measure For Each 
*  Component For Which A Tolerance Is Specified.
*  The Error Measure Is The Usual Relative Error Normalized 
*  By Dividing By The Tolerance.  

*  The Array Etest Specifies The Error Test To Be Applied To Each
*  Error Measure.

*  On Exit, The Logical Flag Errok
*   -- Is .false. If Any Of The Error Measures Exceeds The 
*      Corresponding Value In The Array Etest
*   -- Is .true. If All Error Measures Are Less Than The 
*      Corresponding Values Of Etest.

      Errok = .true.

      Do 10 Im = 1, Nmsh
      Do 10 It = 1, Ntol
         Icmp = Ltol(It)
         Er = U(Icmp,Im) - Uold(Icmp,Im)
         Denom = Max(One, Abs(Uold(Icmp,Im)))
         Errel = Abs(Er/(Tol(It)*Denom)) 
         If (Errel .gt. Etest(It)) Errok = .false.
   10 Continue

      Return
      End 


      Subroutine Getptq( Debug, Mfsrch, Nout, Alfmax, Alfsml, Alfuzz, 
     *                   Epsaf, Epsag, Ftry, Oldf, Oldg,
     *                   Rmu, Tolabs, Tolrel, Toltny, Imprvd,
     *                   Inform, Nfsrch, Alfa, Alfbst, Fbest,
     *                   Braktd, Crampd, Extrap,Vset,Wset,Nsamea,Nsameb,
     *                   A, B, Fa, Factor, Fv, Fw, Xtry, Xv, Xw )

      Implicit Double Precision (A-H,O-Z)
      Logical            Debug, Imprvd
      Logical            Braktd, Crampd, Extrap, Vset, Wset 
      Integer            Mfsrch, Nout, Inform, Nfsrch, Nsamea, Nsameb 
C
C  *********************************************************************
C  Getptq  Is A Step-Length Algorithm For Minimizing A Function Of One
C  Variable.  It Will Be Called Repeatedly By A Search Routine Whose
C  Purpose Is To Estimate A Point  Alfa = Alfbst  That Minimizes Some 
C  Function  F(Alfa)  Over The Closed Interval (0, Alfmax). 
C
C  Getptq  Requires The Function  F(Alfa)  (But Not Its Gradient)
C  To Be Evaluated At Various Points Within The Interval.  New
C  Step-Length Estimates Are Computed Using Quadratic Interpolation With
C  Safeguards.
C
C  Reverse Communication Is Used To Allow The Calling Program To
C  Evaluate  F.  Some Of The Parameters Must Be Set Or Tested
C  By The Calling Program.  The Remainder Would Ordinarily Be Local
C  Variables.
C
C
C  Input Parameters (Relevant To The Calling Program)
C  --------------------------------------------------
C
C  Debug         Specifies Whether Detailed Output Is Wanted.
C
C  Inform        Must Be Nonzero On The First Entry (E.g., -1).
C                It Will Be Altered By  Getptq  For Later Entries.
C
C  Mfsrch        Is An Upper Limit On The Number Of Times  Getptq  Is 
C                To Be Entered Consecutively With  Inform = 0
C                (Following An Initial Entry With  Inform Lt 0).
C
C  Nout          Is The File Number To Be Used For Printed Output
C                If Debug Is True.
C
C  Alfa          Is The First Estimate Of The Step Length.  Alfa  Is
C                Subsequently Altered By  Getptq  (See Below).
C
C  Alfmax        Is The Upper Limit Of The Interval To Be Searched.
C
C  Alfsml        Is Intended To Prevent Inefficiency When The Optimum 
C                Step Is Very Small, For Cases Where The Calling
C                Program Would Prefer To Re-Define  F(Alfa).  Alfsml Is
C                Allowed To Be Zero. Early Termination Will Occur If
C                Getptq  Determines That The Optimum Step Lies
C                Somewhere In The Interval  (0, Alfsml)  (But Not If
C                Alfmax .le. Alfsml).
C
C  Epsaf         Is An Estimate Of The Absolute Precision In The
C                Computed Values Of  F. 
C
C  Oldf          Is The Value Of  F(0). 
C
C  Oldg          Is An Estimate Of The Gradient Of  F  At  Alfa = 0.
C                It Should Be Non-Positive.
C
C  Rmu           Controls What Is Meant By A Significant Decrease In  F.
C                The Final  F(Alfbst)  Should Lie On Or Below The Line
C                      L(Alfa)  =  Oldf + Alfa*Rmu*Oldg.
C                Rmu  Should Be In The Open Interval (0, 0.5).
C                The Value  Rmu = 1.0d-4  Is Good For Most Purposes.
C
C  Tolabs,Tolrel Define A Function  Tol(Alfa) = Tolrel*Alfa + Tolabs
C                Such That If  F  Has Already Been Evaluated At Step
C                Alfa,  Then It Will Not Be Evaluated At Any Point
C                Closer Than  Tol(Alfa).
C                These Values May Be Reduced By  Getptq  If They Seem 
C                To Be Too Large.
C
C  Toltny        Is The Smallest Value That  Tolabs  Is Allowed To Be 
C                Reduced To.
C
C
C  Output Parameters (Relevant To The Calling Program)
C  ---------------------------------------------------
C
C  Imprvd        Is True If The Previous Step  Alfa  Was The Best
C                Point So Far.  Any Related Quantities (E.g., Arrays) 
C                Should Be Saved By The Calling Program Before Paying 
C                Attention To  Inform.
C
C  Inform = 0    Means The Calling Program Should Evaluate
C                           Ftry = F(Alfa)
C                For The New Trial Step  Alfa,  And Then Re-Enter
C                Getptq.
C
C  Inform = 1    Means The Search Has Terminated Successfully
C                With A Step  Alfbst  That Is Less Than The 
C                Upper Bound  Alfmax.
C
C  Inform = 2    Means The Search Has Terminated Successfully
C                With A Step  Alfbst  That Is Equal To The
C                Upper Bound  Alfmax.
C
C  Inform = 3    Means That The Search Failed To Find A Point Of
C                Sufficient Decrease In  Mfsrch  Functions, But An
C                Improved Point Was Found.
C
C  Inform = 4    Means  Alfmax  Is So Small That A Search Should
C                Not Have Been Done.
C
C  Inform = 5    Means That The Search Was Terminated Prematurely
C                Because Of The Value Of  Alfsml  (See Above).
C
C  Inform = 6    Means The Search Has Failed To Find A Useful Step.  If
C                The Subroutine For The Function And Gradient Has Been
C                Programmed Correctly, This Will Usually Occur If The 
C                Minimum Lies Very Close To  Alfa = 0  Or The Gradient
C                Is Not Sufficiently Accurate.
C
C  Inform = 7    Means That The Value Of  G(0) Was Positive On Entry. 
C
C  Alfa          Is The Step At Which The Next Function Value Must Be 
C                Computed.
C
C  Alfbst        Should Be Accepted By The Calling Program As The
C                Required Step-Length Estimate, Whenever  Getptq
C                Returns  Inform = 1,  2  Or  3.
C
C  Fbest         Will Be The Corresponding Value Of  F.
C
C
C  The Following Parameters Retain Information Between Entries
C  -----------------------------------------------------------
C
C  Alfuzz        Is Such That, If The Final  Alfa  Lies In The Interval
C                (0,Alfuzz)  And  Abs( F(Alfa)-Oldf ) Le Epsaf,  Alfa 
C                Cannot Be Guaranteed To Be A Point Of Sufficient
C                Decrease.
C
C  Braktd        Is False If  F  Has Not Been Evaluated At The Far End
C                Of The Interval Of Uncertainty.  In This Case, The
C                Point  B  Will Be At  Alfmax + Tol(Alfmax).
C
C  Crampd        Is True If  Alfmax  Is Very Small (Le Tolabs).
C                If The Search Fails, This Indicates That A Zero
C                Step Should Be Taken.
C
C  Extrap        Is True If Alfbst Has Moved At Least Once And  Xv
C                Lies Outside The Interval Of Uncertainty.  In This
C                Case, Extra Safeguards Are Applied To Allow For
C                Instability In The Polynomial Fit.
C
C  Vset          Records Whether A Third-Best Point Has Been
C                Determined.
C
C  Wset          Records Whether A Second-Best Point Has Been
C                Determined.  It Will Always Be True By The 
C                Time The Convergence Test Is Applied (Label 300).
C
C  Nsamea        Is The Number Of Consecutive Times That The Left-Hand
C                End Of The Interval Of Uncertainty Has Remained The
C                Same.
C
C  Nsameb        Similarly For The Right-Hand End.
C
C  A, B, Alfbst  Define The Current Interval Of Uncertainty.
C                The Required Minimum Lies Somewhere Within The
C                Closed Interval  (Alfbst + A, Alfbst + B). 
C
C  Alfbst        Is The Best Point So Far.  It Is Strictly Within The 
C                The Interval Of Uncertainty Except When It Lies At The
C                Left-Hand End When  Alfbst  Has Not Been Moved.
C                Hence We Have    A Le 0,   B Gt 0.
C
C  Fbest         Is The Value Of  F  At The Point  Alfbst.
C
C  Fa            Is The Value Of  F  At The Point  Alfbst + A.
C
C  Factor        Controls The Rate At Which Extrapolated Estimates Of 
C                Alfa  May Expand Into The Interval Of Uncertainty.
C                Factor Is Not Used If The Minimum Has Been Bracketed 
C                (I.e., When The Variable  Braktd  Is True).
C
C  Fv, Fw        Are The Values Of  F  At The Points  Alfbst + Xv,
C                Alfbst + Xw.  They Are Not Defined Until  Vset
C                Or  Wset  (Respectively) Is True.
C
C  Ftry          Is The Value Of  F  At The New Point  Alfbst + Xtry. 
C
C  Xtry          Is The Trial Point Within The Shifted Interval (A, B).
C                The New Trial Function Value Must Be Computed At The 
C                Point  Alfa  =  Alfbst + Xtry.
C
C  Xv            Is Such That  Alfbst + Xv  Is The Third-Best Point.
C                It Is Not Defined Until  Vset  Is True.
C
C  Xw            Is Such That  Alfbst + Xw  Is The Second-Best Point. 
C                It Is Not Defined Until  Wset  Is True.
C                In Some Cases,  Xw  Will Replace A Previous  Xw  That
C                Has A Lower Function But Has Just Been Excluded From 
C                The Interval Of Uncertainty.
C
C
C  Systems Optimization Laboratory, Stanford University, California.
C  Original Version February 1982.  Rev. May 1983.
C  *********************************************************************
C
      Logical            Closef, Conv1, Conv2, Conv3, Convrg
      Logical            Moved, Sigdec, Xinxw
      Data               Zero, Point1, Half/ 0.0d+0,  0.1d+0, 0.5d+0/ 
      Data                One,  Two,  Five   / 1.0d+0, 2.0d+0,  5.0d+0/
      Data                Ten, Eleven      /10.0d+0, 11.0d+0        / 
C
C
C  Local Variables
C  ---------------
C
C  Closef        Is True If The Worst Function  Fv  Is Within  Epsaf
C                Of  Fbest  (Up Or Down).
C
C  Convrg        Will Be Set To True If At Least One Of The Convergence
C                Conditions Holds At  Alfbst.
C
C  Moved         Is True If A Better Point Has Been Found (Alfbst Gt 0).
C
C  Sigdec        Says Whether  Fbest  Represents A Significant Decrease
C                In The Function, Compared To The Initial Value  Oldf.
C
C  Xinxw         Is True If  Xtry  Is In  (Xw,0)  Or  (0,Xw).
C  ---------------------------------------------------------------------
C
      Imprvd = .false.
      If (Inform .ne. -1) Go To 100
C
C  ---------------------------------------------------------------------
C  First Entry.  Initialize Various Quantities, Check Input Data And
C  Prepare To Evaluate The Function At The Initial Step  Alfa.
C  ---------------------------------------------------------------------
      Nfsrch = 0
      Alfbst = Zero 
      Fbest  = Oldf 
      If (Oldg   .gt.      Zero) Go To 970
      If (Oldg   .ge. (- Epsag)) Go To 960
      If (Alfmax .le.    Toltny) Go To 940
C
      Braktd = .false.
      Crampd = Alfmax .le. Tolabs
      Extrap = .false.
      Vset   = .false.
      Wset   = .false.
      Nsamea = 0
      Nsameb = 0
      Alfuzz = Two*Epsaf/(Rmu*Abs( Oldg ))
      A      = Zero 
      B      = Alfmax + (Tolrel*Alfmax + Tolabs)
      Fa     = Oldf 
      Factor = Five 
      Tol    = Tolabs
      Xtry   = Alfa 
      If (Debug) Write (Nout, 1000) Alfmax, Oldf, Oldg, Tolabs,
     *   Alfuzz, Epsaf, Epsag, Tolrel, Crampd
      Go To 800
C
C  ---------------------------------------------------------------------
C  Subsequent Entries.
C  The Function Has Just Been Evaluated At  Alfa = Alfbst + Xtry,
C  Giving  Ftry.
C  ---------------------------------------------------------------------
  100 Nsamea = Nsamea + 1
      Nsameb = Nsameb + 1
      Xtry   = Alfa - Alfbst
      Moved  = Alfbst .gt. Zero
C
C  Check If  Xtry  Is In The Interval  (Xw,0)  Or  (0,Xw).
C
      Xinxw  = .false.
      If (Wset) Xinxw =       Zero .lt. Xtry  .and.  Xtry .le. Xw
     *                  .or.    Xw .le. Xtry  .and.  Xtry .lt. Zero
C
C  See If The New Step Is Better.
C
      Deltaf = Ftry   - Oldf
      Ctry   = Deltaf - Alfa*Rmu*Oldg
      If (Alfa .le. Alfuzz) Sigdec = Deltaf .le. (- Epsaf)
      If (Alfa .gt. Alfuzz) Sigdec = Ctry   .le.    Epsaf
      Imprvd = Sigdec  .and.  ( Ftry - Fbest ) .le. (- Epsaf)
C
      If (Debug) Write (Nout, 1100) Alfa, Ftry, Ctry
      If (.not. Imprvd) Go To 130
C
C  We Seem To Have An Improvement.  The New Point Becomes The
C  Origin And Other Points Are Shifted Accordingly.
C
      If (.not. Wset) Go To 110
      Xv     = Xw - Xtry
      Fv     = Fw
      Vset   = .true.
  110 Xw     = Zero - Xtry
      Fw     = Fbest
      Wset   = .true.
      Fbest  = Ftry 
      Alfbst = Alfa 
      A      =    A - Xtry
      B      =    B - Xtry
      Moved  = .true.
      Extrap = .not. Xinxw
C
C  Decrease The Length Of The Interval Of Uncertainty.
C
      If (Xtry .lt. Zero) Go To 120
      A      = Xw
      Fa     = Fw
      Nsamea = 0
      Go To 300
  120 B      = Xw
      Nsameb = 0
      Braktd = .true.
      Go To 300
C
C  The New Function Value Is No Better Than The Best Point Found So Far.
C  The Point  Xtry  Must Be A New End Point Of The Interval Of
C  Uncertainty.
C
  130 If (Xtry .ge. Zero) Go To 140
      A      = Xtry 
      Fa     = Ftry 
      Nsamea = 0
      Go To 150
  140 B      = Xtry 
      Nsameb = 0
      Braktd = .true.
C
C  The Origin Remains Unchanged But  Xtry  May Qualify As  Xw.
C
  150 If (.not. Wset)   Go To 160
      If ((Ftry - Fw) .gt. Epsaf) Go To 170
      Xv     = Xw
      Fv     = Fw
      Vset   = .true.
  160 Xw     = Xtry 
      Fw     = Ftry 
      Wset   = .true.
      If (Moved) Extrap = Xinxw
      Go To 300
C
C  Ftry  Is No Better Than  Fbest  Or  Fw.  If The Best Point Has Not 
C  Been Moved, There Must Be More Than One Minimum.
C
  170 If (Moved) Go To 175
      Xw     = Xtry 
      Fw     = Ftry 
      Go To 300
C
C  Ftry  Is No Better Than  Fbest  Or  Fw,  But  Xtry  May Become  Xv.
C  Extrap  Has The Value Set In The Previous Entry.
C
  175 If (.not. Vset) Go To 180
      If ((Ftry - Fv) .gt. Epsaf  .and.  Extrap) Go To 300
  180 If (Xinxw) Go To 190
      Xv     = Xtry 
      Fv     = Ftry 
      Vset   = .true.
      Go To 300
  190 If (Vset) Xw = Xv
      If (Vset) Fw = Fv
      Xv     = Xtry 
      Fv     = Ftry 
      Vset   = .true.
C
C  ---------------------------------------------------------------------
C  Check The Termination Criteria.
C  ---------------------------------------------------------------------
  300 Tol    = Tolrel*Alfbst + Tolabs
      Deltaf = Fbest  - Oldf

      Cbest  = Deltaf - Alfbst*Rmu*Oldg 
      If (Alfbst .le. Alfuzz) Sigdec = Deltaf .le. (- Epsaf)
      If (Alfbst .gt. Alfuzz) Sigdec = Cbest  .le.    Epsaf 
      Closef = .false.
      If (Vset) Closef = Abs( Fbest - Fv ) .le. Epsaf
C
      Conv1  = Max( Abs( A ), B )  .le.  (Tol + Tol)
      Conv2 = Moved .and. Sigdec
      Conv3  = Closef  .and.  (Sigdec   .or. 
     *                        (.not. Moved)  .and.  (B .le. Alfuzz))
      Convrg = Conv1  .or.  Conv2  .or.  Conv3
C
      Atrue  = Alfbst + A
      Btrue  = Alfbst + B
      Alfaw  = Alfbst + Xw
      Gap    = B - A
      If (Debug) Write (Nout, 1200) Atrue, Btrue, Gap, Tol, 
     *   Nsamea, Nsameb, Braktd, Closef, Imprvd, Conv1, Conv2, Conv3, 
     *   Extrap, Alfbst, Fbest, Cbest, Alfaw, Fw
      If (Vset) Alfav  = Alfbst + Xv
      If (Debug  .and.  Vset) Write (Nout, 1300) Alfav, Fv
      If (Convrg  .and.  Moved) Go To 910
C
C  Exit If The Step Is Too Small.
C
      If (Btrue   .lt.  Alfsml) Go To 950
      
      If (Nfsrch  .ge.  Mfsrch) Go To 930
      If (.not. Convrg) Go To 400
C
C  A Better Point Has Not Yet Been Found (The Step  Xw  Is No Better
C  Than Step  Zero).  Check That The Change In  F  Is Consistent With A
C  Perturbation In  X  Of  Tol, The Estimate Of The Minimum Spacing
C  Constant.  If The Change In  F  Is Larger Than  Epsaf,  The Value
C  Of  Tol  Is Reduced.
C
      Tol    = Xw/Ten
      Tolabs = Tol
      If (Abs(Fw - Oldf) .gt. Epsaf  .and.  Tol .gt. Toltny) Go To 400
      If (Crampd) Go To 940
      Go To 960
C
C  ---------------------------------------------------------------------
C  Proceed With The Computation Of A Trial Step Length.
C  The Choices Are...
C  1. Parabolic Fit Using Function Values Only.
C  2. Damped Parabolic Fit If The Regular Fit Appears To Be 
C     Consistently Over-Estimating The Distance To The Minimum.
C  3. Bisection, Geometric Bisection, Or A Step Of  Tol  If The
C     Parabolic Fit Is Unsatisfactory.
C  ---------------------------------------------------------------------
  400 Xmidpt = Half*(A + B)
      Q      = Zero 
      S      = Zero 
C
C  ---------------------------------------------------------------------
C  Fit A Parabola.
C  ---------------------------------------------------------------------
C
C  Check If There Are Two Or Three Points For The Parabolic Fit.
C
      Gw = (Fw - Fbest)/Xw
      If (Vset  .and.  Moved) Go To 450 
C
C  Only Two Points Available.  Use  Fbest,  Fw  And The Derivative
C  Oldg.
C
      If (.not. Moved) S = Oldg
      If (      Moved) S = Oldg - Two*Gw
      Q = Two*(Oldg - Gw)
      If (Debug) Write (Nout, 2100)
      Go To 600
C
C  Three Points Available.  Use  Fbest,  Fw  And  Fv.
C
  450 Gv = (Fv - Fbest)/Xv
      S  = Gv - (Xv/Xw)*Gw
      Q  = Two*(Gv - Gw)
      If (Debug) Write (Nout, 2200)
C
C  ---------------------------------------------------------------------
C  Construct An Artificial Interval  (Artifa, Artifb)  In Which The
C  New Estimate Of The Step Length Must Lie.  Set A Default Value Of
C  Xtry  That Will Be Used If The Polynomial Fit Is Rejected.  In The 
C  Following, The Interval  (A,B)  Is Considered The Sum Of Two
C  Intervals Of Lengths  Dtry  And  Daux, With Common End Point At The
C  Best Point (Zero).  Dtry  Is The Length Of The Interval Into Which 
C  The Default  Xtry  Will Be Placed And  Endpnt  Denotes Its Non-Zero
C  End Point.  The Magnitude Of  Xtry  Is Computed So That The Exponents
C  Of  Dtry  And  Daux  Are Approximately Bisected.
C  ---------------------------------------------------------------------
  600 Artifa = A
      Artifb = B
      If (Braktd) Go To 610
C
C  The Minimum Has Not Been Bracketed.  Set An Artificial Upper Bound 
C  By Expanding The Interval  Xw  By A Suitable Factor.
C
      Xtry   = - Factor*Xw
      Artifb =   Xtry
      If (Alfbst + Xtry .lt. Alfmax) Factor = Five*Factor
      Go To 700
C
C  The Minimum Has Been Bracketed.
C  If The Gradient At The Origin Is Being Used For The
C  Polynomial Fit, The Default  Xtry  Is One Tenth Of  Xw.
C
  610 If (Vset  .and.  Moved) Go To 620 
      Xtry   = Xw/Ten
      If (Debug) Write (Nout, 2400) Xtry
      Go To 700
C
C  Three Points Exist In The Interval Of Uncertainty.  Check Whether
C  The Points Are Configured For An Extrapolation Or Interpolation.
C
  620 If (Extrap) Go To 660
C
C  If The Interpolation Appears To Be Consistently Over-Estimating The
C  Distance To The Minimum,  Damp The Interpolation Step.
C
      If (Nsamea .lt. 3  .and.  Nsameb .lt. 3) Go To 630
      Factor = Factor / Five
      S      = Factor * S
      Go To 640
  630 Factor = One
C
C  The Points Are Configured For An Interpolation.  The Artificial
C  Interval Will Be Just  (A,B).   Set  Endpnt  So That  Xtry
C  Lies In The Larger Of The Intervals  (A,0)  And  (0,B).
C
  640 If (Xmidpt .lt. Zero) Endpnt = A
      If (Xmidpt .gt. Zero) Endpnt = B
C
C  If A Bound Has Remained The Same For Three Iterations, Set  Endpnt 
C  So That  Xtry  Is Likely To Replace The Offending Bound. 
C
      If (Nsamea .ge. 3) Endpnt = A
      If (Nsameb .ge. 3) Endpnt = B
      Go To 680
C
C  The Points Are Configured For An Extrapolation.
C
  660 If (Xw .lt. Zero) Endpnt = B
      If (Xw .gt. Zero) Endpnt = A
C
C  Compute The Default Value Of  Xtry.
C
  680 Dtry = Abs( Endpnt )
      Daux = Gap - Dtry
      If (Daux .ge. Dtry)   Xtry = Five*Dtry*(Point1 + Dtry/Daux)/Eleven
      If (Daux .lt. Dtry)   Xtry = Half*Sqrt( Daux )*Sqrt( Dtry )
      If (Endpnt .lt. Zero) Xtry = - Xtry
      If (Debug) Write (Nout, 2500) Xtry, Daux, Dtry
C
C  If The Points Are Configured For An Extrapolation Set The Artificial
C  Bounds So That The Artificial Interval Lies Strictly Within  (A,B).
C  If The Polynomial Fit Is Rejected,  Xtry  Will Remain At The Relevant
C  Artificial Bound.
C
      If (Extrap  .and.  Xtry .le. Zero) Artifa = Xtry
      If (Extrap  .and.  Xtry .gt. Zero) Artifb = Xtry
C
C  ---------------------------------------------------------------------
C  The Polynomial Fits Give  (S/Q)*Xw  As The New Step.
C  Reject This Step If It Lies Outside  (Artifa, Artifb).
C  ---------------------------------------------------------------------
  700 If (Q .eq. Zero) Go To 800
      If (Q .lt. Zero) S = - S
      If (Q .lt. Zero) Q = - Q
      If (S*Xw .lt. Q*Artifa   .or.   S*Xw .gt. Q*Artifb) Go To 800
C
C  Accept The Polynomial Fit. 
C
      Xtry = Zero
      If (Abs( S*Xw ) .ge. Q*Tol) Xtry = (S/Q)*Xw 
      If (Debug) Write (Nout, 2600) Xtry
C
C  ---------------------------------------------------------------------
C  Test For  Xtry  Being Larger Than  Alfmax  Or Too Close To  A  Or  B.
C  ---------------------------------------------------------------------
  800 If (Braktd) Go To 810
C
C  If The Step Is Close To Or Larger Than  Alfmax,  Replace It By
C  Alfmax  (To Force Evaluation Of The Function At The Boundary).
C
      Alfa   = Alfbst + Xtry
      If (Alfmax - Alfa .gt. (Tolrel*Alfmax + Tolabs)) Go To 810
      Braktd = .true.
      Xtry   = Alfmax - Alfbst
      Alfa   = Alfmax
      Go To 900
C
C  Otherwise, The Function Must Not Be Evaluated At A Point Too Close 
C  To  A  Or  B.  (It Has Already Been Evaluated At Both Those Points.)
C
  810 Xmidpt = Half*(A + B)
      If (Xtry .gt. A + Tol  .and.  Xtry .lt. B - Tol) Go To 820
      If (Xmidpt .gt. Zero) Xtry =   Tol
      If (Xmidpt .le. Zero) Xtry = - Tol
C
C
C  F  Must Not Be Calculated Too Close To  Alfbst.
C
  820 If (Abs( Xtry ) .lt. Tol  .and.  Xmidpt .lt. Zero) Xtry = - Tol 
      If (Abs( Xtry ) .lt. Tol  .and.  Xmidpt .ge. Zero) Xtry =   Tol 
      Alfa   = Alfbst + Xtry
C
C  ---------------------------------------------------------------------
C  Exit.
C  ---------------------------------------------------------------------
C
C  New Function Value Required.
C
  900 Inform = 0
      Go To 990
C
C  Convergence Test Satisfied.
C
  910 Inform = 1
      If (Alfa .eq. Alfmax) Inform = 2
      Go To 990
C
C  Mfsrch  Function Evaluations Without Sufficient Decrease, But An
C  Improved Point Was Found.
C
  930 If (.not. Moved) Go To 960
      Inform = 3
      Go To 990
C
C  Zero Step (Alfmax Too Small).
C
  940 Inform = 4
      Go To 990
C
C  Premature Termination.  The Step Is Smaller Than  Alfsml.
C
  950 Inform = 5
      Go To 990
C
C  Zero Step (A Sufficiently Better Point Could Not Be Found).
C
  960 Inform = 6
      Go To 990
C
C  Zero Step (Positive Gradient At The Starting Point).
C
  970 Inform = 7
C
C  Exit.
C
  990 If (Debug) Write (Nout, 3000)
      Return
C
 1000 Format(/ 31h Alfmax  Oldf    Oldg    Tolabs, 1p2e22.14, 1p2e16.8
     *       / 31h Alfuzz  Epsaf   Epsag   Tolrel, 1p2e22.14, 1p2e16.8
     *       / 31h Crampd                        ,  L6)
 1100 Format(/ 31h Alfa    Ftry    Ctry          , 1p2e22.14, 1pe16.8)
 1200 Format(/ 31h A       B       B - A   Tol   , 1p2e22.14, 1p2e16.8
     *       / 31h Nsamea  Nsameb  Braktd  Closef, 2i3, 2l6 
     *       / 31h Imprvd  Convrg  Extrap        ,  L6, 3x, 3l1, L6
     *       / 31h Alfbst  Fbest   Cbest         , 1p2e22.14, 1pe16.8 
     *       / 31h Alfaw   Fw                    , 1p2e22.14)
 1300 Format(  31h Alfav   Fv                    , 1p2e22.14 /)
 2100 Format(30h Parabolic Fit,    Two Points.)
 2200 Format(30h Parabolic Fit,  Three Points.)
 2400 Format(31h Exponent Reduced.  Trial Point, 1p1e22.14) 
 2500 Format(31h Geo. Bisection. Xtry,Daux,Dtry, 1p3e22.14) 
 2600 Format(31h Polynomial Fit Accepted.  Xtry, 1p1e22.14) 
 3000 Format(53h ---------------------------------------------------- /)
C
C  End Of Getptq
      End 


      Subroutine Interp(Ncomp, Nmsh, Xx, Nudim, U, 
     *                    Nmold, Xxold, Uold) 

      Implicit Double Precision (A-H, O-Z)
      Dimension Xx(*), U(Nudim,*), Xxold(*), Uold(Ncomp,*) 

* Blas: Dcopy

      Parameter (Zero = 0.0d+0)

*  Interp Performs Piecewise Linear Interpolation Of The Old
*  Solution Uold At The Nmold Old Mesh Points Xxold Onto The Nmsh
*  New Mesh Points Xx, Producing An Interpolated Solution U.
*  Note That No Assumption Is Made That The New Mesh Has
*  More Points Than The Old, Nor That The New And Old Mesh
*  Points Are Related In A Specific Way (Except That Their First
*  And Last Points Are Identical).

*  By Construction, Xx(1) = Xxold(1).  Copy The First Ncomp
*  Components Of Uold Into Those Of U.

      Call Dcopy(Ncomp, Uold(1,1), 1, U(1,1), 1)  

      I = 2 
      Do 100 Im = 2, Nmsh-1

   50    Continue
         If (I .gt. Nmold) Return

*  Check Whether The Im-Th Point In The New Mesh Lies Strictly
*  To The Right Of, Or To The Left Of (Or Exactly On) The
*  I-Th Point In The Old Mesh.


         If (Xx(Im) .gt. Xxold(I)) Then
            I = I + 1
            Go To 50
         Else
            Xdif = Xxold(I) - Xx(Im)
            If (Xdif .eq. Zero) Then

*  Xx(Im) And Xxold(I) Are Identical.
  
               Call Dcopy(Ncomp, Uold(1,I), 1, U(1,Im), 1)
               I = I + 1
            Else
               Xint = Xxold(I) - Xxold(I-1)
               Xrat = Xdif/Xint
               Do 70 K = 1, Ncomp
                  U(K,Im) = Uold(K,I) + Xrat*(Uold(K,I-1)-Uold(K,I))
   70          Continue
            Endif
         Endif

  100 Continue
      Call Dcopy(Ncomp, Uold(1,Nmold), 1, U(1,Nmsh), 1)
      Return
      End 


      Double Precision Function D1mach(I)
C
C  Double-Precision Machine Constants
C
C  D1mach( 1) = B**(Emin-1), The Smallest Positive Magnitude.
C
C  D1mach( 2) = B**Emax*(1 - B**(-T)), The Largest Magnitude.
C
C  D1mach( 3) = B**(-T), The Smallest Relative Spacing.
C
C  D1mach( 4) = B**(1-T), The Largest Relative Spacing.
C
C  D1mach( 5) = Log10(B)
C
C  To Alter This Function For A Particular Environment,
C  The Desired Set Of Data Statements Should Be Activated By
C  Removing The C From Column 1.
C  On Rare Machines A Static Statement May Need To Be Added.
C  (But Probably More Systems Prohibit It Than Require It.)
C
C  For Ieee-Arithmetic Machines (Binary Standard), One Of The First
C  Two Sets Of Constants Below Should Be Appropriate.  If You Do Not
C  Know Which Set To Use, Try Both And See Which Gives Plausible
C  Values.
C
C  Where Possible, Decimal, Octal Or Hexadecimal Constants Are Used
C  To Specify The Constants Exactly.  Sometimes This Requires Using
C  Equivalent Integer Arrays.  If Your Compiler Uses Half-Word
C  Integers By Default (Sometimes Called Integer*2), You May Need To
C  Change Integer To Integer*4 Or Otherwise Instruct Your Compiler
C  To Use Full-Word Integers In The Next 5 Declarations.
C
C  Comments Just Before The End Statement (Lines Starting With *)
C  Give C Source For D1mach.
C
      Integer Small(2)
      Integer Large(2)
      Integer Right(2)
      Integer Diver(2)
      Integer Log10(2)
      Integer Sc
C
      Double Precision Dmach(5)
C
      Equivalence (Dmach(1),Small(1))
      Equivalence (Dmach(2),Large(1))
      Equivalence (Dmach(3),Right(1))
      Equivalence (Dmach(4),Diver(1))
      Equivalence (Dmach(5),Log10(1))
C
C     Machine Constants For Big-Endian Ieee Arithmetic (Binary Format)
C     Machines In Which The Most Significant Byte Is Stored First,
C     Such As The At&T 3b Series, Motorola 68000 Based Machines (E.g.
C     Sun 3), And Machines That Use Sparc, Hp, Or Ibm Risc Chips.
C
C      Data Small(1),Small(2) /    1048576,          0 /
C      Data Large(1),Large(2) / 2146435071,         -1 /
C      Data Right(1),Right(2) / 1017118720,          0 /
C      Data Diver(1),Diver(2) / 1018167296,          0 /
C      Data Log10(1),Log10(2) / 1070810131, 1352628735 /, Sc/987/
C
C     Machine Constants For Little-Endian (Binary) Ieee Arithmetic
C     Machines In Which The Least Significant Byte Is Stored First,
C     E.g. Ibm Pcs And Other Machines That Use Intel 80x87 Or Dec
C     Alpha Chips.
C
      Data Small(1),Small(2) /          0,    1048576 /
      Data Large(1),Large(2) /         -1, 2146435071 /
      Data Right(1),Right(2) /          0, 1017118720 /
      Data Diver(1),Diver(2) /          0, 1018167296 /
      Data Log10(1),Log10(2) / 1352628735, 1070810131 /, Sc/987/
C
C     Machine Constants For Amdahl Machines.
C
C      Data Small(1),Small(2) /    1048576,          0 /
C      Data Large(1),Large(2) / 2147483647,         -1 /
C      Data Right(1),Right(2) /  856686592,          0 /
C      Data Diver(1),Diver(2) /  873463808,          0 /
C      Data Log10(1),Log10(2) / 1091781651, 1352628735 /, Sc/987/
C
C     Machine Constants For The Burroughs 1700 System.
C
C      Data Small(1) / Zc00800000 /
C      Data Small(2) / Z000000000 /
C
C      Data Large(1) / Zdffffffff /
C      Data Large(2) / Zfffffffff /
C
C      Data Right(1) / Zcc5800000 /
C      Data Right(2) / Z000000000 /
C
C      Data Diver(1) / Zcc6800000 /
C      Data Diver(2) / Z000000000 /
C
C      Data Log10(1) / Zd00e730e7 /
C      Data Log10(2) / Zc77800dc0 /, Sc/987/
C
C     Machine Constants For The Burroughs 5700 System.
C
C      Data Small(1) / O1771000000000000 /
C      Data Small(2) / O0000000000000000 /
C
C      Data Large(1) / O0777777777777777 /
C      Data Large(2) / O0007777777777777 /
C
C      Data Right(1) / O1461000000000000 /
C      Data Right(2) / O0000000000000000 /
C
C      Data Diver(1) / O1451000000000000 /
C      Data Diver(2) / O0000000000000000 /
C
C      Data Log10(1) / O1157163034761674 /
C      Data Log10(2) / O0006677466732724 /, Sc/987/
C
C     Machine Constants For The Burroughs 6700/7700 Systems.
C
C      Data Small(1) / O1771000000000000 /
C      Data Small(2) / O7770000000000000 /
C
C      Data Large(1) / O0777777777777777 /
C      Data Large(2) / O7777777777777777 /
C
C      Data Right(1) / O1461000000000000 /
C      Data Right(2) / O0000000000000000 /
C
C      Data Diver(1) / O1451000000000000 /
C      Data Diver(2) / O0000000000000000 /
C
C      Data Log10(1) / O1157163034761674 /
C      Data Log10(2) / O0006677466732724 /, Sc/987/
C
C     Machine Constants For Ftn4 On The Cdc 6000/7000 Series.
C
C      Data Small(1) / 00564000000000000000b /
C      Data Small(2) / 00000000000000000000b /
C
C      Data Large(1) / 37757777777777777777b /
C      Data Large(2) / 37157777777777777774b /
C
C      Data Right(1) / 15624000000000000000b /
C      Data Right(2) / 00000000000000000000b /
C
C      Data Diver(1) / 15634000000000000000b /
C      Data Diver(2) / 00000000000000000000b /
C
C      Data Log10(1) / 17164642023241175717b /
C      Data Log10(2) / 16367571421742254654b /, Sc/987/
C
C     Machine Constants For Ftn5 On The Cdc 6000/7000 Series.
C
C      Data Small(1) / O"00564000000000000000" /
C      Data Small(2) / O"00000000000000000000" /
C
C      Data Large(1) / O"37757777777777777777" /
C      Data Large(2) / O"37157777777777777774" /
C
C      Data Right(1) / O"15624000000000000000" /
C      Data Right(2) / O"00000000000000000000" /
C
C      Data Diver(1) / O"15634000000000000000" /
C      Data Diver(2) / O"00000000000000000000" /
C
C      Data Log10(1) / O"17164642023241175717" /
C      Data Log10(2) / O"16367571421742254654" /, Sc/987/
C
C     Machine Constants For Convex C-1
C
C      Data Small(1),Small(2) / '00100000'X, '00000000'X /
C      Data Large(1),Large(2) / '7fffffff'X, 'Ffffffff'X /
C      Data Right(1),Right(2) / '3cc00000'X, '00000000'X /
C      Data Diver(1),Diver(2) / '3cd00000'X, '00000000'X /
C      Data Log10(1),Log10(2) / '3ff34413'X, '509f79ff'X /, Sc/987/
C
C     Machine Constants For The Cray 1, Xmp, 2, And 3.
C
C      Data Small(1) / 201354000000000000000b /
C      Data Small(2) / 000000000000000000000b /
C
C      Data Large(1) / 577767777777777777777b /
C      Data Large(2) / 000007777777777777776b /
C
C      Data Right(1) / 376434000000000000000b /
C      Data Right(2) / 000000000000000000000b /
C
C      Data Diver(1) / 376444000000000000000b /
C      Data Diver(2) / 000000000000000000000b /
C
C      Data Log10(1) / 377774642023241175717b /
C      Data Log10(2) / 000007571421742254654b /, Sc/987/
C
C     Machine Constants For The Data General Eclipse S/200
C
C     Small, Large, Right, Diver, Log10 Should Be Declared
C     Integer Small(4), Large(4), Right(4), Diver(4), Log10(4)
C
C     Note - It May Be Appropriate To Include The Following Line -
C     Static Dmach(5)
C
C      Data Small/20k,3*0/,Large/77777k,3*177777k/
C      Data Right/31420k,3*0/,Diver/32020k,3*0/
C      Data Log10/40423k,42023k,50237k,74776k/, Sc/987/
C
C     Machine Constants For The Harris Slash 6 And Slash 7
C
C      Data Small(1),Small(2) / '20000000, '00000201 /
C      Data Large(1),Large(2) / '37777777, '37777577 /
C      Data Right(1),Right(2) / '20000000, '00000333 /
C      Data Diver(1),Diver(2) / '20000000, '00000334 /
C      Data Log10(1),Log10(2) / '23210115, '10237777 /, Sc/987/
C
C     Machine Constants For The Honeywell Dps 8/70 Series.
C
C      Data Small(1),Small(2) / O402400000000, O000000000000 /
C      Data Large(1),Large(2) / O376777777777, O777777777777 /
C      Data Right(1),Right(2) / O604400000000, O000000000000 /
C      Data Diver(1),Diver(2) / O606400000000, O000000000000 /
C      Data Log10(1),Log10(2) / O776464202324, O117571775714 /, Sc/987/
C
C     Machine Constants For The Ibm 360/370 Series,
C     The Xerox Sigma 5/7/9 And The Sel Systems 85/86.
C
C      Data Small(1),Small(2) / Z00100000, Z00000000 /
C      Data Large(1),Large(2) / Z7fffffff, Zffffffff /
C      Data Right(1),Right(2) / Z33100000, Z00000000 /
C      Data Diver(1),Diver(2) / Z34100000, Z00000000 /
C      Data Log10(1),Log10(2) / Z41134413, Z509f79ff /, Sc/987/
C
C     Machine Constants For The Interdata 8/32
C     With The Unix System Fortran 77 Compiler.
C
C     For The Interdata Fortran Vii Compiler Replace
C     The Z'S Specifying Hex Constants With Y'S.
C
C      Data Small(1),Small(2) / Z'00100000', Z'00000000' /
C      Data Large(1),Large(2) / Z'7effffff', Z'Ffffffff' /
C      Data Right(1),Right(2) / Z'33100000', Z'00000000' /
C      Data Diver(1),Diver(2) / Z'34100000', Z'00000000' /
C      Data Log10(1),Log10(2) / Z'41134413', Z'509f79ff' /, Sc/987/
C
C     Machine Constants For The Pdp-10 (Ka Processor).
C
C      Data Small(1),Small(2) / "033400000000, "000000000000 /
C      Data Large(1),Large(2) / "377777777777, "344777777777 /
C      Data Right(1),Right(2) / "113400000000, "000000000000 /
C      Data Diver(1),Diver(2) / "114400000000, "000000000000 /
C      Data Log10(1),Log10(2) / "177464202324, "144117571776 /, Sc/987/
C
C     Machine Constants For The Pdp-10 (Ki Processor).
C
C      Data Small(1),Small(2) / "000400000000, "000000000000 /
C      Data Large(1),Large(2) / "377777777777, "377777777777 /
C      Data Right(1),Right(2) / "103400000000, "000000000000 /
C      Data Diver(1),Diver(2) / "104400000000, "000000000000 /
C      Data Log10(1),Log10(2) / "177464202324, "047674776746 /, Sc/987/
C
C     Machine Constants For Pdp-11 Fortrans Supporting
C     32-Bit Integers (Expressed In Integer And Octal).
C
C      Data Small(1),Small(2) /    8388608,           0 /
C      Data Large(1),Large(2) / 2147483647,          -1 /
C      Data Right(1),Right(2) /  612368384,           0 /
C      Data Diver(1),Diver(2) /  620756992,           0 /
C      Data Log10(1),Log10(2) / 1067065498, -2063872008 /, Sc/987/
C
C      Data Small(1),Small(2) / O00040000000, O00000000000 /
C      Data Large(1),Large(2) / O17777777777, O37777777777 /
C      Data Right(1),Right(2) / O04440000000, O00000000000 /
C      Data Diver(1),Diver(2) / O04500000000, O00000000000 /
C      Data Log10(1),Log10(2) / O07746420232, O20476747770 /, Sc/987/
C
C     Machine Constants For Pdp-11 Fortrans Supporting
C     16-Bit Integers (Expressed In Integer And Octal).
C
C     Small, Large, Right, Diver, Log10 Should Be Declared
C     Integer Small(4), Large(4), Right(4), Diver(4), Log10(4)
C
C      Data Small(1),Small(2) /    128,      0 /
C      Data Small(3),Small(4) /      0,      0 /
C
C      Data Large(1),Large(2) /  32767,     -1 /
C      Data Large(3),Large(4) /     -1,     -1 /
C
C      Data Right(1),Right(2) /   9344,      0 /
C      Data Right(3),Right(4) /      0,      0 /
C
C      Data Diver(1),Diver(2) /   9472,      0 /
C      Data Diver(3),Diver(4) /      0,      0 /
C
C      Data Log10(1),Log10(2) /  16282,   8346 /
C      Data Log10(3),Log10(4) / -31493, -12296 /, Sc/987/
C
C      Data Small(1),Small(2) / O000200, O000000 /
C      Data Small(3),Small(4) / O000000, O000000 /
C
C      Data Large(1),Large(2) / O077777, O177777 /
C      Data Large(3),Large(4) / O177777, O177777 /
C
C      Data Right(1),Right(2) / O022200, O000000 /
C      Data Right(3),Right(4) / O000000, O000000 /
C
C      Data Diver(1),Diver(2) / O022400, O000000 /
C      Data Diver(3),Diver(4) / O000000, O000000 /
C
C      Data Log10(1),Log10(2) / O037632, O020232 /
C      Data Log10(3),Log10(4) / O102373, O147770 /, Sc/987/
C
C     Machine Constants For The Prime 50 Series Systems
C     With 32-Bit Integers And 64v Mode Instructions,
C     Supplied By Igor Bray.
C
C      Data Small(1),Small(2) / :10000000000, :00000100001 /
C      Data Large(1),Large(2) / :17777777777, :37777677775 /
C      Data Right(1),Right(2) / :10000000000, :00000000122 /
C      Data Diver(1),Diver(2) / :10000000000, :00000000123 /
C      Data Log10(1),Log10(2) / :11504046501, :07674600177 /, Sc/987/
C
C     Machine Constants For The Sequent Balance 8000
C
C      Data Small(1),Small(2) / $00000000,  $00100000 /
C      Data Large(1),Large(2) / $ffffffff,  $7fefffff /
C      Data Right(1),Right(2) / $00000000,  $3ca00000 /
C      Data Diver(1),Diver(2) / $00000000,  $3cb00000 /
C      Data Log10(1),Log10(2) / $509f79ff,  $3fd34413 /, Sc/987/
C
C     Machine Constants For The Univac 1100 Series.
C
C      Data Small(1),Small(2) / O000040000000, O000000000000 /
C      Data Large(1),Large(2) / O377777777777, O777777777777 /
C      Data Right(1),Right(2) / O170540000000, O000000000000 /
C      Data Diver(1),Diver(2) / O170640000000, O000000000000 /
C      Data Log10(1),Log10(2) / O177746420232, O411757177572 /, Sc/987/
C
C     Machine Constants For The Vax Unix F77 Compiler
C
C      Data Small(1),Small(2) /        128,           0 /
C      Data Large(1),Large(2) /     -32769,          -1 /
C      Data Right(1),Right(2) /       9344,           0 /
C      Data Diver(1),Diver(2) /       9472,           0 /
C      Data Log10(1),Log10(2) /  546979738,  -805796613 /, Sc/987/
C
C     Machine Constants For The Vax-11 With
C     Fortran Iv-Plus Compiler
C
C      Data Small(1),Small(2) / Z00000080, Z00000000 /
C      Data Large(1),Large(2) / Zffff7fff, Zffffffff /
C      Data Right(1),Right(2) / Z00002480, Z00000000 /
C      Data Diver(1),Diver(2) / Z00002500, Z00000000 /
C      Data Log10(1),Log10(2) / Z209a3f9a, Zcff884fb /, Sc/987/
C
C     Machine Constants For Vax/Vms Version 2.2
C
C      Data Small(1),Small(2) /       '80'X,        '0'X /
C      Data Large(1),Large(2) / 'Ffff7fff'X, 'Ffffffff'X /
C      Data Right(1),Right(2) /     '2480'X,        '0'X /
C      Data Diver(1),Diver(2) /     '2500'X,        '0'X /
C      Data Log10(1),Log10(2) / '209a3f9a'X, 'Cff884fb'X /, Sc/987/
C
C  ***  Issue Stop 779 If All Data Statements Are Commented...
      If (Sc .ne. 987) Stop 779
C  ***  Issue Stop 778 If All Data Statements Are Obviously Wrong...
      If (Dmach(4) .ge. 1.0d0) Stop 778
      If (I .lt. 1  .or.  I .gt. 5) Goto 999
      D1mach = Dmach(I)
      Return
  999 Write(*,1999) I
 1999 Format(' D1mach - I Out Of Bounds',I10)
      Stop
      End
