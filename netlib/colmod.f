C.... *****************************************************************
C.... Beginning Of Sample Problem Description.
C.... *****************************************************************

C.... This Is A Sample Main Program For Colmod, An Automatic 
C.... Continuation Code For Solving `Stiff' Boundary Value Problems.
C.... The Code Is Fully Documented Within The Subroutine Colmod, Which
C.... Follows Directly After This Main Program And The Sample Problem 
C.... Description.

C.... Any Questions, Reports Of Bugs, Or Comments (Good Or Bad), Will
C.... Gladly Be Received By Ross Wright Or Jeff Cash, Department Of
C.... Mathematics, Imperial College, London Sw7 2bz, United Kingdom. 
C.... [ E-Mail r.wright@ic.ac.uk Or j.cash@ic.ac.uk ]

C.... *****************************************************************

      Program Runcol
      Implicit Double Precision(A-H,O-Z)
      Parameter(Ifdim = 500000)
      Parameter(Iidim = 50000)
      Dimension Fspace(Ifdim),Ispace(Iidim)
      Dimension Zeta(2),Tol(2),Ltol(2),Ipar(13),M(1),Fixpnt(1)
      External Fsub,Dfsub,Gsub,Dgsub,Guess
      Common /Prob/ Linnon
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

C.... Note That Both Test Problems Are Solved In Their Second Order
C.... Form And Have Not Been Reduced To A First Order System.

C.... *****************************************************************

C.... Specify The Desired Final Eps Value.

      Epsmin = 1.d-9

C.... Notol Is The Number Of Tolerances

      Notol = 2
      Tol(1) = 1.d-8
      Tol(2) = 1.d-4
      Ltol(1) = 1
      Ltol(2) = 2

C.... Specify The Number Of Components

      Ncomp = 1
      M(1) = 2

C.... The User Chooses A Linear Or Nonlinear Test Problem (0 = Linear,
C.... 1 = Nonlinear).

 10   Write(*,1000) 
      Read*, Linnon

C.... Specify The Boundary Points

      If (Linnon .eq. 1) Then
        Aleft = 0.d0
        Aright = 1.0d0
        Zeta(1) = 0.d0
        Zeta(2) = 1.0d0
      Else If (Linnon .eq. 0) Then
        Aleft = -1.0d0
        Aright = 1.0d0
        Zeta(1) = -1.0d0
        Zeta(2) = 1.0d0
      Else
         Write(*,2000)
         Goto 10
      Endif

C.... Initialise The Variables In The Ipar Array

      Ipar(1) = Linnon
      Ipar(2) = 0
      Ipar(3) = 5
      Ipar(4) = Notol
      Ipar(5) = Ifdim
      Ipar(6) = Iidim
      Ipar(7) = 0
      Ipar(8) = 0
      Ipar(9) = 0
      Ipar(10) = 0
      Ipar(11) = 0
      Ipar(12) = 0
      Ipar(13) = 0

C.... Solve The Problem

      Call Colmod( Ncomp, M, Aleft, Aright, Zeta, Ipar, Ltol,
     +     Tol, Fixpnt, Ispace, Fspace, Iflag, Eps, Epsmin,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess)

C----------------------------------------------------------------------
 1000 Format(//,'Choose A Test Problem (0 = Linear, 1 = Nonlinear).')
 2000 Format(/,'** Only 0 Or 1 Is Allowed.')
C----------------------------------------------------------------------

      End 

      Subroutine Fsub(X,Z,F,Eps)
      Implicit Double Precision(A-H,O-Z)
      Dimension Z(*),F(*)
      Common /Prob/ Linnon
      Common /Valpi/ Pi

      If (Linnon .eq. 0) Then
         Pix = Pi*X
         F(1) = (Z(1)-Eps*Pi*Pi*Cos(Pix)-Cos(Pix))/Eps
      Else
         F(1) = (Z(1)*(1.d0-Z(2)))/Eps
      Endif

      Return
      End

      Subroutine Dfsub(X,Z,Df,Ncomp,Eps)
      Implicit Double Precision(A-H,O-Z)
      Dimension Z(*),Df(Ncomp,*)
      Common /Prob/ Linnon

      If (Linnon .eq. 0) Then
         Df(1,1) = 1.0d0/Eps
         Df(1,2) = 0.0d0
      Else 
         Df(1,1) = (1.d0-Z(2))/Eps
         Df(1,2) = -Z(1)/Eps
      Endif

      Return
      End 

      Subroutine Gsub(I,Z,G,Eps)
      Implicit Double Precision(A-H,O-Z)
      Dimension Z(*)
      Common /Prob/ Linnon

      If (Linnon .eq. 0) Then
         G = Z(1)
      Else 
         If (I .eq. 1) G = Z(1)-1.d0
         If (I .eq. 2) G = Z(1)+1.d0/3.d0
      Endif

      Return 
      End 

      Subroutine Dgsub(I,Z,Dg,Eps)
      Implicit Double Precision(A-H,O-Z)
      Dimension Z(*),Dg(*)

      Dg(1) = 1.d0
      Dg(2) = 0.d0

      Return
      End 

      Subroutine Guess(X,Z,Dmval,Eps)
      Return
      End

C.... *****************************************************************
C.... End Of Sample Problem Description.
C.... *****************************************************************


      Subroutine Colmod(Ncomp, M, Aleft, Aright, Zeta, Ipar, Ltol,
     +     Tol, Fixpnt, Ispace, Fspace, Iflag, Eps, Epsmin,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess)

C.... ******************************************************************

C.... Colmod - Automatic Continuation With Collocation.

C.... This Package Solves `Stiff' Boundary Value Problems For Ordinary 
C.... Differential Equations By Using Continuation, As Described Below.

C.... Colmod Is A Revised Version Of The Package Colnew By Ascher 
C.... And Bader [2], Which In Turn Is A Modification Of The Package 
C.... Colsys By Ascher, Christiansen And Russell [1]. Colmod Has Been
C.... Adapted To Allow An Automatic Continuation Strategy (See [7]) To 
C.... Be Used. The Mesh Selection Algorithm Used In Colmod Differs From 
C.... That Used In Colnew. The Details Of The New Mesh Selection 
C.... Algorithm Are Explained In [6]

C.... ******************************************************************

C...  A LARGE PART OF THE ORIGINAL COLSYS CODE REMAINS UNCHANGED IN
C.... COLMOD, AND WE GRATEFULLY ACKNOWLEDGE THE INPUT OF THE AUTHORS OF
C.... COLSYS AND COLNEW.

C.... ******************************************************************

C----------------------------------------------------------------------
C                           P A R T  1
C   Main Storage Allocation And Program Control Subroutines -- See [7] 
C----------------------------------------------------------------------

C.... ******************************************************************

C.... Adapted From Colsys And Colnew By:

C....          R. W. Wright
C....          Department Of Mathematics,
C....          Imperial College Of Science, Technology And Medicine,
C....          South Kensington,
C....          London Sw7 2bz,
C....          England.                    

C....          J. R. Cash
C....          Department Of Mathematics,
C....          Imperial College Of Science, Technology And Medicine,
C....          South Kensington,
C....          London Sw7 2bz,
C....          England.                    

C....          E-Mail To r.wright@ic.ac.uk Or j.cash@ic.ac.uk

C.... ******************************************************************

C.... Purpose

C.... This Package Solves A Multi-Point Boundary Value Problem For 
C.... A Mixed Order System Of Odes Involving A Small Positive 
C...  Parameter (Which We Denote By Eps). Such A System Is Given By

C....         (M(I))
C....        U       =  F  ( X; Z(U(X)), Eps )     I = 1, ... ,Ncomp
C....         I          I

C....                                        Aleft .lt. X .lt. Aright,


C....        G  ( Zeta(J); Z(U(Zeta(J))), Eps ) = 0,  J = 1, ... ,Mstar,
C....         J
C....                                    Mstar = M(1)+M(2)+...+M(Ncomp),


C.... Where                           T
C....        U = (U , U , ... ,U     )  Is The Exact Solution Vector
C....              1   2        Ncomp

C....         (Mi)
C....        U     Is The Mi (= M(I)) Th Derivative Of U  W.r.t. X
C....         I                                         I

C....                           (1)         (M1-1)   
C....        Z(U(X)) = ( U (X),U (X), ... ,U (X), ... , U (X), ...
C....                     1     1           1            2

C....                                          (1)         (Mncomp-1)
C....                            ... , U (X), U (X), ... ,U (X)       )
C....                                   Ncomp  Ncomp       Ncomp


C....        F (X,Z(U),Eps)  Is A (Generally) Nonlinear Function Of
C....         I
C....                            Z(U) = Z(U(X)).

C....        G (Zeta(J);z(U), Eps)  Is A (Generally) Nonlinear Function
C....         J
C....                          Used To Represent A Boundary Condition.

C.... The Boundary Points Satisfy
C....          Aleft .le. Zeta(1) .le. .. .le. Zeta(Mstar) .le. Aright

C.... The Orders, M(I), Of The Differential Equations Satisfy
C....                   1 .le. M(I) .le. 4,    For All 1 <= I <= Ncomp.

C.... ******************************************************************

C.... Method

C.... The Method Used To Approximate The Solution U Is
C.... Collocation At Gaussian Points, With The Requirement That There 
C.... Are M(I)-1 Continuous Derivatives In The I-Th Component, 
C.... I = 1, ..., Ncomp. Here, K Is The Number Of Collocation Points 
C.... (Stages) Per Subinterval And Is Chosen Such That K .ge. Max M(I).
C.... A Runge-Kutta-Monomial Basis Representation Is Utilized.

C.... The Continuation Algorithm Used Is Based Upon Predicting The 
C.... Behaviour Of The Monitor Function (Used For Mesh Selection). The
C.... Algorithm Is Fully Described In [7]. The Mesh Selection Algorithm
C.... Used In Our Continuation Framework Is Described In [6]. 

C.... ******************************************************************

C.... Restriction On The Choice Of Eps

C.... The Type Of Problems Which Colmod Is Designed To Solve Typically 
C.... Involve A Small Positive Parameter 0 < Eps << 1. As Eps Becomes
C.... Progressively Smaller, The Problem Normally Becomes Increasingly 
C.... Difficult To Approximate Numerically (For Example, Due To The 
C.... Appearance Of Narrow Transition Layers In The Profile Of The 
C.... Analytic Solution). The Idea Of Continuation Is To Solve A Chain 
C.... Of Problems In Which The Parameter Eps Decreases Monotonically
C.... Towards Some Desired Value. That Is, We Attempt To Solve A 
C.... Sequence Of Problems

C....    Eps_s  > Eps_1  > Eps_2  > Eps_3  >  .....  > Eps_min  > 0

C.... Where Eps_s Is A User Provided Starting Value And Eps_min Is A
C.... User Desired Final Value For The Parameter. Given Eps_s And 
C.... Eps_min, The Code Automatically Selects The Intermediate Eps
C.... Values And Passes On Meshes And Solutions (If Solving A Nonlinear
C.... Problem) At Each Step. We Note That The Intermediate Problems Are
C.... Not Necessarily Solved To The User Requested Accuracy. However,
C.... It Is Expected That The Intermediate Problems Will, In The Worst
C.... Cases, Fail The Requested Accuracy By Only A Relatively Small
C.... Amount. For A Detailed Discussion Of This, See [8].

C.... The Dependence On A Small Positive Parameter Is Not Unduly
C.... Restrictive, Since Problems Involving A Large Parameter
C.... Can Simply Be Re-Expressed As Problems Involving The Reciprocal
C.... Of A Small Number. Our Experience Is That It Is Possible To 
C.... Identify Such A Parameter In Most Stiff Problems Of Practical
C.... Interest.

C.... ******************************************************************

C.... References

C....    [1] U. Ascher, J. Christiansen And R. D. Russell,
C....        Collocation Software For Boundary-Value Odes,
C....        Acm Trans. Math Software 7 (1981), 209-222.
C....        This Paper Contains Examples Where Use Of The Code
C....        Is Demonstrated.

C....    [2] G. Bader And U. Ascher,
C....        A New Basis Implementation For A Mixed Order
C....        Boundary Value Ode Solver,
C....        Siam J. Scient. Stat. Comput. (1987).

C....    [3] U. Ascher, J. Christiansen And R. D. Russell,
C....        A Collocation Solver For Mixed Order
C....        Systems Of Boundary Value Problems,
C....        Math. Comp. 33 (1979), 659-679.

C....    [4] U. Ascher, J. Christiansen And R. D. Russell,
C....        Colsys - A Collocation Code For Boundary
C....        Value Problems,
C....        Lecture Notes Comp.sc. 76, Springer Verlag,
C....        B. Childs Et. Al. (Eds.) (1979), 164-185.

C....    [5] C. Deboor And R. Weiss,
C....        Solveblok: A Package For Solving Almost Block Diagonal
C....        Linear Systems,
C....        Acm Trans. Math. Software 6 (1980), 80-87.

C....    [6] R. Wright, J. Cash And G. Moore, 
C....        Mesh Selection For Stiff Two-Point Boundary Value 
C....        Problems,
C....        Numerical Algorithms 7 (1994), 205-224.

C....    [7] J. R. Cash, G. Moore And R. W. Wright, 
C....        An Automatic Continuation Strategy For The Solution Of
C....        Singularly Perturbed Linear Two-Point Boundary Value
C....        Problems,
C....        J. Comp. Phys. 122 (1995), 266-279.

C....    [8] J. R. Cash, G. Moore And R. W. Wright, 
C....        An Automatic Continuation Strategy For The Solution Of
C....        Singularly Perturbed Nonlinear Two-Point Boundary Value
C....        Problems,
C....        To Appear.

C.... ******************************************************************

C.... Input To Colmod

C....   Variables

C....   Ncomp - No. Of Differential Equations            (Ncomp .le. 20)

C....   M(J) - Order Of The J-Th Differential Equation     (M(J) .le. 4)
C....             ( Mstar = M(1) + ... + M(Ncomp) .le. 40 ) 

C....   Aleft - Left End Of Interval

C....   Aright - Right End Of Interval

C....   Zeta(J) - J-Th Side Condition Point (Boundary Point). Must
C....             Have  Zeta(J) .le. Zeta(J+1). All Side Condition
C....             Points Must Be Mesh Points In All Meshes Used,
C....             See Description Of Ipar(10) And Fixpnt Below.

C....   Ipar - An Integer Array Dimensioned At Least 13.
C....          A List Of The Parameters In Ipar And Their Meaning 
C....          Follows. Some Parameters Are Renamed In Colmod; These
C....          New Names Are Given In Parentheses.

C....   Ipar(1)                                             ( = Nonlin )
C....           = 0 If The Problem Is Linear
C....           = 1 If The Problem Is Nonlinear

C....   Ipar(2) = No. Of Collocation Points Per Subinterval      ( = K )
C....             Where Max M(I) .le.  K .le. 7 . If Ipar(2) = 0 Then
C....             Colmod Sets  K = Max ( Max M(I)+1, 5-Max M(I) )

C....   Ipar(3) = No. Of Subintervals In The Initial Mesh.       ( = N )
C....             If Ipar(3) = 0 Then Colmod Arbitrarily Sets N = 5.

C....   Ipar(4) = No. Of Solution And Derivative Tolerances.  ( = Ntol )
C....             We Require  0 .lt. Ntol .le. Mstar.

C....   Ipar(5) = Dimension Of Fspace.                       ( = Ndimf )

C....   Ipar(6) = Dimension Of Ispace.                       ( = Ndimi )

C....   Ipar(7) -  Output Control                           ( = Iprint )
C....           = -1 For Full Diagnostic Printout
C....           = 0 For Selected Printout
C....           = 1 For No Printout

C....   Ipar(8)                                              ( = Iread )
C....           = 0 Causes Colmod To Generate A Uniform Initial Mesh.
C....           = 1 If The Initial Mesh Is Provided By The User.  It
C....               Is Defined In Fspace As Follows:  The Mesh
C....               Aleft = X(1) .lt. X(2) .lt.  ...  .lt. X(N) 
C....                                            .lt. X(N+1) = Aright
C....               Will Occupy  Fspace(1), ..., Fspace(N+1). The
C....               User Needs To Supply Only The Interior Mesh
C....               Points  Fspace(J) = X(J), J = 2, ..., N.

C....   Ipar(9)                                             ( = Iguess )
C....           = 0 If No Initial Guess For The Solution Is
C....               Provided.
C....           = 1 If An Initial Guess Is Provided By The User
C....               In Subroutine  Guess.
C....           = 2 If An Initial Mesh And Approximate Solution
C....               Coefficients Are Provided By The User In  Fspace.
C....               (The Former And New Mesh Are The Same).
C....           = 3 If A Former Mesh And Approximate Solution
C....               Coefficients Are Provided By The User In Fspace,
C....               And The New Mesh Is To Be Taken Twice As Coarse;
C....               I.e.,Every Second Point From The Former Mesh.
C....           = 4 If In Addition To A Former Initial Mesh And
C....               Approximate Solution Coefficients, A New Mesh
C....               Is Provided In Fspace As Well.

C....               The Discussion Directly Below Describes How The
C....               Variable Ipar(9) Might Be Used. However, As We
C....               Are Already Working In A Continuation Framework, It 
C....               Is Unlikely That The User Will Want To Use Anything
C....               Other Than Ipar(9) = 0. An Exception To This Might
C....               Be The Case Where The User Wants To Solve The Same
C....               Problem With A Different Value Of Eps_min. In This 
C....               Case, The Mesh And Solution For One Value Of Eps_min
C....               Might Be Provided As An Initial Mesh And Initial  
C....               Guess For The Next (Smaller) Value Of Eps_min (By
C....               Setting Ipar(9) = 2 Or 3).

C....   [ A Formerly Obtained Solution Can Be Used As The
C....     First Approximation For The Nonlinear Iteration For A
C....     New Problem By Setting Ipar(9) = 2, 3 Or 4.

C....     If The Former Solution Has Just Been Obtained Then The
C....     Values Needed To Define The First Approximation Are
C....     Already In Ispace And Fspace.
C....     Alternatively, If The Former Solution Was Obtained In A
C....     Previous Run, And Its Coefficients Were Saved, Then Those
C....     Coefficients Must Be Put Back Into
C....     Ispace(1),..., Ispace(7+Ncomp)    And
C....     Fspace(1),..., Fspace(Ispace(7)).

C....     For Ipar(9) = 2 Or 3 Set Ipar(3) = Ispace(1) ( = The
C....     Size Of The Previous Mesh ).

C....     For Ipar(9) = 4 The User Specifies A New Mesh Of N 
C....     Subintervals As Follows.
C....     The Values In  Fspace(1),...,Fspace(Ispace(7)) Have To Be
C....     Moved By N+1 Locations To Fspace(N+2),..,Fspace(Ispace(7)+N+1)
C....     And The New Mesh Is Specified In Fspace(1),..., Fspace(N+1).
C....     Also Set Ipar(3) = N. ]

C....   Ipar(10) =  No. Of Fixed Points In The Mesh Other Than Aleft
C....               And Aright. ( = Nfxpnt , The Dimension Of Fixpnt)
C....               The Code Requires That All Side Condition Points
C....               Other Than Aleft And Aright (See Description Of
C....               Zeta ) Be Included As Fixed Points In Fixpnt.

C....   Ipar(11) -                                         ( = Maxcon )
C....               The Maximum Number Of Continuation Steps Taken When
C....               Attempting To Solve A Particular Problem. The 
C....               Default Value Is Maxcon = 50 (When Ipar(11) = 0). 

C....   Ipar(12) -                                         ( = Itsaim ) 
C....               Used Only For Nonlinear Problems. The Continuation
C....               Steps Are Chosen With The Aim That The Number Of 
C....               Newton Iterations Required For Convergence On The
C....               First Mesh, For Each Continuation Problem, Is Less
C....               Than Or Equal To Itsaim (See [8] For More Details).
C....               The Default Value Is Itsaim = 7 (If The User 
C....               Specifies Ipar(12) = 0). The User Should Note That 
C....               Small Values Of Itsaim (E.g. Itsaim = 1,2) May 
C....               Restrict The Size Of The Continuation Steps 
C....               Significantly.

C....   Ipar(13) -  Indicates Whether The Initial Value Of Eps Is 
C....               Provided By The User. If Ipar(13) .eq. 1 The User 
C....               Must Specify The Initial Eps Value Before Calling 
C....               Subroutine Colmod. If Ipar(13) .eq. 0, Eps Takes 
C....               The Default Starting Value Of Eps = 0.5.

C....   Ltol - An Array Of Dimension Ipar(4). Ltol(J) = L  Specifies
C....          That The J-Th Tolerance In  Tol  Controls The Error
C....          In The L-Th Component Of Z(U). We Also Require That
C....          1 .le. Ltol(1) .lt. Ltol(2) .lt.  ...  
C....                                     .lt. Ltol(Ntol) .le. Mstar

C....   Tol - An Array Of Dimension Ipar(4). Tol(J) Is The Mixed 
C....         Relative/Absolute Error Tolerance On The Ltol(J)-Th
C....         Component Of Z(U). Thus, The Code Attempts To Satisfy
C....         For J = 1,...,Ntol  On Each Subinterval
C....         Abs(Z(V)-Z(U))       .le. Tol(J)*Max(1,Abs(Z(U))       )
C....                       Ltol(J)                           Ltol(J)
C....
C....         If V(X) Is The Approximate Solution Vector.

C....   Fixpnt - An Array Of Dimension Ipar(10).   It Contains
C....            The Points, Other Than Aleft And Aright, Which
C....            Are To Be Included In Every Mesh.

C....   Ispace - An Integer Work Array Of Dimension Ipar(6).
C....            Its Size Provides A Constraint On Nmax,
C....            The Maximum Number Of Subintervals. Choose
C....            Ipar(6) According To The Formula
C....                  Ipar(6)  .ge.  Nmax*Nsizei
C....            Where
C....                  Nsizei = 3 + Kdm
C....            With
C....                  Kdm = Kd + Mstar  ;  Kd = K * Ncomp ;
C....                  Nrec = No. Of Right End Boundary Conditions.

C....   Fspace - A Real Work Array Of Dimension Ipar(5).
C....            Its Size Provides A Constraint On Nmax.
C....            Choose Ipar(5) According To The Formula
C....                  Ipar(5)  .ge.  Nmax*Nsizef
C....            Where, If We Are Solving A Linear Problem,
C....                   Nsizef = 7 + 3 * Mstar + (5+Kd) * Kdm +
C....                                   (2*Mstar-Nrec) * 2*Mstar,
C....            And If We Are Solving A Nonlinear Problem
C....                   Nsizef = 8 + 4 * Mstar + (5+Kd) * Kdm +
C....                                   (2*Mstar-Nrec) * 2*Mstar + Kd.

C....   Iflag - The Mode Of Return From Colmod.
C....         =  1  For Normal Return
C....         =  0  If The Collocation Matrix Is Singular For The Final
C....               Continuation Problem.
C....         = -1  If The Expected No. Of Subintervals Exceeds Storage
C....               Specifications.
C....         = -2  If The Nonlinear Iteration Has Not Converged For The
C....               Final Continuation Problem.
C....         = -3  If There Is An Input Data Error.

C....   Eps - The Initial Value Of The Continuation Parameter. Eps Must
C....         Be Set Before Calling Colmod If Ipar(13) .eq. 1. (If
C....         Ipar(13) .eq. 0 Eps Takes The Default Starting Value Of 
C....         Eps = 0.5.) For Many Singular Perturbation Type Problems, 
C....         The Choice Of 0.1 < Eps < 1 Represents A (Fairly) Easy 
C....         Problem. The User Should Attempt To Specify An Initial 
C....         Problem That Is Not `Too' Challenging. The Code Assumes 
C....         That Eps = 1 Represents An `Easy' Problem. On Return From 
C....         Colmod The Variable Eps Contains The Value Used When 
C....         Solving The Final Continuation Problem (Usually Epsmin).
C....         If Ipar(13) .eq. 1, Eps Must Be Initialised Strictly Less
C....         Than 1 And Greater Than 0.

C....   Epsmin - The Desired Value Of Eps For Which The User Would 
C....            Like To Solve The Problem. Epsmin Must Be Less Than
C....            Or Equal To Eps.

C.... ******************************************************************

C.... User Supplied Subroutines

C.... The Following Subroutines Must Be Declared External In The
C.... Main Program Which Calls Colmod.

C.... Fsub  - Name Of Subroutine For Evaluating F(X,Z(U(X)),Eps)  = 
C....                        T
C....         (F ,...,F     )  At A Point X In (Aleft,Aright).  It
C....           1      Ncomp
C....         Should Have The Heading

C....                   Subroutine Fsub (X , Z , F, Eps)

C....         Where F Is The Vector Containing The Value Of Fi(X,Z(U))
C....         In The I-Th Component And                              T
C....                                   Z(U(X)) = (Z(1),...,Z(Mstar))
C....         Is Defined As Above Under Purpose .

C.... Dfsub - Name Of Subroutine For Evaluating The Jacobian Of
C....         F(X,Z(U),Eps) At A Point X.  It Should Have The Heading

C....                   Subroutine Dfsub (X , Z , Df, Ncomp, Eps)

C....         Where Z(U(X)) Is Defined As For Fsub And The (Ncomp) By
C....         (Mstar) Array Df Should Be Filled By The Partial Deriv-
C....         Atives Of F, Viz, For A Particular Call One Calculates
C....                            Df(I,J) = Dfi / Dzj, I = 1,...,Ncomp
C....                                                 J = 1,...,Mstar.

C.... Gsub  - Name Of Subroutine For Evaluating The I-Th Component Of
C....         G(X,Z(U(X)),Eps) = G (Zeta(I),Z(U(Zeta(I))),Eps) At Point
C....                             I
C....         X = Zeta(I) Where 1 .le. I .le. Mstar. It Should Have The 
C....         Heading

C....                   Subroutine Gsub (I , Z , G, Eps)

C....         Where Z(U) Is As For Fsub, And I And G = G  Are As Above.
C....                                                   I
C....         Note That In Contrast To F In  Fsub , Here Only One 
C....         Value Per Call Is Returned In G.

C....    Dgsub - Name Of Subroutine For Evaluating The I-Th Row Of
C....            The Jacobian Of G(X,U(X),Eps). It Has The Heading

C....                   Subroutine Dgsub (I , Z , Dg, Eps)

C....            Where Z(U) Is As For Fsub, I As For Gsub And The Mstar-
C....            Vector Dg Should Be Filled With The Partial Derivatives
C....            Of G, Viz, For A Particular Call One Calculates
C....            Dg(I,J) = Dgi / Dzj      J = 1,...,Mstar.

C....    Guess - Name Of Subroutine To Evaluate The Initial 
C....            Approximation For  Z(U(X)) And For Dmval(U(X)) = Vector
C....            Of The Mj-Th Derivatives Of U(X). It Should Have The
C....            Heading

C....                   Subroutine Guess (X , Z , Dmval, Eps)

C....            Note That This Subroutine Is Needed Only If Using
C....            Ipar(9) = 1, And Then All  Mstar  Components Of Z
C....            And  Ncomp  Components Of  Dmval  Should Be Specified
C....            For Any X,  Aleft .le. X .le. Aright .

C.... ******************************************************************

C.... Use Of Output From Colmod - Solution Evaluation   

C.... On Return From Colmod, The Arrays Fspace And Ispace
C.... Contain Information Specifying The Approximate Solution.
C.... The User Can Produce The Solution Vector  Z( U(X) )  At
C.... Any Point X, Aleft .le. X .le. Aright, By The Statement,

C....          Call Appsln (X, Z, Fspace, Ispace)

C.... When Saving The Coefficients For Later Reference, Only
C.... Ispace(1),...,Ispace(7+Ncomp)    And
C.... Fspace(1),...,Fspace(Ispace(7))    Need To Be Saved As
C.... These Are The Quantities Used By Appsln.

C.... ******************************************************************

C.... Package Subroutines

C.... The Following Description Gives A Brief Overview Of How The
C.... Procedure Is Broken Down Into The Subroutines Which Make Up
C.... The Package Called Colmod . For Further Details The
C.... User Should Refer To Documentation In The Various Subroutines
C.... And To The References Cited Above.

C.... The Subroutines Fall Into Four Groups:

C Part 1 - The Main Storage Allocation And Program Control Subroutines

C.... Colmod -  Tests Input Values, Does Initialization And Breaks Up
C....          The Work Areas, Fspace And Ispace, Into The Arrays
C....          Used By The Program. After Initialisation, This
C....          Subroutine Contains The Actual Continuation Algorithm.
C....          The Continuation Step Sizes Are Selected In Colmod.

C.... Contrl - Is The Actual Driver Of The Package. This Routine
C....          Contains The Strategy For Nonlinear Equation Solving.

C.... Skale  - Provides Scaling For The Control
C....          Of Convergence In The Nonlinear Iteration.

C Part 2 - Mesh Selection And Error Estimation Subroutines

C.... Consts - Is Called Once By Colmod  To Initialize Constants
C....          Used For Error Estimation And Mesh Selection.

C.... Newmsh - Generates Meshes. 

C.... Errchk - Produces Error Estimates And Checks Against The
C....          Tolerances At Each Subinterval

C Part 3 - Collocation System Set-Up Subroutines

C.... Lsyslv - Controls The Set-Up And Solution Of The Linear
C....          Algebraic Systems Of Collocation Equations Which
C....          Arise At Each Newton Iteration.

C.... Gderiv - Is Used By Lsyslv To Set Up The Equation Associated
C....          With A Side Condition Point.

C.... Vwblok - Is Used By Lsyslv To Set Up The Equation(S) Associated
C....          With A Collocation Point.

C.... Gblock - Is Used By Lsyslv To Construct A Block Of The Global
C....          Collocation Matrix Or The Corresponding Right Hand
C....          Side.

C Part 4 - Service Subroutines

C.... Appsln - Sets Up A Standard Call To  Approx .

C.... Approx - Evaluates A Piecewise Polynomial Solution.
      
C.... Rkbas  - Evaluates The Mesh Independent Runge-Kutta Basis

C.... Vmonde - Solves A Vandermonde System For Given Right Hand
C....          Side

C.... Horder - Evaluates The Highest Order Derivatives Of The
C....          Current Collocation Solution Used For Mesh Refinement.

C Part 5 - Linear Algebra  Subroutines

C.... To Solve The Global Linear Systems Of Collocation Equations
C.... Constructed In Part 3, Colmod Uses A Column Oriented Version
C.... Of The Package  Solveblok Originally Due To De Boor And Weiss.
C
C.... To Solve The Linear Systems For Static Parameter Condensation
C.... In Each Block Of The Collocation Equations, The Linpack
C.... Routines  Dgefa And  Dgesl  Are Included. But These
C.... May Be Replaced When Solving Problems On Vector Processors
C.... Or When Solving Large Scale Sparse Jacobian Problems.
C----------------------------------------------------------------------

      Implicit Double Precision (A-H,O-Z)
      Dimension M(*), Zeta(*), Ipar(*), Ltol(*), Tol(*), Dummy(1),
     +          Fixpnt(*), Ispace(*), Fspace(*)
      Dimension Phi(3), E(3), Pmax(2), Hord(2)

      Parameter ( Zero = 0.0d+0, One = 1.0d+0, Two = 2.0d+0 )
      Parameter ( Three = 3.0d+0, Third = 1.0d+0/3.0d+0, Huge = 1.d+30 )

      Common /Flags/ Ifinal,Iatt,Iback,Iprec
      Common /Convg/ Nits
      Common /Mshvar/ Pmax,Hord,Hsml
      Common /Colout/ Precis, Iout, Iprint
      Common /Colloc/ Rho(7), Coef(49)
      Common /Colord/ K, Nc, Mstar, Kd, Mmax, Mt(20)
      Common /Colapr/ N, Nold, Nmax, Nz, Ndmz, Mshflg
      Common /Colsid/ Tzeta(40), Tleft, Tright, Izeta, Idum
      Common /Colnln/ Nonlin, Iter, Limit, Iguess
      Common /Colest/ Wgtmsh(40), Wgterr(40), Tolin(40),
     +                Root(40), Jtol(40), Lttol(40), Ntol

      External Fsub, Dfsub, Gsub, Dgsub, Guess

C.... *****************************************************************

C.... The Subroutine Colmod Serves As An Interface With
C.... The Package Of Subroutines Referred To Collectively As
C.... Colmod. The Subroutine Serves To Test Some Of The Input
C.... Parameters, Rename Some Of The Parameters (To Make Under-
C.... Standing Of The Coding Easier), To Do Some Initialization,
C.... And To Break The Work Areas Fspace And Ispace Up Into The
C.... Arrays Needed By The Program.

C.... The Second Part Of This Subroutine Contains The Actual Algorithm
C.... For Selecting The Continuation Parameter.

C.... ******************************************************************

C.... Specify Machine Dependent Output Unit  Iout  And Compute Machine
C.... Dependent Constant  Precis = 100 * Machine Unit Roundoff

      Iout = 6
      Precis = One
   10 Precis = Precis / Two
      Precp1 = Precis + One
      If ( Precp1 .gt. One )                       Go To 10
      Precis = Precis * 100.d0
C
C.... In Case Incorrect Input Data Is Detected, The Program Returns
C.... Immediately With Iflag = -3.
C
      Iflag = -3
      If ( Ncomp .lt. 1 .or. Ncomp .gt. 20 )        Return
      Do 20 I = 1,Ncomp
         If ( M(I) .lt. 1 .or. M(I) .gt. 4 )        Return
   20 Continue
C
C.... Rename Some Of The Parameters And Set Default Values.
C
      Nonlin = Ipar(1)
      K = Ipar(2)
      N = Ipar(3)
      If ( N .eq. 0 )  N = 5
      Iread = Ipar(8)
      Iguess = Ipar(9)
      If ( Nonlin .eq. 0 .and. Iguess .eq. 1 )  Iguess = 0
      If ( Iguess .ge. 2 .and. Iread .eq. 0 )   Iread = 1
      Ntol = Ipar(4)
      Ndimf = Ipar(5)
      Ndimi = Ipar(6)
      Nfxpnt = Ipar(10)
      Iprint = Ipar(7)
      Maxcon = Ipar(11)
      If ( Maxcon .eq. 0 ) Maxcon = 50
      Itsaim = Ipar(12)
      If ( Itsaim .eq. 0 ) Itsaim = 7
      If ( Ipar(13) .eq. 0 ) Eps = 0.5d0
      Mstar = 0
      Mmax = 0
      Do  30 I = 1, Ncomp
         Mmax = Max ( Mmax, M(I) )
         Mstar = Mstar + M(I)
         Mt(I) = M(I)
   30 Continue
      If ( K .eq. 0 )   K = Max( Mmax + 1 , 5 - Mmax )
      Do 40 I = 1, Mstar
   40 Tzeta(I) = Zeta(I)
      Do 50 I = 1, Ntol
         Lttol(I) = Ltol(I)
   50 Tolin(I) = Tol(I)
      Tleft = Aleft
      Tright = Aright
      Nc = Ncomp
      Kd = K * Ncomp

C.... Print The Input Data For Checking.

      If ( Iprint .gt. 0 )                         Go To 80
      If ( Nonlin .gt. 0 )                          Go To 60
      Write (Iout,1000) Ncomp, (M(Ip), Ip = 1,Ncomp)
      Go To 70
   60 Write (Iout,1001) Ncomp, (M(Ip), Ip = 1,Ncomp)
   70 Write (Iout,1002) (Zeta(Ip), Ip = 1,Mstar)
      If ( Nfxpnt .gt. 0 )
     +   Write (Iout,1003) Nfxpnt, (Fixpnt(Ip), Ip = 1,Nfxpnt)
      Write (Iout,1004) K
      Write (Iout,1005) (Ltol(Ip), Ip = 1,Ntol)
      Write (Iout,1006) (Tol(Ip), Ip = 1,Ntol)
      If (Iguess .ge. 2) Write (Iout,1007)
      Write (Iout,1008) Eps, Epsmin
   80 Continue

C.... Check For Correctness Of Data

      If ( K .lt. 0 .or. K .gt. 7 )                 Return
      If ( N .lt. 0 )                               Return
      If ( Maxcon .lt. 0 )                          Return
      If ( Iread .lt. 0 .or. Iread .gt. 2 )         Return
      If ( Iguess .lt. 0 .or. Iguess .gt. 4 )       Return
      If ( Ntol .lt. 0 .or. Ntol .gt. Mstar )       Return
      If ( Nfxpnt .lt. 0 )                          Return
      If ( Iprint .lt. (-1) .or. Iprint .gt. 1 )    Return
      If ( Mstar .lt. 0 .or. Mstar .gt. 40 )        Return
      If ( Eps .le. Zero .or. Eps .ge. One )        Return
      If ( Eps .lt. Epsmin .or. Epsmin .le. Zero)   Return
      If ( Nonlin .eq. 1 .and. Itsaim .lt. 0 )      Return

      Ip = 1
      Do 100 I = 1, Mstar
      If ( Abs(Zeta(I) - Aleft) .lt. Precis  .or. 
     +     Abs(Zeta(I) - Aright) .lt. Precis )     Go To 100
   90 If ( Ip .gt. Nfxpnt )                         Return
        If ( Zeta(I) - Precis .lt. Fixpnt(Ip) )     Go To 95
        Ip = Ip + 1
      Go To 90
   95 If ( Zeta(I) + Precis .lt. Fixpnt(Ip) )       Return
  100 Continue

C.... Set Limits On Iterations And Initialize Counters.
C.... Limit = Maximum Number Of Newton Iterations Per Mesh.
C.... See Subroutine Newmsh For The Role Of Mshflg.

      Mshflg = 0
      Limit = 40

C.... Compute The Maxium Possible N For The Given Sizes Of
C.... Ispace  And  Fspace.

      Nrec = 0
      Do 110 I = 1, Mstar
           Ib = Mstar + 1 - I
           If ( Zeta(Ib) .ge. Aright )  Nrec = I
  110 Continue
      Nfixi = Mstar 
      Nsizei = 3 + Kd + Mstar
      If (Nonlin .eq. 1) Then
         Nfixf = Nrec * (2*Mstar) + 6 * Mstar + 6
         Nsizef = 8 + 4 * Mstar + (Kd+5) * (Kd+Mstar) + 
     +        (2*Mstar-Nrec) * 2*Mstar + Kd
      Else
         Nfixf = Nrec * (2*Mstar) + 5 * Mstar + 5
         Nsizef = 7 + 3 * Mstar + (Kd+5) * (Kd+Mstar) + 
     +        (2*Mstar-Nrec) * 2*Mstar
      Endif
      Nmaxf = (Ndimf - Nfixf) / Nsizef
      Nmaxi = (Ndimi - Nfixi) / Nsizei
      If ( Iprint .lt. 0 )  Write(Iout,1009) Nmaxf, Nmaxi
      Nmax = Min( Nmaxf, Nmaxi )
      If ( Nmax .lt. N )                            Return
      If ( Nmax .lt. Nfxpnt+1 )                     Return
      If (Nmax .lt. 2*Nfxpnt+2 .and. Iprint .lt. 1)  Write(Iout,1010)

C.... Generate Pointers To Break Up Fspace And Ispace.

      Lxi    = 1
      Lg     = Lxi + Nmax + 1
      Lxiold = Lg + 2*Mstar * (Nmax * (2*Mstar-Nrec) + Nrec)
      Lw     = Lxiold + Nmax + 1
      Lv     = Lw + Kd**2 * Nmax
      Lz     = Lv + Mstar * Kd * Nmax
      Ldmz   = Lz + Mstar * (Nmax + 1)
      Ldelz  = Ldmz + Kd * Nmax
      Ldeldz = Ldelz + Mstar * (Nmax + 1)
      Ldqz   = Ldeldz + Kd * Nmax
      Ldqdmz = Ldqz + Mstar * (Nmax + 1)
      Lrhs   = Ldqdmz + Kd * Nmax
      Lvalst = Lrhs   + Kd * Nmax + Mstar
      Lslope = Lvalst + 4 * Mstar * Nmax
      Laccum = Lslope + Nmax
      Lscl   = Laccum + Nmax + 1
      Ldscl  = Lscl + Mstar * (Nmax + 1)
      Lext1  = Ldscl + Kd * Nmax
      Lext2  = Lext1 + Nmax 
      Lext3  = Lext2 + Nmax + 1
      Lpvtg  = 1
      Lpvtw  = Lpvtg + Mstar * (Nmax + 1)
      Linteg = Lpvtw + Kd * Nmax

C.... If  Iguess .ge. 2, Move  Xiold, Z, And  Dmz  To Their Proper
C.... Locations In  Fspace.

      If ( Iguess .lt. 2 )                          Go To 160
      Nold = N
      If (Iguess .eq. 4)  Nold = Ispace(1)
      Nz = Mstar * (Nold + 1)
      Ndmz = Kd * Nold
      Np1 = N + 1
      If ( Iguess .eq. 4 )  Np1 = Np1 + Nold + 1
      Do 120 I = 1,Nz
  120 Fspace( Lz+I-1 )  =  Fspace( Np1+I )
      Idmz = Np1 + Nz
      Do 125 I = 1,Ndmz
  125 Fspace( Ldmz+I-1 )  =  Fspace( Idmz+I )
      Np1 = Nold + 1
      If ( Iguess .eq. 4 )                          Go To 140
      Do 130 I = 1,Np1
  130 Fspace( Lxiold+I-1 )  =  Fspace( Lxi+I-1 )
      Go To 160
  140 Do 150 I = 1,Np1
  150 Fspace( Lxiold+I-1 )  =  Fspace( N+1+I )
  160 Continue

C.... Initialize Collocation Points And Constants.

      Call Consts ( K, Rho, Coef )

C.... Initialize Mesh.

      Call Newmsh (3+Iread, Fspace(Lxi), Fspace(Lxiold), Dummy,
     +  Dummy, Dummy, Dummy, Dummy, Nfxpnt, Fixpnt, Dummy, Dummy)

C.... Determine First Approximation, If The Problem Is Nonlinear.

      If (Iguess .ge. 2)                            Go To 200
      Np1 = N + 1
      Do 170 I = 1, Np1
 170     Fspace( I + Lxiold - 1 ) = Fspace( I + Lxi - 1 )
      Nold = N
      If ( Nonlin .eq. 0  .or. Iguess .eq. 1 )      Go To 200

C.... System Provides First Approximation Of The Solution.
C.... Choose Z(J) = 0  For J = 1,...,Mstar.

      Do 180 I = 1, Nz
 180     Fspace( Lz-1+I ) = 0.5d0
      Do 190 I = 1, Ndmz
 190     Fspace( Ldmz-1+I ) = 0.5d0
 200  Continue
      If (Iguess .ge. 2)  Iguess = 0

C.... *****************************************************************
C.... Initialisation And Explanation Of Variables Used In The 
C.... Continuation Strategy Begin Here.
C.... *****************************************************************

C.... Save Details Of Initial Mesh And Initial Guess In Case We Need 
C.... To Backtrack

      Nbk1 = N
      Nbk2 = Nold
      Do 202, J = 1,Nbk1+1
         Fspace(Lext3 + J - 1) = Fspace(J)
 202  Continue
      If (Nonlin .eq. 1) Then
         Indx = Lext3 + Nbk1
         Do 204, J = 1,Nbk2+1
            Fspace(Indx + J) = Fspace(Lxiold+J-1)
 204     Continue
         Indx = Indx + Nbk2 + 1
         Nz = Mstar * (Nbk2 + 1)
         Do 206, J = 1,Nz
            Fspace(Indx + J) = Fspace(Lz+J-1)
 206     Continue
         Indx = Indx + Nz
         Ndmz = Kd * Nbk2
         Do 208, J = 1,Ndmz
            Fspace(Indx + J) = Fspace(Ldmz+J-1)
 208     Continue
      Endif

C.... Initialise Extrapolation Variables

      Do 210 J = 1,3
        Phi(J) = Zero
        E(J) = Zero
 210  Continue

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
C.... Attempted. If Iback = 1 Then We Do Not Backtrack.

      Ifinal = 0
      Iback = 0

C.... Istep Is A Flag Which Limits The Size Of The Continuation Steps
C.... After Experiencing A Failure For Some Step. 

      Istep = 0

C.... Ncs Counts The Total Number Of Continuation Steps Taken.
C.... Nss Counts The Number Of Successful Continuation Steps Taken.

      Ncs = 0
      Nss = 0

C.... Idc Counts The Number Of Consecutive Steps For Which The Maximum
C.... Value Of The Monitor Function Is Decreasing.
C.... Iextrap Is A Flag That Indicates Whether Monitor Function
C.... Extrapolation Is Possible.
C.... Istuk Is A Flag Which Is Set To One When The Continuation 
C.... Algorithm Is Backtracking.
C.... Iprec Is A Flag Used To Indicate Whether The Bounds Imposed By 
C.... Machine Precision Are Being Approached Or Passed.
C.... Epsp Is The Value Of Eps Beyond Which The Machine Precision Is
C.... (Possibly) Not Sufficient To Solve The Given Problem.

      Idc = 0
      Iextrap = 0
      Istuk = 0
      Iprec = 0
      Epsp = Zero

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
      Phimax = 1.0d+2/Precis
      Do 220 J = 1,Ipar(4)
        Tolmax = Min(Tol(J),Tolmax)
 220  Continue
      Phimax = Phimax*Tolmax 
      Phiaim = Zero
      Phialt = Zero
      Epold = Zero
      Hsmlpv = Zero
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
C.... The Algorithm Loops Back To Line 230 For Every Continuation Step.
C.... *****************************************************************

 230  Eps = Max(One/Ep,Epsmin)
      Ep = One/Eps

C.... Solve For Epsmin Exactly

      If (Eps .lt. (1.00001d0*Epsmin)) Ifinal = 1

C.... Increase The Continuation Step Counter Ncs.

      Ncs = Ncs+1
      If (Iprint .lt. 1) Write(Iout,1011) Ncs, Eps 
      Iflag = 1

C.... If We Have Reached The Maximum Number Of Continuation Steps
C.... Then This Will Be The Last Problem We Attempt.

      If (Ncs .eq. Maxcon) Then
         If (Iprint .lt. 1) Write(Iout,1012) Maxcon  
         Iback = 1
         Ifinal = 1
      Endif

C.... Iatt Is A Flag That Governs Whether Mesh Selection Is Performed 
C.... Or Control Is Returned To The Driver To Select A New Parameter.

      Iatt = -1

C.... Attempt To Solve The Latest Continuation Problem.

      Call Contrl (Fspace(Lxi),Fspace(Lxiold),Fspace(Lz),Fspace(Ldmz),
     +     Fspace(Lrhs),Fspace(Ldelz),Fspace(Ldeldz),Fspace(Ldqz),
     +     Fspace(Ldqdmz),Fspace(Lg),Fspace(Lw),Fspace(Lv),
     +     Fspace(Lvalst),Fspace(Lslope),Fspace(Lscl),Fspace(Ldscl),
     +     Fspace(Laccum),Ispace(Lpvtg),Ispace(Linteg),Ispace(Lpvtw),
     +     Nfxpnt,Fixpnt,Iflag,Fsub,Dfsub,Gsub,Dgsub,Guess,
     +     Fspace(Lext1),Fspace(Lext2),Eps)

C.... *****************************************************************
C.... The Logic For A Successful Continuation Step Starts Here
C.... *****************************************************************

      If (Iflag .eq. 1) Then

         Nss = Nss + 1

C.... If The Problem Eps = Epsmin Has Been Solved Succesfully Then
C.... We May Finish. Ifinal = 1 When Eps = Epsmin.

         If (Ifinal .eq. 1) Then 
            If (Iprint .lt. 1) Write(Iout,1013) Eps
            Goto 340
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
            If (Iprint .lt. 1) Write(Iout,1014) Eps
            Epsp = Eps
            Iprec = 1
         Endif
 
C.... Save Details Of Last Problem Solved Successfully In Case We
C.... Need To Backtrack

         Nbk1 = N
         Nbk2 = Nold
         Do 240, J = 1,Nbk1+1
            Fspace(Lext3 + J - 1) = Fspace(J)
 240     Continue
         If (Nonlin .eq. 1) Then
            Indx = Lext3 + Nbk1
            Do 250, J = 1,Nbk2+1
               Fspace(Indx + J) = Fspace(Lxiold+J-1)
 250        Continue
            Indx = Indx + Nbk2 + 1
            Nz = Mstar * (Nbk2 + 1)
            Do 260, J = 1,Nz
               Fspace(Indx + J) = Fspace(Lz+J-1)
 260        Continue
            Indx = Indx + Nz
            Ndmz = Kd * Nbk2
            Do 270, J = 1,Ndmz
               Fspace(Indx + J) = Fspace(Ldmz+J-1)
 270        Continue
         Endif

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
            Goto 290
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
               Goto 290
            Else
               Idc = 1
               E(3) = Ep
               Phi(3) = Phit
            Endif
        
C.... Otherwise Update Extrapolation Data

         Else
            Idc = 0
            Do 280 J = 1,2
               E(J) = E(J+1)
               Phi(J) = Phi(J+1)
 280        Continue
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
         If (Dele .lt. 0.01d0*E(3) .or. (N .gt. (Nmax/2))) Then
            Ep = E(3)
            Epsmin = One/Ep
            If (Iprint .lt. 1) Then
               If (Dele .lt. 0.01d0*E(3)) Then
                  Write(Iout,1015) Epsmin
               Else
                  Write(Iout,1016) Epsmin
               Endif
            Endif
            Emin = Ep
            Iback = 1
            Goto 290      
         Endif
         
C.... The Following Section Of Code Calculates The Desired Value 
C.... (Phiaim) That We Would Like The Maximum Value Of The Monitor 
C.... Function To Take.
  
         Itru = 0
         If (Phit .gt. Pmax(1) .and. Pmax(2) .ne. Phit) Itru = 1
         If (Itru .eq. 1) Amax = -Phit*Phit/(Two*C1h)
         Bmax = Phit*H1
         If (Nold .gt. 3*N/2) Then
            Ccm = Max(Bmax-Three,Bmax/1.5d0)
            Cmax = Min(Cmax,Ccm)
         Endif
         If (Itru .eq. 0) Amax = Max(1.5d0*Bmax,Bmax+Three)
         If (Nss .ge. 2) Hrat = Max(H1/Hsmlpv,One)
         Hsmlpv = Hsml
         Bbm = Max(1.5d0*Bmax,Bmax+Three)
         If (Nonlin .eq. 1) Then
            Fmax = Dble(Itsaim)*Bmax/Dble(Nits)
            If (Fmax .gt. Bmax) Fmax = Max(Fmax,Bbm)
            If (Fmax .lt. Bmax) Fmax = Max(Bmax-Three,Bmax/1.5d0,Fmax)
         Endif
         Htot = Min(Amax,Bbm,Cmax)
         If (Nonlin .eq. 0) Fmax = Htot
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
         
         If ( (Irest .eq. 1 .or. Istep .eq. 1) .and. Nss .ne. 1) Then
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
         
 290     Continue

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
            If (Iprint .lt. 1) Then
               If (Epsp .ne. Zero) Then
                  Write(Iout,1017) Eps,Iflag,Eps,Epsp
               Else
                  Write(Iout,1018) Eps,Iflag,Eps
               Endif
            Endif
            Goto 340
         Endif
         
C.... If Iprec = 2, Then We Know That We Cannot Define A Mesh On
C.... Which The Current Eps Value Can Be Solved To The Requested
C.... Tolerances. We Alter Epsmin Accordingly.

         If (Iprec .eq. 2) Then
            Epsmin = One/(Max((Ep+E(3))/Two,0.9d0*Ep))
         Endif

C.... Insert Details For Backtracking

         If (Iprint .lt. 1) Write(Iout,1019) 
         Ifinal = 0
         N = Nbk1
         Ep = (Ep+E(3))/Two
         If ((Ep-E(3)) .lt. 0.01d0*E(3)) Then
            Ep = E(3)
            Epsmin = One/Ep
            If (Iprint .lt. 1) Write(Iout,1015) Epsmin
            Emin = Ep
            Iback = 1
            Iflag = 1
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

         N = Nbk1 
         Nold = Nbk2
         Do 300, J = 1,N+1
            Fspace(J) = Fspace(Lext3+J-1)
 300     Continue
         If (Nonlin .eq. 1) Then
            Indx = Lext3 + N
            Do 310, J = 1,Nbk2+1
               Fspace(Lxiold+J-1) = Fspace(Indx+J)
 310        Continue
            Indx = Indx + Nbk2 + 1
            Nz = Mstar * (Nbk2 + 1)
            Do 320, J = 1,Nz
               Fspace(Lz+J-1) = Fspace(Indx + J)
 320        Continue
            Indx = Indx + Nz
            Ndmz = Kd * Nbk2
            Do 330, J = 1,Ndmz
               Fspace(Ldmz+J-1) = Fspace(Indx + J)
 330        Continue
         Endif
         
      Endif   
C
C.... *****************************************************************
C.... End Of Logic For An Unsuccessful Continuation Step.
C.... *****************************************************************
C
C.... Set Continuation Variables For Colmod
C
      Iguess = 0
      Iread = 1
      Nz = Mstar * (N + 1)
      Ndmz = Kd * N

C.... *****************************************************************
C.... The Program Loops Back To Line 230 For Every Continuation Step.
C.... *****************************************************************

      Goto 230

C.... *****************************************************************
C.... Before Exiting, Prepare Output.
C.... *****************************************************************

 340  Ispace(1) = N
      Ispace(2) = K
      Ispace(3) = Ncomp
      Ispace(4) = Mstar
      Ispace(5) = Mmax
      Ispace(6) = Nz + Ndmz + N + 2
      K2 = K * K
      Ispace(7) = Ispace(6) + K2 - 1
      Do 350 I = 1, Ncomp
 350     Ispace(7+I) = M(I)
      Do 360 I = 1, Nz
 360     Fspace( N+1+I ) = Fspace( Lz-1+I )
      Idmz = N + 1 + Nz
      Do 370 I = 1, Ndmz
 370     Fspace( Idmz+I ) = Fspace( Ldmz-1+I )
      Ic = Idmz + Ndmz
      Do 380 I = 1, K2
 380     Fspace( Ic+I ) = Coef(I)

      If (Iflag .eq. 1 .and. Iprint .eq. 0) Then
         Write(Iout,1020) ('Z',Ltol(J), J=1,Ntol)
         Write(Iout,*)
         Jstep = Max(N/30,1)
         Iadd = Jstep*Mstar
         Ind = N + 1 
         Do 390 I = 1, N, Jstep
            Write(Iout,1021) I,Fspace(I),(Fspace(Ind+Ltol(J)),J=1,Ntol)
            Ind = Ind + Iadd
 390     Continue
         Ind = N+1 + Nz - Mstar
         Write(Iout,1021) N+1,Fspace(N+1),(Fspace(Ind+Ltol(J)),J=1,Ntol)
         Write(Iout,1022)
      Endif

      Return
C-----------------------------------------------------------------------
 1000 Format(/ 37h The Number Of (Linear) Diff Eqns Is , I3/,
     +     17h Their Orders Are, 20i3)
 1001 Format(/ 40h The Number Of (Nonlinear) Diff Eqns Is , I3/,
     +     17h Their Orders Are, 20i3)
 1002 Format(27h Side Condition Points Zeta, 8f10.6, 4( / 27x, 8f10.6))
 1003 Format(10h There Are ,I5,27h Fixed Points In The Mesh - ,
     +     10(6f10.6/))
 1004 Format(37h Number Of Colloc Pts Per Interval Is, I3)
 1005 Format(39h Components Of Z Requiring Tolerances -,8(7x,I2,1x),
     +     4(/38x,8i10))
 1006 Format(33h Corresponding Error Tolerances -,6x,8d10.2,
     +     4(/39x,8d10.2))
 1007 Format(44h Initial Mesh(Es) And Z,Dmz Provided By User)
 1008 Format (33h The Initial Value Of Epsilon Is ,D11.4,/,
     +     39h The Desired Final Value Of Epsilon Is ,D11.4)
 1009 Format(44h The Maximum Number Of Subintervals Is Min (, I4,
     +     23h (Allowed From Fspace),,I4, 24h (Allowed From Ispace) ))
 1010 Format(/53h Insufficient Space To Double Mesh For Error Estimate)
 1011 Format (/,1x,82('*'),/,' Continuation Step ',I2,
     +        ', Epsilon = ',D9.4,/,1x,82('*'))
 1012 Format (/,1x,'This Is The Final Continuation Step Since ',
     +        'Maxcon =',I4)
 1013 Format (/,1x,82('$'),//,' Tolerances Satisfied For Final Problem',
     +        ' Epsilon = ',D9.4,//,1x,82('$'),/)
 1014 Format(/,' ** Machine Precision (Possibly) Not Sufficient For',
     +        ' Epsilon Less Than ',D9.4) 
 1015 Format (/,' Continuation Steps Too Small, Change Epsmin To ',D9.4)
 1016 Format (/,' Storage Limit Being Approached, Change Epsmin To ',
     +     D9.4)
 1017 Format (/,1x,82('$'),//,' Final Problem Epsilon = ',D10.4,/,
     +        ' Not Solved, Iflag =',I3,/, ' Try Running Problem Again',
     +        ' With Epsmin Greater Than ',D9.4,/, ' Machine Precision',
     +        ' (Possibly) Not Sufficient For Eps Less Than ',D9.4,//,
     +        1x,82('$'),/) 
 1018 Format (/,1x,82('$'),//,' Final Problem Epsilon = ',D10.4,/,
     +        ' Not Solved, Iflag =',I3,/, ' Try Running Problem Again',
     +        ' With Epsmin Greater Than ',D9.4,//,1x,82('$'),/) 
 1019 Format (/,' ** Failed Step - Bactracking For Larger Value ',
     +          'Of Epsilon')
 1020 Format(' The Final Mesh And Solution Components Are:',//,
     +       5x,'I',10x,'X(I) ',40(14x,A,'(',I2,')'))
 1021 Format(I6,41d19.8)
 1022 Format(/,1x,82('$'),/)
C-----------------------------------------------------------------------
      End 


      Subroutine Contrl (Xi, Xiold, Z, Dmz, Rhs, Delz, Deldmz,
     +           Dqz, Dqdmz, G, W, V, Valstr, Slope, Scale, Dscale,
     +           Accum, Ipvtg, Integs, Ipvtw, Nfxpnt, Fixpnt, Iflag,
     +           Fsub, Dfsub, Gsub, Dgsub, Guess, Slpold, Voldmsh, Eps )
C
C**********************************************************************
C
C   Purpose
C     This Subroutine Is The Actual Driver.  The Nonlinear Iteration
C     Strategy Is Controlled Here ( See [4] ). Upon Convergence, Errchk
C     Is Called To Test For Satisfaction Of The Requested Tolerances.
C
C   Variables
C
C     Check  - Maximum Tolerance Value, Used As Part Of Criteria For
C              Checking For Nonlinear Iteration Convergence
C     Relax  - The Relaxation Factor For Damped Newton Iteration
C     Relmin - Minimum Allowable Value For Relax  (Otherwise The
C              Jacobian Is Considered Singular).
C     Rlxold - Previous Relax
C     Rstart - Initial Value For Relax When Problem Is Sensitive
C     Ifrz   - Number Of Fixed Jacobian Iterations
C     Lmtfrz - Maximum Value For Ifrz Before Performing A Reinversion
C     Iter   - Number Of Iterations (Counted Only When Jacobian
C              Reinversions Are Performed).
C     Xi     - Current Mesh
C     Xiold  - Previous Mesh
C     Ipred  = 0  If Relax Is Determined By A Correction
C            = 1  If Relax Is Determined By A Prediction
C     Ifreez = 0  If The Jacobian Is To Be Updated
C            = 1  If The Jacobian Is Currently Fixed (Frozen)
C     Iconv  = 0  If No Previous Convergence Has Been Obtained
C            = 1  If Convergence On A Previous Mesh Has Been Obtained
C     Rnorm  - Norm Of Rhs (Right Hand Side) For Current Iteration
C     Rnold  - Norm Of Rhs For Previous Iteration
C     Anscl  - Scaled Norm Of Newton Correction
C     Anfix  - Scaled Norm Of Newton Correction At Next Step
C     Anorm  - Scaled Norm Of A Correction Obtained With Jacobian Fixed
C     Nz     - Number Of Components Of  Z  (See Subroutine Approx)
C     Ndmz   - Number Of Components Of  Dmz  (See Subroutine Approx)
C     Imesh  - A Control Variable For Subroutines Newmsh And Errchk
C            = 1  The Current Mesh Resulted From Mesh Selection
C                 Or Is The Initial Mesh.
C            = 2  The Current Mesh Resulted From Doubling The
C                 Previous Mesh
C
C**********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Xi(*), Xiold(*), Z(*), Dmz(*), Rhs(*)
      Dimension G(*), W(*), V(*), Valstr(*), Slope(*), Accum(*)
      Dimension Delz(*), Deldmz(*), Dqz(*), Dqdmz(*) , Fixpnt(*)
      Dimension Dummy(1), Scale(*), Dscale(*)
      Dimension Integs(*), Ipvtg(*), Ipvtw(*)
      Dimension Slpold(*), Voldmsh(*)
C
      Common /Colout/ Precis, Iout, Iprint
      Common /Colord/ K, Ncomp, Mstar, Kd, Mmax, M(20)
      Common /Colapr/ N, Nold, Nmax, Nz, Ndmz, Mshflg
      Common /Colnln/ Nonlin, Iter, Limit, Iguess
      Common /Colest/ Wgtmsh(40), Wgterr(40), Tolin(40),
     +                Root(40), Jtol(40), Ltol(40), Ntol
C
      Common /Flags/ Ifinal,Iatt,Iback,Iprec
      Common /Convg/ Nits
C
      External Fsub, Dfsub, Gsub, Dgsub, Guess

*  The Parameter Inumb Is A Counter Used To Limit To Three The Number 
*  Of Mesh Selections Performed For The Final Continuation Problem. 

      Inumb = 0

      Relmin = 1.d-3
      Rstart = 1.d-2
      Lmtfrz = 4
C
C.... Compute The Maximum Tolerance
C
      Check = 0.d0
      Do 10 I = 1, Ntol
   10   Check = Max ( Tolin(I), Check )
      Imesh = 1
      Iconv = 0
      If ( Nonlin .eq. 0 ) Iconv = 1
      Icor = 0
      Msing = 0
C
C.... The Main Iteration Begins Here .
C.... Loop 20 Is Executed Until Error Tolerances Are Satisfied Or
C.... The Code Fails (Due To A Singular Matrix Or Storage Limitations)
C
 20   Continue
C
C.... Initialization For A New Mesh
C     
      Iter = 0
      If ( Nonlin .gt. 0 )                     Go To 50
C     
C.... The Linear Case.
C.... Set Up And Solve Equations
C     
      Call Lsyslv (Msing, Xi, Xiold, Dummy, Dummy, Z, Dmz, G,
     +     W, V, Rhs, Dummy, Integs, Ipvtg, Ipvtw, Rnorm, 0,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess, Eps )
C     
C.... Check For A Singular Matrix
C     
      If ( Msing .eq. 0 )                      Go To 400
 30   If ( Msing .lt. 0 )                      Go To 40
      If ( Iprint .lt. 1 )  Write (Iout,1000)
      Go To 460
 40   If ( Iprint .lt. 1 )  Write (Iout,1001)
      Iflag = 0
      Return
C     
C.... Iteration Loop For Nonlinear Case
C.... Define The Initial Relaxation Parameter ( =  Relax)
C     
 50   Relax = 1.d0
C     
C.... Check For Previous Convergence
C     
      If ( Iconv .eq. 0 )                      Go To 160
C     
C.... Convergence On A Previous Mesh Has Been Obtained.    Thus
C.... We Have A Very Good Initial Approximation For The Newton
C.... Process.    Proceed With One Full Newton And Then Iterate
C.... With A Fixed Jacobian.
C     
      Ifreez = 0
C     
C.... Evaluate Right Hand Side And Its Norm  And
C.... Find The First Newton Correction
C     
      Call Lsyslv (Msing, Xi, Xiold, Z, Dmz, Delz, Deldmz, G,
     +     W, V, Rhs, Dqdmz, Integs, Ipvtg, Ipvtw, Rnold, 1,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess, Eps )
C     
      If ( Iprint .lt. 0 )  Write(Iout,1002)
      If ( Iprint .lt. 0 )  Write (Iout,1003) Iter, Rnold
      Go To 70
C     
C.... Solve For The Next Iterate .
C.... The Value Of Ifreez Determines Whether This Is A Full
C.... Newton Step ( = 0) Or A Fixed Jacobian Iteration ( = 1).
C     
 60   If ( Iprint .lt. 0 )  Write (Iout,1003) Iter, Rnorm
      Rnold = Rnorm
      Call Lsyslv (Msing, Xi, Xiold, Z, Dmz, Delz, Deldmz, G,
     +     W, V, Rhs, Dummy, Integs, Ipvtg, Ipvtw, Rnorm,
     +     3+Ifreez, Fsub, Dfsub, Gsub, Dgsub, Guess, Eps )
C     
C.... Check For A Singular Matrix
C     
 70   If ( Msing .ne. 0 )                      Go To 30
      If ( Ifreez .eq. 1 )                     Go To 80
C     
C.... A Full Newton Step
C     
      Iter = Iter + 1
      Ifrz = 0
 80   Continue
C     
C.... Update   Z And Dmz , Compute New  Rhs  And Its Norm
C     
      Do 90 I = 1, Nz
         Z(I) = Z(I) + Delz(I)
 90   Continue
      Do 100 I = 1, Ndmz
         Dmz(I) = Dmz(I) + Deldmz(I)
 100  Continue
      Call Lsyslv (Msing, Xi, Xiold, Z, Dmz, Delz, Deldmz, G,
     +     W, V, Rhs, Dummy, Integs, Ipvtg, Ipvtw, Rnorm, 2,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess, Eps )
C     
C.... Check Monotonicity. If The Norm Of  Rhs  Gets Smaller,
C.... Proceed With A Fixed Jacobian; Else Proceed Cautiously,
C.... As If Convergence Has Not Been Obtained Before (Iconv = 0).
C     
      If ( Rnorm .lt. Precis )                 Go To 390
      If ( Rnorm .gt. Rnold )                  Go To 130
      If ( Ifreez .eq. 1 )                     Go To 110
      Ifreez = 1
      Go To 60
C     
C.... Verify That The Linear Convergence With Fixed Jacobian
C.... Is Fast Enough.
C     
 110  Ifrz = Ifrz + 1
      If ( Ifrz .ge. Lmtfrz )       Ifreez = 0
      If ( Rnold .lt. 4.d0*Rnorm )  Ifreez = 0
C     
C.... Check Convergence (Iconv = 1).
C     
      Do 120 It = 1, Ntol
         Inz = Ltol(It)
         Do 120 Iz = Inz, Nz, Mstar
            If ( Abs(Delz(Iz))  .gt. 
     +           Tolin(It) * (Abs(Z(Iz)) + 1.d0))  Go To 60
 120  Continue
C     
C.... Convergence Obtained
C     
      If ( Iprint .eq. -1 )  Write (Iout,1004) Iter
      If ( Iprint .eq. 0 )  Write (Iout,1016) Iter
      If (Iatt .eq. -1) Nits = Iter
      Go To 400
C     
C.... Convergence Of Fixed Jacobian Iteration Failed.
C     
 130  If ( Iprint .lt. 0 )  Write (Iout,1003) Iter, Rnorm
      If ( Iprint .lt. 0 )  Write (Iout,1005)
      Iconv = 0
      Relax = Rstart
      Do 140 I = 1, Nz
         Z(I) = Z(I) - Delz(I)
 140  Continue
      Do 150 I = 1, Ndmz
         Dmz(I) = Dmz(I) - Deldmz(I)
 150  Continue
C     
C.... Update Old Mesh
C     
      Np1 = N + 1
      Do 155 I = 1, Np1
         Xiold(I) = Xi(I)
 155  Continue
      Nold = N
C     
      Iter = 0
C     
C.... No Previous Convergence Has Been Obtained. Proceed
C.... With The Damped Newton Method.
C.... Evaluate Rhs And Find The First Newton Correction.
C     
 160  If(Iprint .lt. 0)  Write (Iout,1006)
      Call Lsyslv (Msing, Xi, Xiold, Z, Dmz, Delz, Deldmz, G,
     +     W, V, Rhs, Dqdmz, Integs, Ipvtg, Ipvtw, Rnold, 1,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess, Eps )
C     
C.... Check For A Singular Matrix
C     
      If ( Msing .ne. 0 )                       Go To 30
C     
C.... Bookkeeping For First Mesh
C     
      If ( Iguess .eq. 1 )  Iguess = 0
C     
C.... Find Initial Scaling
C     
      Call Skale (N, Mstar, Kd, Z, Xi, Scale, Dscale)
      Go To 220
C     
C.... Main Iteration Loop
C     
 170  Rnold = Rnorm
      If ( Iter .ge. Limit )                   Go To 430
C     
C.... Update Scaling
C     
      Call Skale (N, Mstar, Kd, Z, Xi, Scale, Dscale)
C     
C.... Compute Norm Of Newton Correction With New Scaling
C     
      Anscl = 0.d0
      Do 180 I = 1, Nz
         Anscl = Anscl + (Delz(I) * Scale(I))**2
 180  Continue
      Do 190 I = 1, Ndmz
         Anscl = Anscl + (Deldmz(I) * Dscale(I))**2
 190  Continue
      Anscl = Sqrt(Anscl / Dfloat(Nz+Ndmz))
C     
C.... Find A Newton Direction
C     
      Call Lsyslv (Msing, Xi, Xiold, Z, Dmz, Delz, Deldmz, G,
     +     W, V, Rhs, Dummy, Integs, Ipvtg, Ipvtw, Rnorm, 3,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess, Eps )
C     
C.... Check For A Singular Matrix
C     
      If ( Msing .ne. 0 )                      Go To 30
C     
C.... Predict Relaxation Factor For Newton Step.
C     
      Andif = 0.d0
      Do 200 I = 1, Nz
         Andif = Andif + ((Dqz(I) - Delz(I)) * Scale(I))**2
 200  Continue
      Do 210 I = 1, Ndmz
         Andif = Andif + ((Dqdmz(I) - Deldmz(I)) * Dscale(I))**2
 210  Continue
      Andif = Sqrt(Andif/Dfloat(Nz+Ndmz) + Precis)
      Relax = Relax * Anscl / Andif
      If ( Relax .gt. 1.d0 )  Relax = 1.d0
 220  Rlxold = Relax
      Ipred = 1
      Iter = Iter + 1
C     
C.... Determine A New  Z And Dmz  And Find New  Rhs  And Its Norm
C     
      Do 230 I = 1, Nz
         Z(I) = Z(I)  +  Relax * Delz(I)
 230  Continue
      Do 240 I = 1, Ndmz
         Dmz(I) = Dmz(I)  +  Relax * Deldmz(I)
 240  Continue
 250  Call Lsyslv (Msing, Xi, Xiold, Z, Dmz, Dqz, Dqdmz, G,
     +     W, V, Rhs, Dummy, Integs, Ipvtg, Ipvtw, Rnorm, 2,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess, Eps )
C     
C.... Compute A Fixed Jacobian Iterate (Used To Control Relax)
C     
      Call Lsyslv (Msing, Xi, Xiold, Z, Dmz, Dqz, Dqdmz, G,
     +     W, V, Rhs, Dummy, Integs, Ipvtg, Ipvtw, Rnorm, 4,
     +     Fsub, Dfsub, Gsub, Dgsub, Guess, Eps )
C     
C.... Find Scaled Norms Of Various Terms Used To Correct Relax
C     
      Anorm = 0.d0
      Anfix = 0.d0
      Do 260 I = 1, Nz
         Anorm = Anorm  +  (Delz(I) * Scale(I))**2
         Anfix = Anfix  +  (Dqz(I) * Scale(I))**2
 260  Continue
      Do 270 I = 1, Ndmz
         Anorm = Anorm  +  (Deldmz(I) * Dscale(I))**2
         Anfix = Anfix  +  (Dqdmz(I) * Dscale(I))**2
 270  Continue
      Anorm = Sqrt(Anorm / Dfloat(Nz+Ndmz))
      Anfix = Sqrt(Anfix / Dfloat(Nz+Ndmz))
      If ( Icor .eq. 1 )                         Go To 280
      If (Iprint .lt. 0)  Write (Iout,1007) Iter, Relax, Anorm,
     +     Anfix, Rnold, Rnorm
      Go To 290
 280  If (Iprint .lt. 0) Write (Iout,1008) Relax, Anorm, Anfix,
     +     Rnold, Rnorm
 290  Icor = 0
C     
C.... Check For Monotonic Decrease In  Delz And Deldmz.
C     
      If (Anfix .lt. Precis .or. Rnorm .lt. Precis)  Go To 390
      If ( Anfix .gt. Anorm )                    Go To 300
C     
C.... We Have A Decrease.
C.... If  Dqz  And Dqdmz  Small, Check For Convergence
C     
      If ( Anfix .le. Check )                    Go To 350
C     
C.... Correct The Predicted  Relax  Unless The Corrected
C.... Value Is Within 10 Percent Of The Predicted One.
C     
      If ( Ipred .ne. 1 )                        Go To 170
 300  If ( Iter .ge. Limit )                     Go To 430
C     
C.... Correct The Relaxation Factor.
C     
      Ipred = 0
      Arg = (Anfix/Anorm - 1.d0) / Relax + 1.d0
      If ( Arg .lt. 0.d0 )                       Go To 170
      If (Arg .le. .25d0*Relax+.125d0*Relax**2)  Go To 310
      Factor = -1.d0 + Sqrt (1.d0+8.d0 * Arg)
      If ( Abs(Factor-1.d0) .lt. .1d0*Factor )  Go To 170
      If ( Factor .lt. 0.5d0 )  Factor = 0.5d0
      Relax = Relax / Factor
      Go To 320
 310  If ( Relax .ge. .9d0 )                     Go To 170
      Relax = 1.d0
 320  Icor = 1
      If ( Relax .lt. Relmin )                   Go To 440
      Fact = Relax - Rlxold
      Do 330 I = 1, Nz
         Z(I) = Z(I)  +  Fact * Delz(I)
 330  Continue
      Do 340 I = 1, Ndmz
         Dmz(I) = Dmz(I)  +  Fact * Deldmz(I)
 340  Continue
      Rlxold = Relax
      Go To 250
C     
C.... Check Convergence (Iconv = 0).
C     
 350  Continue
      Do 360 It = 1, Ntol
         Inz = Ltol(It)
         Do 360 Iz = Inz, Nz, Mstar
            If ( Abs(Dqz(Iz))  .gt. 
     +           Tolin(It) * (Abs(Z(Iz)) + 1.d0) )   Go To 170
 360  Continue
C     
C.... Convergence Obtained
C     
      If ( Iprint .eq. -1 )  Write (Iout,1004) Iter
      If ( Iprint .eq. 0 )  Write (Iout,1016) Iter
      If (Iatt .eq. -1) Nits = Iter
C     
C.... Since Convergence Obtained, Update  Z And Dmz  With Term
C.... From The Fixed Jacobian Iteration.
C     
      Do 370 I = 1, Nz
         Z(I) = Z(I)  +  Dqz(I)
 370  Continue
      Do 380 I = 1, Ndmz
         Dmz(I) = Dmz(I)  +  Dqdmz(I)
 380  Continue
 390  If (Anfix .lt. Precis .or. Rnorm .lt. Precis) Then
         If (Iatt .eq. -1) Nits = Iter
         If ( Iprint .eq. -1 )  Write (Iout,1004) Iter
         If ( Iprint .eq. 0 )  Write (Iout,1016) Iter
      Endif

      Iconv = 1
C     
C.... If Full Output Has Been Requested, Print Values Of The
C.... Solution Components   Z  At The Meshpoints.
C     
 400  If ( Iprint .ge. 0 )                     Go To 420
      Do 410 J = 1, Mstar
         Write(Iout,1009) J
         Write(Iout,1010) (Z(Lj), Lj = J, Nz, Mstar)
 410  Continue
      
C.... Check For Error Tolerance Satisfaction

 420  Ifin = 1
      If (Imesh .eq. 2) Call Errchk (Xi, Z, Dmz, Valstr, Ifin)
      If ( Imesh .eq. 1  .or. Ifin .eq. 0 )     Go To 460
C     
      Iflag = 1
      Return
C     
C.... Diagnostics For Failure Of Nonlinear Iteration.
C     
 430  If ( Iprint .lt. 1 )  Write (Iout,1011) Iter
      Go To 450
 440  If( Iprint .lt. 1 )  Write(Iout,1012) Relax, Relmin
 450  Iflag = -2
      Return
C     
C.... Update Old Mesh
C     
 460  Np1 = N + 1
      Do 470 I = 1, Np1
         Xiold(I) = Xi(I)
 470  Continue
      Nold = N
         
C.... Pick A New Mesh
      
      Imesh = 1
      
      If (Ifinal .eq. 1) Then
         If (Iatt .ge. 1) Then
            If (N .gt. Nmax/2) Then
               Iflag = -1
               Return
            Else
               Imesh = 2
               Inumb = Inumb+1
            Endif
            If (Inumb .eq. 4 .and. Iback .eq. 0) Then
               Iflag = -1
               Return
            Endif
         Endif
      Else
         If (Iatt .eq. 0) Then
            Call Newmsh (Imesh, Xi, Xiold, Z, Dmz, Valstr,
     +           Slope, Accum, Nfxpnt, Fixpnt, Slpold, Voldmsh)
            If (Iprec .eq. 2) Then
               If (Iprint .lt. 1) Write(Iout,1013)
               Iflag = -1
            Else
               Iflag = 1
            Endif
            Return
         Endif
      Endif
      
      Call Newmsh (Imesh, Xi, Xiold, Z, Dmz, Valstr,
     +     Slope, Accum, Nfxpnt, Fixpnt, Slpold, Voldmsh)
      If (Iprec .eq. 2) Then
         If (Iprint .lt. 1) Write(Iout,1013)
         Iflag = -1
         Return
      Endif
      
      Iatt = Iatt+1
      
C.... Exit If Expected N Is Too Large (But May Try N = Nmax Once)
      
      If ( N .le. Nmax )                       Go To 480
      N = N / 2
      Iflag = -1
      If ( Iconv .eq. 0 .and. Iprint .lt. 1 )  Write (Iout,1014)
      If ( Iconv .eq. 1 .and. Iprint .lt. 1 )  Write (Iout,1015)
      Return
 480  If ( Iconv .eq. 0 )  Imesh = 1
      Go To 20
C     ---------------------------------------------------------------
 1000 Format(/40h A Local Elimination Matrix Is Singular )
 1001 Format(/35h The Global Bvp-Matrix Is Singular )
 1002 Format(/27h Fixed Jacobian Iterations,)
 1003 Format(13h Iteration = , I3, 15h  Norm (Rhs) = , D10.2)
 1004 Format(/18h Convergence After , I3,11h Iterations )
 1005 Format(/35h Switch To Damped Newton Iteration,)
 1006 Format(/30h Full Damped Newton Iteration,)
 1007 Format(13h Iteration = ,I3,22h  Relaxation Factor = ,D10.2
     +       /33h Norm Of Scaled Rhs Changes From ,D10.2,3h To,D10.2
     +       /33h Norm   Of   Rhs  Changes  From  ,D10.2,3h To,D10.2,
     +       D10.2)
 1008 Format(40h Relaxation Factor Corrected To Relax = , D10.2
     +       /33h Norm Of Scaled Rhs Changes From ,D10.2,3h To,D10.2
     +       /33h Norm   Of   Rhs  Changes  From  ,D10.2,3h To,D10.2
     +       ,D10.2)
 1009 Format(/19h Mesh Values For Z(, I2, 2h), )
 1010 Format(1h , 8d15.7)
 1011 Format(/22h No Convergence After , I3, 11h Iterations)
 1012 Format(/37h No Convergence.  Relaxation Factor =,D10.3
     +       ,24h Is Too Small (Less Than, D10.3, 1h))
 1013 Format(/,1x,'** Mesh Cannot Be Defined Within The Bounds Imposed',
     +       ' By The Machine Precision')
 1014 Format(18h  (No Convergence) )
 1015 Format(50h  (Probably Tolerances Too Stringent, Or Nmax Too
     +       ,6hsmall) )
 1016 Format(18h Convergence After , I3,11h Iterations )
      End

      Subroutine Skale (N, Mstar, Kd, Z, Xi, Scale, Dscale)
C
C**********************************************************************
C
C   Purpose
C            Provide A Proper Scaling Of The State Variables, Used
C            To Control The Damping Factor For A Newton Iteration [2].
C
C   Variables
C
C            N      = Number Of Mesh Subintervals
C            Mstar  = Number Of Unknomns In Z(U(X))
C            Kd     = Number Of Unknowns In Dmz
C            Z      = The Global Unknown Vector
C            Xi     = The Current Mesh
C            Scale  = Scaling Vector For Z
C            Dscale = Scaling Vector For Dmz
C
C**********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Z(Mstar,*), Scale(Mstar,*), Dscale(Kd,*)
      Dimension Xi(*), Basm(5)
C
      Common /Colord/ K, Ncomp, Id1, Id2, Mmax, M(20)
C
      Basm(1) = 1.d0
      Do 50 J = 1,N
        Iz = 1
        H = Xi(J+1) - Xi(J)
        Do 10 L = 1, Mmax
          Basm(L+1) = Basm(L) * H / Dfloat(L)
  10    Continue
        Do 40 Icomp = 1, Ncomp
          Scal = (Abs(Z(Iz,J)) + Abs(Z(Iz,J+1))) * .5d0 + 1.d0
          Mj = M(Icomp)
          Do 20 L = 1, Mj
            Scale(Iz,J) = Basm(L) / Scal
            Iz = Iz + 1
  20      Continue
          Scal = Basm(Mj+1) / Scal
          Do 30 Idmz = Icomp, Kd, Ncomp
            Dscale(Idmz,J) = Scal
  30      Continue
  40    Continue
  50  Continue
      Np1 = N + 1
      Do 60 Iz = 1, Mstar
        Scale(Iz,Np1) = Scale(Iz,N)
  60  Continue
      Return
      End
C----------------------------------------------------------------------
C                            P A R T  2
C          Mesh Selection, Error Estimation, (And Related
C          Constant Assignment) Routines -- See [3], [4], [6]
C----------------------------------------------------------------------
C
      Subroutine Newmsh (Mode, Xi, Xiold, Z, Dmz, Valstr,
     +                   Slope, Accum, Nfxpnt, Fixpnt, Slpold, Voldmsh)
C
C**********************************************************************
C
C   Purpose
C            Select A Mesh On Which A Collocation Solution Is To Be
C            Determined
C
C                           There Are 4 Possible Modes Of Action:
C            Mode = 4,3 - Deal Mainly With Definition Of An Initial
C                           Mesh For The Current Boundary Value Problem
C                 = 2,1   - Deal With Definition Of A New Mesh, Either
C                           By Simple Mesh Halving Or By Mesh Selection
C            More Specifically, For
C            Mode = 4  An Initial (Generally Nonuniform) Mesh Is
C                      Defined By The User
C                 = 3  A Simple Uniform Mesh (Except Possibly For Some
C                      Fixed Points) Is Defined; N = No. Of Subintervals
C                 = 1  The Automatic Mesh Selection Procedure Is Used
C                      (See [6] For Details)
C                 = 2  A Simple Mesh Halving Is Performed
C
C**********************************************************************
C
C   Variables
C
C            N      = Number Of Mesh Subintervals
C            Nold   = Number Of Subintervals For Former Mesh
C            Xi     - Mesh Point Array
C            Xiold  - Former Mesh Point Array
C            Mshflg = 1  The Mesh Is A Halving Of Its Former Mesh
C                       (So An Error Estimate Has Been Calculated)
C                   = 0  Otherwise
C            Iguess - Ipar(9) In Subroutine Colmod.  It Is Used
C                     Here Only For Mode = 5 And 4, Where
C                   = 2 The Subroutine Sets Xi = Xiold.  This Is
C                       Used E.g. If Continuation Is Being Per-
C                       Formed, And A Mesh For The Old Differen-
C                       Tial Equation Is Being Used
C                   = 3 Same As For  = 2, Except Xi Uses Every Other
C                       Point Of Xiold (So Mesh Xiold Is Mesh Xi
C                       Halved)
C                   = 4 Xi Has Been Defined By The User, And An Old
C                       Mesh Xiold Is Also Available
C                       Otherwise, Xi Has Been Defined By The User
C                       And We Set Xiold = Xi In This Subroutine
C            Slope  - An Approximate Quantity To Be Equidistributed For
C                     Mesh Selection (See [3]), Viz,
C                             .                        (K+Mj)
C                     Slope(I) =      Max   (Weight(L) *U      (Xi(I)))
C                               1 .le. L .le. Ntol         J
C
C                     Where J = Jtol(L)
C            Slphmx - Maximum Of Slope(I)*(Xiold(I+1)-Xiold(I)) For
C                     I = 1 ,..., Nold.
C            Accum  - Accum(I) Is The Integral Of  Slope  From  Aleft
C                     To  Xiold(I).
C            Valstr - Is Assigned Values Needed In  Errchk  For The
C                     Error Estimate.
C**********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension D1(40), D2(40), Slope(*), Accum(*), Valstr(*), D3(40)
      Dimension Xi(*), Xiold(*), Z(*), Dmz(*), Fixpnt(*), Dummy(1)
      Dimension Slpold(*), Voldmsh(*)
      Dimension Pmax(2), Hord(2)
C
      Common /Colout/ Precis, Iout, Iprint
      Common /Colord/ K, Ncomp, Mstar, Kd, Mmax, M(20)
      Common /Colapr/ N, Nold, Nmax, Nz, Ndmz, Mshflg
      Common /Colnln/ Nonlin, Iter, Limit, Iguess
      Common /Colsid/  Zeta(40), Aleft, Aright, Izeta, Idum
      Common /Colbas/ B(28), Acol(28,7), Asave(28,4)
      Common /Colest/ Wgtmsh(40), Wgterr(40), Tolin(40),
     +                Root(40), Jtol(40), Ltol(40), Ntol
      Common /Flags/ Ifinal,Iatt,Iback,Iprec
      Common /Mshvar/ Pmax,Hord,Hsml
C
      Nfxp1 = Nfxpnt +1
      Iprec = Min(Iprec,1)
      Go To (180, 100, 50, 20), Mode
C
C.... Mode = 4   The User-Specified Initial Mesh Is Already In Place.
C
   20 If ( Iguess .lt. 2 )                          Go To 40
C
C.... Iguess = 2, 3 Or 4.
C
      Noldp1 = Nold + 1
      If (Iprint .lt. 0)  Write(Iout,1000) Nold,(Xiold(I), I = 1,Noldp1)
      If ( Iguess .ne. 3 )                          Go To 40
C
C.... If Iread ( Ipar(8) ) .eq. 1 And Iguess ( Ipar(9) ) .eq. 3
C.... Then The First Mesh Is Every Second Point Of The
C.... Mesh In  Xiold .
C
      N = Nold /2
      I = 0
      Do 30 J = 1, Nold, 2
           I = I + 1
   30 Xi(I) = Xiold(J)
   40 Continue
      Np1 = N + 1
      Xi(1) = Aleft
      Xi(Np1) = Aright
      Go To 320
C
C.... Mode = 3   Generate A (Piecewise) Uniform Mesh. If There Are
C.... Fixed Points Then Ensure That The N Being Used Is Large Enough.
C
   50 If ( N .lt. Nfxp1 )  N = Nfxp1
      Np1 = N + 1
      Xi(1) = Aleft
      Ileft = 1
      Xleft = Aleft
C
C.... Loop Over The Subregions Between Fixed Points.
C
      Do 90 J = 1, Nfxp1
           Xright = Aright
           Iright = Np1
           If ( J .eq. Nfxp1 )                      Go To 60
           Xright = Fixpnt(J)
C
C....     Determine Where The J-Th Fixed Point Should Fall In The
C....     New Mesh - This Is Xi(Iright) And The (J-1)St Fixed
C....     Point Is In Xi(Ileft)
C
           Nmin = Int((Xright-Aleft) / (Aright-Aleft) * Dble(N) + 1.5d0)
           If (Nmin .gt. N-Nfxpnt+J)  Nmin = N - Nfxpnt + J
           Iright = Max (Ileft+1, Nmin)
   60      Xi(Iright) = Xright
C
C....     Generate Equally Spaced Points Between The J-1st And The
C....     J-Th Fixed Points.
C
           Nregn = Iright - Ileft - 1
           If ( Nregn .eq. 0 )                      Go To 80
           Dx = (Xright - Xleft) / Dfloat(Nregn+1)
           Do 70 I = 1, Nregn
   70      Xi(Ileft+I) = Xleft  +  Dfloat(I) * Dx
   80      Ileft = Iright
           Xleft = Xright
   90 Continue
      Go To 320
C
C.... Mode = 2  Halve The Current Mesh (I.e. Double Its Size)
C
  100 N2 = 2 * N
C
C.... Check That N Does Not Exceed Storage Limitations
C
 112  If ( N2 .le. Nmax )                           Go To 120
C
C.... If Possible, Try With N = Nmax. Redistribute First.
C
      If ( Mode .eq. 2 )                            Go To 110
      N = Nmax / 2
      Go To 220
  110 If ( Iprint .lt. 1 )  Write(Iout,1001)
      N = N2
      Return
C
C.... Calculate The Old Approximate Solution Values At
C.... Points To Be Used In  Errchk  For Error Estimates.
C.... If  Mshflg   = 1 An Error Estimate Was Obtained For
C.... For The Old Approximation So Half The Needed Values
C.... Will Already Be In  Valstr .
C
  120 If ( Mshflg .eq. 0 )                          Go To 140
C
C.... Save In  Valstr  The Values Of The Old Solution
C.... At The Relative Positions 1/6 And 5/6 In Each Subinterval.
C
      Kstore = 1
      Do 130 I = 1, Nold
          Hd6 = (Xiold(I+1) - Xiold(I)) / 6.d0
          X = Xiold(I) + Hd6
          Call Approx (I, X, Valstr(Kstore), Asave(1,1), Dummy, Xiold,
     +         Nold, Z, Dmz, K, Ncomp, Mmax, M, Mstar, 4, Dummy, 0)
          X = X + 4.d0 * Hd6
          Kstore = Kstore  +  3 * Mstar
          Call Approx (I, X, Valstr(Kstore), Asave(1,4), Dummy, Xiold,
     +         Nold, Z, Dmz, K, Ncomp, Mmax, M, Mstar, 4, Dummy, 0)
          Kstore = Kstore  +  Mstar
  130 Continue
      Go To 160
C
C.... Save In  Valstr  The Values Of The Old Solution
C.... At The Relative Positions 1/6, 2/6, 4/6 And 5/6 In
C.... Each Subinterval.
C
  140 Kstore = 1
      Do 150 I = 1, N
         X = Xi(I)
         Hd6 = (Xi(I+1) - Xi(I)) / 6.d0
         Do 150 J = 1, 4
           X = X + Hd6
           If ( J .eq. 3 )  X = X + Hd6
           Call Approx (I, X, Valstr(Kstore), Asave(1,J), Dummy, Xiold,
     +          Nold, Z, Dmz, K, Ncomp, Mmax, M, Mstar, 4, Dummy, 0)
           Kstore = Kstore  +  Mstar
  150 Continue
  160 Mshflg = 0
C
C.... Generate The Halved Mesh.
C
      J = 2
      Do 170 I = 1, N
           Xi(J) = (Xiold(I) + Xiold(I+1)) / 2.d0
           Xi(J+1) = Xiold(I+1)
  170 J = J + 2
      N = N2
      Go To 320
C
C.... Mode = 1  We Do Mesh Selection If It Is Deemed Worthwhile
C
  180 If ( Nold .eq. 1 )                            Go To 100
      If ( Nold .le. 2*Nfxpnt )                     Go To 100
C
C.... The First Interval Has To Be Treated Separately From The
C.... Other Intervals (Generally The Solution On The (I-1)St And Ith
C.... Intervals Will Be Used To Approximate The Needed Derivative, But
C.... Here The 1st And Second Intervals Are Used.)
C
      I = 1
      Hiold = Xiold(2) - Xiold(1)
      Call Horder (1, D1, Hiold, Dmz, Ncomp, K)
      Hiold = Xiold(3) - Xiold(2)
      Call Horder (2, D2, Hiold, Dmz, Ncomp, K)
      Accum(1) = 0.d0
      Slope(1) = 0.d0
      Oneovh = 2.d0 / ( Xiold(3) - Xiold(1) )
      Do 190 J = 1, Ntol
         Jj = Jtol(J)
         Jz = Ltol(J)
         Jz1 = Jz + Mstar
         Zap = Min(Abs(Z(Jz)),Abs(Z(Jz1)))
         D3(Jj) = D2(Jj)
         Slope(1) = Max(Slope(1),(Abs(D2(Jj)-D1(Jj))*Wgtmsh(J)*
     +        Oneovh/Max(Zap,1.d0)) **Root(J))
 190  Continue
      Iflip = 1
C
C.... Go Through The Remaining Intervals Generating  Slope
C.... And  Accum.
C
      X2 = Xiold(3)
      Hiold = Xiold(2) - Xiold(1)
      Hiold1 = X2 - Xiold(2)
      Do 210 I = 2, Nold-1
         X1 = X2
         X2 = Xiold(I+2)
         Hioldm1 = Hiold
         Hiold = Hiold1
         Hiold1 = X2 - X1
         Do 195 J = 1,Ntol
            Jj = Jtol(J)
            If ( Iflip .eq. -1 ) D1(Jj) = D3(Jj)
            If ( Iflip .eq. 1 )  D2(Jj) = D3(Jj)
 195     Continue
         Call Horder ( I+1, D3, Hiold1, Dmz, Ncomp, K)
         Rat1 = Max(Hiold/Hioldm1,Hioldm1/Hiold)
         Rat2 = Max(Hiold/Hiold1,Hiold1/Hiold)
         If (Rat1 .le. 2.d0 .and. Rat2 .le. 2.d0) Then
            Lsrs = -1
            Oneovh = 2.d0 / ( Hiold + Hioldm1 )
            Oneovh1 = 2.d0 / ( Hiold1 + Hiold )   
         Else If (Rat1 .lt. Rat2) Then
            Lsrs = 0
            Oneovh = 2.d0 / ( Hiold + Hioldm1 )
         Else
            Lsrs = 1
            Oneovh = 2.d0 / ( Hiold1 + Hiold )
         Endif
         Slope(I) = 0.d0
         
C.... Evaluate Function To Be Equidistributed

         Do 200 J = 1, Ntol
            Jj = Jtol(J)
            Jz = Ltol(J) + (I-1)*Mstar
            Jz1 = Ltol(J) + I*Mstar
            Zap = Min(Abs(Z(Jz)),Abs(Z(Jz1)))
            Zap = Max(Zap,1.d0)
            If (Lsrs .eq. 0) Then
               Slope(I) = Max(Slope(I),(Abs(D2(Jj)-D1(Jj))*Wgtmsh(J)*
     +              Oneovh/Zap) **Root(J))
            Else If ((Lsrs .eq. 1) .and. (Iflip .eq. 1)) Then
               Slope(I) = Max(Slope(I),(Abs(D2(Jj)-D3(Jj))*Wgtmsh(J)*
     +              Oneovh/Zap) **Root(J))
            Else If ((Lsrs .eq. 1) .and. (Iflip .eq. -1)) Then
               Slope(I) = Max(Slope(I),(Abs(D3(Jj)-D1(Jj))*Wgtmsh(J)*
     +              Oneovh/Zap) **Root(J))
            Else If ((Lsrs .eq. -1) .and. (Iflip .eq. 1)) Then
               Slope(I) = Max(Slope(I),(Abs(D2(Jj)-D1(Jj))*Wgtmsh(J)*
     +              Oneovh/Zap) **Root(J))
               Slope(I) = Max(Slope(I),(Abs(D2(Jj)-D3(Jj))*Wgtmsh(J)*
     +              Oneovh1/Zap) **Root(J))
            Else If ((Lsrs .eq. -1) .and. (Iflip .eq. -1)) Then
               Slope(I) = Max(Slope(I),(Abs(D2(Jj)-D1(Jj))*Wgtmsh(J)*
     +              Oneovh/Zap) **Root(J))
               Slope(I) = Max(Slope(I),(Abs(D3(Jj)-D1(Jj))*Wgtmsh(J)*
     +              Oneovh1/Zap) **Root(J))
            Endif
 200     Continue
         Iflip = -Iflip
 210  Continue

C.... Similarly To First Subinterval, Treat Last Subinterval Separately,

      I = Nold
      Slope(I) = 0.d0 
      Oneovh = 2.d0/(X2-Xiold(I-1))
      Do 215 J = 1, Ntol 
         Jj = Jtol(J)
         Jz = Ltol(J) + (I-1)*Mstar
         Jz1 = Ltol(J) + I * Mstar
         Zap = Min(Abs(Z(Jz)),Abs(Z(Jz1)))
         Zap = Max(Zap,1.d0)
         If (Iflip .eq. -1) Then
            Slope(I) = Max(Slope(I),(Abs(D2(Jj)-D3(Jj))*Wgtmsh(J)*
     +           Oneovh/Max(Zap,1.d0)) **Root(J))
         Else
            Slope(I) = Max(Slope(I),(Abs(D3(Jj)-D1(Jj))*Wgtmsh(J)*
     +           Oneovh/Max(Zap,1.d0)) **Root(J))
         Endif
 215  Continue
C
C.... Accumulate Approximate Integral Of Function To Be Equidistributed
C 
      Accum(1) = 0.0d0
      Slphmx = 0.0d0
      Philrg = 0.d0
      Hordlrg = 0.d0
      Do 225 J = 1,Nold
         Temp = Slope(J)*(Xiold(J+1)-Xiold(J)) 
         Slphmx = Max(Slphmx,Temp)
         Accum(J+1) = Accum(J) + Temp
         If (Iatt .eq. -1) Then
            Slpold(J) = Slope(J)
            Voldmsh(J) = Xiold(J)
         Else If (Iatt .eq. 0 .and. Slope(J) .ge. Philrg) Then
            Philrg = Slope(J)
            Hordlrg = Max(Hordlrg,Temp)
            Imreg = J
         Endif
 225  Continue

      If (Iatt .eq. -1) Then
         Voldmsh(Nold+1) = Xiold(Nold+1)
         Nvold = Nold
      Endif

      If (Iatt .eq. 0) Then
         Pmax(2) = Philrg
         Hord(2) = Hordlrg
         Xloc1 = Xiold(Imreg)
         Xloc2 = Xiold(Imreg+1)
         Ichkpt = Nvold/2
         If (Xloc1 .lt. Voldmsh(Ichkpt)) Ichkpt = 1
         Xb = Voldmsh(Ichkpt)
         Do 227 J = Ichkpt, Nvold
            Xa = Xb
            Xb = Voldmsh(J+1)
            If (Xloc1 .ge. Xa .and. Xloc1 .lt. Xb) Then
               If (Xloc2-Xb .lt. (Xloc2-Xloc1)/2.d0) Then
                  Pmax(1) = Slpold(J)
                  Hord(1) = (Xb-Xa)*Pmax(1)
               Else
                  Pmax(1) = Slpold(J+1)
                  Hord(1) = (Voldmsh(J+2)-Xb)*Pmax(1)
               Endif
               Goto 228
            Endif
 227     Continue
 228  Endif

      Avrg = Accum(Nold+1)/Dfloat(Nold)
      Degequ = Avrg/Max(Slphmx,Precis)
C
C.... Naccum = Expected N To Achieve .1x User Requested Tolerances
C
      Naccum = Int(Accum(Nold+1) + 1.d0)
      If ( Iprint .lt. 0 )  Write(Iout,1002) Degequ, Naccum
C
C.... This Assures That Halving Will Be Possible Later (For Error Est)
C
      Nmax2 = Nmax / 2
C
C.... We Do Not Know In Advance Exactly How Many Points Will Be In 
C.... The New Mesh. The Parameter Nsafety Has Been Introduced So That 
C.... We Do Not Violate The Maximum Mesh Size.
C
      If (Ifinal .eq. 1 .and. Iatt .eq. 0) Then
         If ( Iprint .lt. 1 )  Write(Iout,1007) 
         N = Min(Nmax2-10,Naccum/2)
      Else
         N = Min(Nmax-20,Naccum)
      Endif
  220 Noldp1 = Nold + 1
      Nsafety = N
      If ( N .lt. Nfxp1 )  N = Nfxp1
      Mshflg = 0
C
C.... Having Decided To Generate A New Mesh With N Subintervals We Now
C.... Do So, Taking Into Account That The Nfxpnt Points In The Array
C.... Fixpnt Must Be Included In The New Mesh.
C
      In = 1
      Accl = 0.d0
      Lold = 2
      Xi(1) = Aleft
      Xi(N+1) = Aright
      Do 310 I = 1, Nfxp1
           If ( I .eq. Nfxp1 )                      Go To 250
           Do 230 J = Lold, Noldp1
             Lnew = J
             If ( Fixpnt(I) .le. Xiold(J) )         Go To 240
  230      Continue
  240      Continue
           Accr = Accum(Lnew) + (Fixpnt(I)-Xiold(Lnew))*Slope(Lnew-1)
           Nregn = Int((Accr-Accl) / Accum(Noldp1) * Dble(N) - 0.5d0)
           Nregn = Min(Nregn, N - In - Nfxp1 + I)
           Xi(In + Nregn + 1) = Fixpnt(I)
           Go To 260
  250      Accr = Accum(Noldp1)
           Lnew = Noldp1
           Nregn = N - In
  260      If ( Nregn .eq. 0 )                      Go To 300
           Temp = Accl
           Tsum = (Accr - Accl) / Dfloat(Nregn+1)
           Do 290 J = 1, Nregn
             In = In + 1
             Temp = Temp + Tsum
             Do 270 L = Lold, Lnew
               Lcarry = L
               If ( Temp .le. Accum(L) )            Go To 280
  270        Continue
  280        Continue
             Lold = Lcarry
  290      Xi(In) = Xiold(Lold-1) + (Temp - Accum(Lold-1)) /
     +     Slope(Lold-1)
  300      In = In + 1
           Accl = Accr
           Lold = Lnew
  310 Continue
      If (N .eq. Nmax2) Goto 320
      If (N .gt. Nmax2-30) Goto 321
C
C.... Check That This New Mesh Will Not Create Any (Innacurately) Large
C.... Subinterval Errors. If It Does Then Extra Points Must Be Added To
C.... The Mesh.
C
      I = 1
      Jpt1 = 1
      X2 = Xi(1)
      Slprgt = Slope(1)
      Equcon = Accum(Noldp1)/Dble(N)
      Equmax = 1.5d0*Equcon
      If (Iatt .eq. 0) Hsml = Equcon/Slope(Imreg)
C
 331  X1 = X2
      Jpt = Jpt1
      Slplft = Slprgt
 335  X2 = Xi(I+1)
      Xrpt = Xiold(Jpt)
      Do 332 J = Jpt,Nold
        Xlpt = Xrpt
        Xrpt = Xiold(J+1)
        If (X2 .ge. Xlpt .and. X2 .le. Xrpt) Then
          Slprgt = Slope(J)
          Jpt1 = J
          Goto 333
        Endif
 332  Continue

 333  If (Slplft .gt. Slprgt) Then
        Iend = 1
        Appint = Slplft*(X2-X1)
      Else
        Iend = 2
        Appint = Slprgt*(X2-X1)
      Endif

      If (Appint .gt. Equmax) Then
        N = N+1
        Do 334 Jj = N,I+1,-1
          Xi(Jj+1) = Xi(Jj)
 334    Continue
        Requ = Max(2.d0,Appint/Equcon)
        Equaim = Appint/Requ
        If (Iend .eq. 1) Xi(I+1) = X1+Equaim/Slplft
        If (Iend .eq. 2) Xi(I+1) = X2-Equaim/Slprgt
        If (X2 .eq. Xi(I+1) .or. X1 .eq. Xi(I+1)) Then 
           Iprec = 2
           Goto 320
        Endif
        Goto 335
      Else
        I = I+1
        If (I .lt. N+1) Goto 331
      Endif
 321  Continue
C
C.... Check Two Subintervals Approximately Same Size
C
      If (N .eq. 1) Then
        Xpt = (Xi(1)+Xi(2))/2.d0
        Xi(3) = Xi(2)
        Xi(2) = Xpt
        N = N+1
      Endif 
                
      I = 1
      X1 = Xi(I)
      X2 = Xi(I+1)
      X3 = Xi(I+2)
      Hiold = X2 - X1
      Hiold1 = X3 - X2
      If (X1 .eq. X2 .or. X2 .eq. X3) Then
         Iprec = 2
         Goto 320
      Endif
      Rat2 = Max(Hiold1/Hiold,Hiold/Hiold1)
      If (Rat2 .gt. 2.d0) Then
        Xpt = (X1+X2)/2.d0
        N = N+1
        Do 400 Jj = N,2,-1
          Xi(Jj+1) = Xi(Jj)
 400    Continue
        Xi(2) = Xpt
        Hiold = Hiold/2.d0
        I = I+1
      Endif
      If (Rat2 .le. 2.d0 .and. N .eq. 2) Goto 320
 410  I = I+1
      X1 = X2
      X2 = X3
      X3 = Xi(I+2)
      If (X2 .eq. X3) Then
         Iprec = 2
         Goto 320
      Endif
      Hioldm1 = Hiold
      Hiold = Hiold1
      Hiold1 = X3 - X2
      Rat1 = Max(Hiold/Hioldm1,Hioldm1/Hiold)
      Rat2 = Max(Hiold/Hiold1,Hiold1/Hiold)
      If ((Rat1 .gt. 2.d0) .and. (Rat2 .gt. 2.d0)) Then
        Xpt = (X1+X2)/2.d0
        N = N+1
        Do 420 Jj = N,I+1,-1
          Xi(Jj+1) = Xi(Jj)
 420    Continue
        Xi(I+1) = Xpt
        Hiold = Hiold/2.d0
        I = I+1
      Endif
      If (N .eq. Nmax2) Goto 320
      If (I .lt. N-1) Goto 410
      Hioldm1 = Hiold
      Hiold = Hiold1
      Rat1 = Max(Hioldm1/Hiold,Hiold/Hioldm1)
      If (Rat1 .gt. 2.d0) Then
        Xpt = (X2+X3)/2.d0
        N = N+1
        Xi(N+1) = X3
        Xi(N) = Xpt
      Endif
C
      If (N .gt. Nmax) Then
        N = Nsafety - 20
        Goto 220
      Endif

  320 Continue
      Np1 = N + 1
      If (Mode .le. 2) Then
         If ( Iprint .eq. -1 )  Write(Iout,1003) N, (Xi(I),I = 1,Np1)
         If ( Iprint .eq. 0 )  Write(Iout,1004) N
      Else
         If ( Iprint .eq. -1 )  Write(Iout,1005) N, (Xi(I),I = 1,Np1)
         If ( Iprint .eq. 0 )  Write(Iout,1006) N
      Endif

      Nz   = Mstar * (N + 1)
      Ndmz = Kd * N
     
      Return
C----------------------------------------------------------------
 1000 Format(/20h The Former Mesh (Of,I5,15h Subintervals),,
     +       600(/8f12.6))
 1001 Format (/23h  Expected N Too Large  )
 1002 Format(/21h Mesh Selection Info,/30h Degree Of Equidistribution = 
     +     , F8.5, 28h Prediction For Required N = , I8)
 1003 Format(/17h The New Mesh (Of,I5,16h Subintervals), ,600(/8f12.6))
 1004 Format(17h The New Mesh Has,I5,13h Subintervals)
 1005 Format(/21h The Initial Mesh (Of,I5,16h Subintervals), ,
     +     600(/8f12.6))
 1006 Format(/21h The Initial Mesh Has,I5,13h Subintervals)
 1007 Format(' Halving Mesh (And Then Doubling) In Order To',
     +     ' Calculate An Error Estimate')
      End

      Subroutine Consts (K, Rho, Coef)
C
C**********************************************************************
C
C   Purpose
C            Assign (Once) Values To Various Array Constants.
C
C   Arrays Assigned During Compilation:
C     Cnsts1 - Weights For Extrapolation Error Estimate
C     Cnsts2 - Weights For Mesh Selection
C              (The Above Weights Come From The Theoretical Form For
C              The Collocation Error -- See [3])
C
C   Arrays Assigned During Execution:
C     Wgterr - The Particular Values Of Cnsts1 Used For Current Run
C              (Depending On K, M)
C     Wgtmsh - Gotten From The Values Of Cnsts2 Which In Turn Are
C              The Constants In The Theoretical Expression For The
C              Errors. The Quantities In Wgtmsh Are 10x The Values
C              In Cnsts2 So That The Mesh Selection Algorithm
C              Is Aiming For Errors .1x As Large As The User
C              Requested Tolerances.
C     Jtol   - Components Of Differential System To Which Tolerances
C              Refer (Viz, If Ltol(I) Refers To A Derivative Of U(J),
C              Then Jtol(I) = J)
C     Root   - Reciprocals Of Expected Rates Of Convergence Of Compo-
C              Nents Of Z(J) For Which Tolerances Are Specified
C     Rho    - The K Collocation Points On (0,1)
C     Coef   -
C     Acol  -  The Runge-Kutta Coefficients Values At Collocation
C              Points
C
C**********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Rho(7), Coef(K,*), Cnsts1(28), Cnsts2(28), Dummy(1)
C
      Common /Colord/ Kdum, Ncomp, Mstar, Kd, Mmax, M(20)
      Common /Colbas/ B(28), Acol(28,7), Asave(28,4)
      Common /Colest/ Wgtmsh(40), Wgterr(40), Tolin(40),
     +                Root(40), Jtol(40), Ltol(40), Ntol
C
      Data Cnsts1 /    .25d0,     .625d-1,  7.2169d-2, 1.8342d-2,
     +     1.9065d-2, 5.8190d-2, 5.4658d-3, 5.3370d-3, 1.8890d-2,
     +     2.7792d-2, 1.6095d-3, 1.4964d-3, 7.5938d-3, 5.7573d-3,
     +     1.8342d-2, 4.673d-3,  4.150d-4,  1.919d-3,  1.468d-3,
     +     6.371d-3,  4.610d-3,  1.342d-4,  1.138d-4,  4.889d-4,
     +     4.177d-4,  1.374d-3,  1.654d-3,  2.863d-3  /
      Data Cnsts2 /   1.25d-1,   2.604d-3,  8.019d-3,  2.170d-5,
     +     7.453d-5,  5.208d-4,  9.689d-8,  3.689d-7,  3.100d-6,
     +     2.451d-5,  2.691d-10, 1.120d-9,  1.076d-8,  9.405d-8,
     +     1.033d-6,  5.097d-13, 2.290d-12, 2.446d-11, 2.331d-10,
     +     2.936d-9,  3.593d-8,  7.001d-16, 3.363d-15, 3.921d-14,
     +     4.028d-13, 5.646d-12, 7.531d-11, 1.129d-9  /
C
C.... Assign Weights For Error Estimate
C
      Koff = K * ( K + 1 ) / 2
      Iz = 1
      Do 10 J = 1, Ncomp
           Mj = M(J)
           Do 10 L = 1, Mj
             Wgterr(Iz) = Cnsts1(Koff - Mj + L)
             Iz = Iz + 1
   10 Continue
C
C.... Assign Array Values For Mesh Selection: Wgtmsh, Jtol, And Root
C
      Jcomp = 1
      Mtot = M(1)
      Do 40 I = 1, Ntol
           Ltoli = Ltol(I)
   20      Continue
           If ( Ltoli .le. Mtot )                   Go To 30
           Jcomp = Jcomp + 1
           Mtot = Mtot + M(Jcomp)
           Go To 20
   30      Continue
           Jtol(I) = Jcomp
           Wgtmsh(I) = 1.d1 * Cnsts2(Koff+Ltoli-Mtot) / Tolin(I)
           Root(I) = 1.d0 / Dfloat(K+Mtot-Ltoli+1)
   40 Continue
C
C.... Specify Collocation Points
C
      Go To (50,60,70,80,90,100,110), K
   50 Rho(1) = 0.d0
      Go To 120
   60 Rho(2) = .57735026918962576451d0
      Rho(1) = - Rho(2)
      Go To 120
   70 Rho(3) = .77459666924148337704d0
      Rho(2) = .0d0
      Rho(1) = - Rho(3)
      Go To 120
   80 Rho(4) = .86113631159405257523d0
      Rho(3) = .33998104358485626480d0
      Rho(2) = - Rho(3)
      Rho(1) = - Rho(4)
      Go To 120
   90 Rho(5) = .90617984593866399280d0
      Rho(4) = .53846931010568309104d0
      Rho(3) = .0d0
      Rho(2) = - Rho(4)
      Rho(1) = - Rho(5)
      Go To 120
  100 Rho(6) = .93246951420315202781d0
      Rho(5) = .66120938646626451366d0
      Rho(4) = .23861918608319690863d0
      Rho(3) = -Rho(4)
      Rho(2) = -Rho(5)
      Rho(1) = -Rho(6)
      Go To 120
  110 Rho(7) = .949107991234275852452d0
      Rho(6) = .74153118559939443986d0
      Rho(5) = .40584515137739716690d0
      Rho(4) = 0.d0
      Rho(3) = -Rho(5)
      Rho(2) = -Rho(6)
      Rho(1) = -Rho(7)
  120 Continue
C
C.... Map (-1,1) To (0,1) By  T = .5 * (1. + X)
C
      Do 130 J = 1, K
         Rho(J) = .5d0 * (1.d0 + Rho(J))
  130 Continue
C
C.... Now Find Runge-Kutta Coeffitients B, Acol And Asave
C.... The Values Of Asave Are To Be Used In  Newmsh  And Errchk .
C
      Do 140 J = 1, K
         Do 135 I = 1, K
  135      Coef(I,J) = 0.d0
         Coef(J,J) = 1.d0
         Call Vmonde (Rho, Coef(1,J), K)
  140 Continue
      Call Rkbas ( 1.d0, Coef, K, Mmax, B, Dummy, 0)
      Do 150 I = 1, K
         Call Rkbas ( Rho(I), Coef, K, Mmax, Acol(1,I), Dummy, 0)
  150 Continue
      Call Rkbas ( 1.d0/6.d0, Coef, K, Mmax, Asave(1,1), Dummy, 0)
      Call Rkbas ( 1.d0/3.d0, Coef, K, Mmax, Asave(1,2), Dummy, 0)
      Call Rkbas ( 2.d0/3.d0, Coef, K, Mmax, Asave(1,3), Dummy, 0)
      Call Rkbas ( 5.d0/6.d0, Coef, K, Mmax, Asave(1,4), Dummy, 0)
      Return
      End
      Subroutine Errchk (Xi, Z, Dmz, Valstr, Ifin)
C
C**********************************************************************
C
C      Purpose
C               Determine The Error Estimates And Test To See If The
C               Error Tolerances Are Satisfied.
C
C      Variables
C        Xi     - Current Mesh Points
C        Valstr - Values Of The Previous Solution Which Are Needed
C                 For The Extrapolation- Like Error Estimate.
C        Wgterr - Weights Used In The Extrapolation-Like Error
C                 Estimate. The Array Values Are Assigned In
C                 Subroutine  Consts.
C        Errest - Storage For Error Estimates
C        Err    - Temporary Storage Used For Error Estimates
C        Z      - Approximate Solution On Mesh Xi
C        Ifin   - A 0-1 Variable. On Return It Indicates Whether
C                 The Error Tolerances Were Satisfied
C        Mshflg - Is Set By Errchk To Indicate To Newmsh Whether
C                 Any Values Of The Current Solution Are Stored In
C                 The Array Valstr. (0 For No, 1 For Yes)
C
C**********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Err(40), Errest(40), Dummy(1)
      Dimension Xi(*), Z(*), Dmz(*), Valstr(*)
C
      Common /Colout/ Precis, Iout, Iprint
      Common /Colord/ K, Ncomp, Mstar, Kd, Mmax, M(20)
      Common /Colapr/ N, Nold, Nmax, Nz, Ndmz, Mshflg
      Common /Colbas/ B(28), Acol(28,7), Asave(28,4)
      Common /Colest/ Wgtmsh(40), Wgterr(40), Tolin(40),
     +                Root(40), Jtol(40), Ltol(40), Ntol
C
C.... Error Estimates Are To Be Generated And Tested
C.... To See If The Tolerance Requirements Are Satisfied.
C
      Ifin = 1
      Mshflg = 1
      Do 10 J = 1, Mstar
   10   Errest(J) = 0.d0
      Do 60 Iback = 1, N
           I = N + 1 - Iback
C
C....     The Error Estimates Are Obtained By Combining Values Of
C....     The Numerical Solutions For Two Meshes.
C....     For Each Value Of Iback We Will Consider The Two
C....     Approximations At 2 Points In Each Of
C....     The New Subintervals.  We Work Backwards Through
C....     The Subinterval So That New Values Can Be Stored
C....     In Valstr In Case They Prove To Be Needed Later
C....     For An Error Estimate. The Routine  Newmsh
C....     Filled In The Needed Values Of The Old Solution
C....     In Valstr.
C
           Knew = ( 4 * (I-1) + 2 ) * Mstar + 1
           Kstore = ( 2 * (I-1) + 1 ) * Mstar + 1
           X = Xi(I) +  (Xi(I+1)-Xi(I)) * 2.d0 / 3.d0
           Call Approx (I, X, Valstr(Knew), Asave(1,3), Dummy, Xi,
     +            N, Z, Dmz, K, Ncomp, Mmax, M, Mstar, 4, Dummy, 0)
           Do 20 L = 1,Mstar
             Err(L) = Wgterr(L) * Abs(Valstr(Knew) -
     +       Valstr(Kstore))
             Knew = Knew + 1
             Kstore = Kstore + 1
   20      Continue
           Knew = ( 4 * (I-1) + 1 ) * Mstar + 1
           Kstore = 2 * (I-1) * Mstar + 1
           X = Xi(I) +  (Xi(I+1)-Xi(I)) / 3.d0
           Call Approx (I, X, Valstr(Knew), Asave(1,2), Dummy, Xi,
     +            N, Z, Dmz, K, Ncomp, Mmax, M, Mstar, 4, Dummy, 0)
           Do 30 L = 1,Mstar
             Err(L) = Err(L) + Wgterr(L) * Abs(Valstr(Knew) -
     +       Valstr(Kstore))
             Knew = Knew + 1
             Kstore = Kstore + 1
   30      Continue
C
C....     Find Component-Wise Maximum Error Estimate
C
           Do 40 L = 1,Mstar
             Lzpt = L + (I-1) * Mstar
             Errest(L) = Max(Errest(L),Err(L)/Max(Abs(Z(Lzpt)),1.d0))
   40      Continue
C


C....     Test Whether The Tolerance Requirements Are Satisfied
C....     In The I-Th Interval.
C
           If ( Ifin .eq. 0 )                       Go To 60
           Do 50 J = 1, Ntol
             Ltolj = Ltol(J)
             Ltjz = Ltolj  +  (I-1) * Mstar
           If ( Err(Ltolj)  .gt. 
     +          Tolin(J) * Max(Abs(Z(Ltjz)),1.d0) )  Ifin = 0
   50      Continue
   60 Continue
      If ( Iprint .gt. 0 )                          Return
      If ( Ifin .eq. 0 .and. Iprint .eq. 0) Return
      Write(Iout,1000)
      Lj = 1
      Do 70 J = 1,Ncomp
           Mj = Lj - 1 + M(J)
           Write(Iout,1001) J, (Errest(L), L =  Lj, Mj)
           Lj = Mj + 1
   70 Continue
      Return
C---------------------------------------------------------------------
 1000 Format (/52h The Estimated (Mixed Relative/Absolute) Errors Are,)
 1001 Format (3h U(, I2, 3h) -,4d12.4)
      End
C---------------------------------------------------------------------
C                            P A R T  3
C          Collocation System Setup Routines
C---------------------------------------------------------------------
C
      Subroutine Lsyslv (Msing, Xi, Xiold, Z, Dmz, Delz, Deldmz,
     +           G, W, V, Rhs, Dmzo, Integs, Ipvtg, Ipvtw, Rnorm,
     +           Mode, Fsub, Dfsub, Gsub, Dgsub, Guess, Eps )
C*********************************************************************
C
C   Purpose
C         This Routine Controls The Set Up And Solution Of A Linear
C      System Of Collocation Equations.
C         The Matrix  G  Is Cast Into An Almost Block Diagonal
C      Form By An Appropriate Ordering Of The Columns And Solved
C      Using The Package Of De Boor-Weiss [5]. The Matrix Is Composed
C      Of N Blocks. The I-Th Block Has The Size
C                  Integs(1,I) * Integs(2,I).
C      It Contains In Its Last Rows The Linearized Collocation
C      Equations, Condensed As Described In [2],
C      And The Linearized Side Conditions Corresponding To
C      The I-Th Subinterval.  Integs(3,I)  Steps Of Gaussian
C      Elimination Are Applied To It To Achieve A  Partial Plu
C      Decomposition.  The Right Hand Side Vector Is Put Into  Rhs
C      And The Solution Vector Is Returned In  Delz And Deldmz.
C
C         Lsyslv Operates According To One Of 5 Modes:
C      Mode = 0 - Set Up The Collocation Matrices  V , W , G
C                 And The Right Hand Side  Rhs ,  And Solve.
C                 (For Linear Problems Only.)
C      Mode = 1 - Set Up The Collocation Matrices  V , W , G
C                 And The Right Hand Sides  Rhs  And  Dmzo ,
C                 And Solve. Also Set Up  Integs .
C                 (First Iteration Of Nonlinear Problems Only).
C      Mode = 2 - Set Up  Rhs  Only And Compute Its Norm.
C      Mode = 3 - Set Up  V, W, G  Only And Solve System.
C      Mode = 4 - Perform Forward And Backward Substitution Only
C                 (Do Not Set Up The Matrices Nor Form The Rhs).
C
C   Variables
C
C      Ig,Izeta  - Pointers To G,Zeta Respectively
C                       (Necessary To Keep Track Of Blocks Of G
C                       During Matrix Manipulations)
C      Idmz,Irhs,Iv,Iw - Pointers To  Rhs,V,W Rspectively
C      Df    - Partial Derivatives Of F From Dfsub
C      Rnorm - Euclidean Norm Of Rhs
C      Lside - Number Of Side Conditions In Current And Previous Blocks
C      Iguess = 1 When Current Soln Is User Specified Via  Guess
C             = 0 Otherwise
C
C*********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension  Z(*), Dmz(*), Delz(*), Deldmz(*), Xi(*), Xiold(*)
      Dimension  G(*), W(*), V(*),  Rhs(*) , Dmzo(*), Dummy(1)
      Dimension  Integs(3,*), Ipvtg(*), Ipvtw(*)
      Dimension  Zval(40), F(40), Dgz(40), Dmval(20), Df(800), At(28)
C
      Common /Colout/ Precis, Iout, Iprint
      Common /Colloc/ Rho(7), Coef(49)
      Common /Colord/ K, Ncomp, Mstar, Kd,  Mmax, M(20)
      Common /Colsid/ Zeta(40), Aleft, Aright, Izeta, Izsave
      Common /Colapr/ N, Nold, Nmax, Nz, Ndmz, Mshflg
      Common /Colnln/ Nonlin, Iter, Limit, Iguess
      Common /Colbas/ B(28), Acol(28,7), Asave(28,4)
C
      External Dfsub, Dgsub
C
      M1 = Mode + 1
      Go To (10, 30, 30, 30, 310), M1
C
C.... Linear Problem Initialization
C
   10 Do 20 I = 1,Mstar
   20 Zval(I) = 0.d0
C
C.... Initialization
C
   30 Idmz = 1
      Idmzo = 1
      Irhs = 1
      Ig = 1
      Iw = 1
      Iv = 1
      Izeta = 1
      Lside = 0
      Iold = 1
      Ncol = 2 * Mstar
      Rnorm = 0.d0
      If ( Mode .gt. 1 )                            Go To 80
C
C.... Build Integs (Describing Block Structure Of Matrix)
C
      Do 70 I = 1,N
           Integs(2,I) = Ncol
           If (I .lt. N)                            Go To 40
           Integs(3,N) = Ncol
           Lside = Mstar
           Go To 60
   40      Integs(3,I) = Mstar
   50      If( Lside .eq. Mstar )                   Go To 60
           If ( Zeta(Lside+1) .ge. Xi(I)+Precis )   Go To 60
           Lside = Lside + 1
           Go To 50
   60      Nrow = Mstar + Lside
   70 Integs(1,I) = Nrow
   80 Continue
      If ( Mode .eq. 2 )                            Go To 90
C
C.... Zero The Matrices To Be Computed
C
      Lw = Kd * Kd * N
      Do 84 L = 1, Lw
   84   W(L) = 0.d0
C
C.... The Do Loop 290 Sets Up The Linear System Of Equations.
C
  90  Continue
      Do 290 I = 1, N
C
C....     Construct A Block Of  A  And A Corresponding Piece Of  Rhs.
C
           Xii = Xi(I)
           H = Xi(I+1) - Xi(I)
           Nrow = Integs(1,I)
C
C....     Go Thru The Ncomp Collocation Equations And Side Conditions
C....     In The I-Th Subinterval
C
  100      If ( Izeta .gt. Mstar )                  Go To 140
           If ( Zeta(Izeta) .gt. Xii + Precis )      Go To 140
C
C....     Build Equation For A Side Condition.
C
           If ( Mode .eq. 0 )                       Go To 110
           If ( Iguess .ne. 1 )                     Go To 102
C
C....     Case Where User Provided Current Approximation
C
           Call Guess (Xii, Zval, Dmval, Eps)
           Go To 110
C
C....     Other Nonlinear Case
C
  102      If ( Mode .ne. 1 )                       Go To 106
           Call Approx (Iold, Xii, Zval, At, Coef, Xiold, Nold,
     +          Z, Dmz, K, Ncomp, Mmax, M, Mstar, 2, Dummy, 0)
           Go To 110
  106      Call Approx (I, Xii, Zval, At, Dummy, Xi, N, Z, Dmz,
     +                  K, Ncomp, Mmax, M, Mstar, 1, Dummy, 0)
  108      If ( Mode .eq. 3 )                       Go To 120
C
C....     Find  Rhs  Boundary Value.
C
  110      Call Gsub (Izeta, Zval, Gval, Eps)
           Rhs(Ndmz+Izeta) = -Gval
           Rnorm = Rnorm + Gval**2
           If ( Mode .eq. 2 )                       Go To 130
C
C....     Build A Row Of  A  Corresponding To A Boundary Point
C
  120      Call Gderiv (G(Ig), Nrow, Izeta, Zval, Dgz, 1, Dgsub, Eps)
  130      Izeta = Izeta + 1
           Go To 100
C
C....     Assemble Collocation Equations
C
  140      Do 220 J = 1, K
             Hrho = H * Rho(J)
             Xcol = Xii + Hrho
C
C....       This Value Corresponds To A Collocation (Interior)
C....       Point. Build The Corresponding  Ncomp  Equations.
C
             If ( Mode .eq. 0 )                     Go To 200
             If ( Iguess .ne. 1 )                   Go To 160
C
C....       Use Initial Approximation Provided By The User.
C
             Call Guess (Xcol, Zval, Dmzo(Irhs), Eps )
             Go To 170
C
C....       Find  Rhs  Values
C
  160        If ( Mode .ne. 1 )                     Go To 190
             Call Approx (Iold, Xcol, Zval, At, Coef, Xiold, Nold,
     +            Z, Dmz, K, Ncomp, Mmax, M, Mstar, 2, Dmzo(Irhs), 1)
C
  170        Call Fsub (Xcol, Zval, F, Eps)
             Do 180 Jj = 1, Ncomp
               Value = Dmzo(Irhs) - F(Jj)
               Rhs(Irhs) = - Value
               Rnorm = Rnorm + Value**2
               Irhs = Irhs + 1
  180        Continue
             Go To 210
C
C....       Evaluate Former Collocation Solution
C
  190        Call Approx (I, Xcol, Zval, Acol(1,J), Coef, Xi, N,
     +            Z, Dmz, K, Ncomp, Mmax, M, Mstar, 4, Dummy, 0)
             If ( Mode .eq. 3 )                     Go To 210
C
C....       Fill In  Rhs  Values (And Accumulate Its Norm).
C
             Call Fsub (Xcol, Zval, F, Eps)
             Do 195 Jj = 1, Ncomp
               Value = Dmz(Irhs) - F(Jj)
               Rhs(Irhs) = - Value
               Rnorm = Rnorm + Value**2
               Irhs = Irhs + 1
  195        Continue
             Go To 220
C
C....       The Linear Case
C
  200        Call Fsub (Xcol, Zval, Rhs(Irhs), Eps)
             Irhs = Irhs + Ncomp
C
C....       Fill In Ncomp Rows Of  W And V
C
  210        Call Vwblok (Xcol, Hrho, J, W(Iw), V(Iv), Ipvtw(Idmz), Kd,
     +       Zval, Df, Acol(1,J), Dmzo(Idmzo), Ncomp, Dfsub, Msing, Eps)
             If ( Msing .ne. 0 )                    Return
  220      Continue
C
C....     Build Global Bvp Matrix  G
C
           If ( Mode .ne. 2 )
     +      Call Gblock (H, G(Ig), Nrow, Izeta, W(Iw), V(Iv), Kd,
     +                  Dummy, Deldmz(Idmz), Ipvtw(Idmz), 1 )
           If ( I .lt. N )                          Go To 280
           Izsave = Izeta
  240      If ( Izeta .gt. Mstar )                  Go To 290
C
C....     Build Equation For A Side Condition.
C
           If ( Mode .eq. 0 )                       Go To 250
           If ( Iguess .ne. 1 )                     Go To 245
C
C....     Case Where User Provided Current Approximation
C
           Call Guess (Aright, Zval, Dmval, Eps)
           Go To 250
C
C....     Other Nonlinear Case
C
  245      If ( Mode .ne. 1 )                       Go To 246
           Call Approx (Nold+1, Aright, Zval, At, Coef, Xiold, Nold,
     +          Z, Dmz, K, Ncomp, Mmax, M, Mstar, 1, Dummy, 0)
           Go To 250
  246      Call Approx (N+1, Aright, Zval, At, Coef, Xi, N,
     +          Z, Dmz, K, Ncomp, Mmax, M, Mstar, 1, Dummy, 0)
  248      If ( Mode .eq. 3 )                       Go To 260
C
C....     Find  Rhs  Boundary Value.
C
 250       Call Gsub (Izeta, Zval, Gval, Eps)
           Rhs(Ndmz+Izeta) = - Gval
           Rnorm = Rnorm + Gval**2
           If ( Mode .eq. 2 )                       Go To 270
C
C....     Build A Row Of  A  Corresponding To A Boundary Point
C
 260       Izm = Izeta+Mstar
           Call Gderiv (G(Ig), Nrow, Izm, Zval, Dgz, 2, Dgsub, Eps)
 270       Izeta = Izeta + 1
           Go To 240
C
C....     Update Counters -- I-Th Block Completed
C
 280       Ig = Ig + Nrow * Ncol
           Iv = Iv + Kd * Mstar
           Iw = Iw + Kd * Kd
           Idmz = Idmz + Kd
           If ( Mode .eq. 1 )  Idmzo = Idmzo + Kd
  290 Continue
C
C....     Assembly Process Completed
C
      If ( Mode .eq. 0 .or. Mode .eq. 3 )           Go To 300
      Rnorm = Sqrt(Rnorm / Dfloat(Nz+Ndmz))
      If ( Mode .ne. 2 )                            Go To 300
      Return
C
C.... Solve The Linear System.
C
C.... Matrix Decomposition
C
  300 Call Fcblok (G, Integs, N, Ipvtg, Df, Msing)
C
C.... Check For Singular Matrix
C
      Msing = - Msing
      If( Msing .ne. 0 )                            Return
C
C.... Perform Forward And Backward Substitution For Mode = 4 Only.
C
  310 Continue
      Do 311 L = 1, Ndmz
        Deldmz(L) = Rhs(L)
  311 Continue
      Iz = 1
      Idmz = 1
      Iw = 1
      Izet = 1
      Do 320 I = 1, N
         Nrow = Integs(1,I)
         Izeta = Nrow + 1 - Mstar
         If ( I .eq. N ) Izeta = Izsave
  322    If ( Izet .eq. Izeta )                     Go To 324
           Delz(Iz-1+Izet) = Rhs(Ndmz+Izet)
           Izet = Izet + 1
         Go To 322
  324    H = Xi(I+1) - Xi(I)
         Call Gblock (H, G(1), Nrow, Izeta, W(Iw), V(1), Kd,
     +                Delz(Iz), Deldmz(Idmz), Ipvtw(Idmz), 2 )
         Iz = Iz + Mstar
         Idmz = Idmz + Kd
         Iw = Iw + Kd * Kd
         If ( I .lt. N )                            Go To 320
  326    If ( Izet .gt. Mstar )                     Go To 320
           Delz(Iz-1+Izet) = Rhs(Ndmz+Izet)
           Izet = Izet + 1
         Go To 326
  320 Continue
C
C.... Perform Forward And Backward Substitution For Mode = 0,2, Or 3.
C
      Call Sbblok (G, Integs, N, Ipvtg, Delz)
C
C.... Finaly Find Deldmz
C
      Call Dmzsol (Kd, Mstar, N, V, Delz, Deldmz)
C
      If ( Mode .ne. 1 )                            Return
      Do 321 L = 1, Ndmz
        Dmz(L) = Dmzo(L)
  321 Continue
      Iz = 1
      Idmz = 1
      Iw = 1
      Izet = 1
      Do 350 I = 1, N
         Nrow = Integs(1,I)
         Izeta = Nrow + 1 - Mstar
         If ( I .eq. N ) Izeta = Izsave
  330    If ( Izet .eq. Izeta )                     Go To 340
           Z(Iz-1+Izet) = Dgz(Izet)
           Izet = Izet + 1
         Go To 330
  340    H = Xi(I+1) - Xi(I)
         Call Gblock (H, G(1), Nrow, Izeta, W(Iw), Df, Kd,
     +                Z(Iz), Dmz(Idmz), Ipvtw(Idmz), 2 )
         Iz = Iz + Mstar
         Idmz = Idmz + Kd
         Iw = Iw + Kd * Kd
         If ( I .lt. N )                            Go To 350
  342    If ( Izet .gt. Mstar )                     Go To 350
            Z(Iz-1+Izet) = Dgz(Izet)
            Izet = Izet + 1
         Go To 342
  350 Continue
      Call Sbblok (G, Integs, N, Ipvtg, Z)
C
C.... Finaly Find Dmz
C
      Call Dmzsol (Kd, Mstar, N, V, Z, Dmz)
C
      Return
      End
      Subroutine Gderiv ( Gi, Nrow, Irow, Zval, Dgz, Mode, Dgsub, Eps)
C
C**********************************************************************
C
C   Purpose:
C
C      Construct A Collocation Matrix Row According To Mode:
C      Mode = 1  -  A Row Corresponding To A Initial Condition
C                   (I.e. At The Left End Of The Subinterval).
C      Mode = 2  -  A Row Corresponding To A Final Condition.
C
C   Variables:
C
C      Gi     - The Sub-Block Of The Global Bvp Matrix In
C               Which The Equations Are To Be Formed.
C      Nrow   - No. Of Rows In Gi.
C      Irow   - The Row In Gi To Be Used For Equations.
C      Zval   - Z(Xi)
C      Dg     - The Derivatives Of The Side Condition.
C
C**********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Gi(Nrow,*), Zval(*), Dgz(*), Dg(40)
C
      Common /Colord/ Kdum, Ndum, Mstar, Kd, Mmax, M(20)
      Common /Colsid/ Zeta(40), Aleft, Aright, Izeta, Idum
      Common /Colnln/ Nonlin, Iter, Limit, Iguess
C
C.... Zero Jacobian Dg
C
      Do 10 J = 1,Mstar
   10   Dg(J) = 0.d0
C
C.... Evaluate Jacobian Dg
C
      Call Dgsub (Izeta, Zval, Dg, Eps)
C
C.... Evaluate  Dgz = Dg * Zval  Once For A New Mesh
C
      If (Nonlin .eq. 0 .or. Iter .gt. 0)           Go To 30
      Dot = 0.d0
      Do 20 J = 1, Mstar
   20   Dot = Dot  +  Dg(J) * Zval(J)
      Dgz(Izeta) = Dot
C
C.... Branch According To  M O D E
C
   30 If ( Mode .eq. 2 )                            Go To 50
C
C.... Provide Coefficients Of The J-Th Linearized Side Condition.
C.... Specifically, At X = Zeta(J) The J-Th Side Condition Reads
C.... Dg(1)*Z(1) + ... +Dg(Mstar)*Z(Mstar) + G = 0
C
C
C.... Handle An Initial Condition
C
      Do 40 J = 1, Mstar
        Gi(Irow,J) =  Dg(J)
   40 Gi(Irow,Mstar+J) = 0.d0
      Return
C
C.... Handle A Final Condition
C
   50 Do 60 J =  1, Mstar
        Gi(Irow,J) = 0.d0
   60 Gi(Irow,Mstar+J) = Dg(J)
      Return
      End
      Subroutine Vwblok (Xcol, Hrho, Jj, Wi, Vi, Ipvtw, Kd, Zval,
     +                   Df, Acol, Dmzo, Ncomp, Dfsub, Msing, Eps)
C
C**********************************************************************
C
C   Purpose:
C
C      Construct A Group Of  Ncomp  Rows Of The Matrices  Wi  And  Vi.
C      Corresponding To An Interior Collocation Point.
C
C
C   Variables:
C
C      Xcol   - The Location Of The Collocation Point.
C      Jj     - Xcol Is The Jj-Th Of K Collocation Points
C               In The I-Th Subinterval.
C      Wi,Vi  - The I-Th Block Of The Collocation Matrix
C               Before Parameter Condensation.
C      Kd     - No. Of Rows In Vi And Wi .
C      Zval   - Z(Xcol)
C      Df     - The Jacobian At Xcol .
C      Jcomp  - Counter For The Component Being Dealt With.
C
C**********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Wi(Kd,*), Vi(Kd,*), Zval(*), Dmzo(*), Df(Ncomp,*)
      Dimension Ipvtw(*),  Ha(7,4), Acol(7,4), Basm(5)
C
      Common /Colord/ K, Ncdum, Mstar, Kdum, Mmax, M(20)
      Common /Colnln/ Nonlin, Iter, Limit, Iguess
C
C.... If Jj = 1 Initialize  Wi .
C
      If ( Jj .gt. 1 )                              Go To 30
      Do 10 Id = 1, Kd
        Wi(Id,Id) = 1.d0
   10 Continue
C
C.... Calculate Local Basis
C
   30        Fact = 1.d0
             Do 150 L = 1,Mmax
                Fact = Fact * Hrho / Dfloat(L)
                Basm(L) = Fact
                Do 150 J = 1,K
                   Ha(J,L) = Fact * Acol(J,L)
  150        Continue
C
C....zero Jacobian
C
      Do 40 Jcol = 1, Mstar
        Do 40 Ir = 1, Ncomp
   40 Df(Ir,Jcol) = 0.d0
C
C.... Build Ncomp Rows For Interior Collocation Point X.
C.... The Linear Expressions To Be Constructed Are:
C.... (M(Id))
C.... U     -  Df(Id,1)*Z(1) - ... - Df(Id,Mstar)*Z(Mstar)
C.... Id
C.... For Id = 1 To Ncomp.
C
      Call Dfsub (Xcol, Zval, Df, Ncomp, Eps)
      I0 = (Jj-1) * Ncomp
      I1 = I0 + 1
      I2 = I0 + Ncomp
C
C.... Evaluate  Dmzo = Dmz - Df * Zval  Once For A New Mesh
C
      If (Nonlin .eq. 0 .or. Iter .gt. 0)          Go To 60
      Do 50 J = 1, Mstar
        Fact = - Zval(J)
        Do 50 Id = 1, Ncomp
          Dmzo(I0+Id) = Dmzo(I0+Id)  +  Fact * Df(Id,J)
  50  Continue
C
C.... Loop Over The  Ncomp  Expressions To Be Set Up For The
C.... Current Collocation Point.
C
   60 Do 70 J = 1, Mstar
        Do 70 Id = 1, Ncomp
          Vi(I0+Id,J) = Df(Id,J)
   70 Continue
      Jn = 1
      Do 140 Jcomp = 1, Ncomp
         Mj = M(Jcomp)
         Jn = Jn + Mj
         Do 130 L = 1, Mj
            Jv = Jn - L
            Jw = Jcomp
            Do 90 J = 1, K
              Ajl = - Ha(J,L)
              Do 80 Iw = I1, I2
                 Wi(Iw,Jw) = Wi(Iw,Jw)  +  Ajl * Vi(Iw,Jv)
   80         Continue
   90       Jw = Jw + Ncomp
            Lp1 = L + 1
            If ( L .eq. Mj )                        Go To 130
            Do 110 Ll = Lp1, Mj
              Jdf = Jn - Ll
              Bl = Basm(Ll-L)
              Do 100 Iw = I1, I2
                Vi(Iw,Jv) = Vi(Iw,Jv)  +  Bl * Vi(Iw,Jdf)
  100         Continue
  110       Continue
  130    Continue
  140 Continue
      If ( Jj .lt. K )                          Return
C
C   ...decompose The Wi Block And Solve For The Mstar Columns Of Vi
C
C
C.... Do Parameter Condensation
C
      Msing = 0
      Call Dgefa  (Wi, Kd, Kd, Ipvtw, Msing)
C
C.... Check For Singularity
C
      If ( Msing .ne. 0 )                         Return
      Do 250 J =  1,Mstar
         Call Dgesl  (Wi, Kd, Kd, Ipvtw, Vi(1,J), 0)
  250 Continue
      Return
      End
      Subroutine Gblock (H, Gi, Nrow, Irow, Wi, Vi, Kd,
     +                   Rhsz, Rhsdmz, Ipvtw, Mode)
C
C**********************************************************************
C
C   Purpose:
C
C      Construct Collocation Matrix Rows According To Mode:
C      Mode = 1  -  A Group Of  Mstar    Rows Corresponding
C                   An Interior Mesh Interval.
C           = 2  -  Corresponding Right Hand Side
C
C   Variables:
C
C      H      - The  Local Stepsize.
C      Gi     - The Sub-Block Of The Collocation Matrix In
C               Which The Equations Are To Be Formed.
C      Wi     - The Sub-Block Of Noncondensed Collocation Equations,
C               Left-Hand Side Part.
C      Vi     - The Sub-Block Of Noncondensed Collocation Equations,
C               Right-Hand Side Part.
C      Rhsdmz - The Inhomogenous Term Of The Uncondensed Collocation
C               Equations.
C      Rhsz   - The Inhomogenous Term Of The Condensed Collocation
C               Equations.
C      Nrow   - No. Of Rows In Gi.
C      Irow   - The First Row In Gi To Be Used For Equations.
C
C**********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Hb(7,4), Basm(5)
      Dimension Gi(Nrow,*), Wi(*), Vi(Kd,*)
      Dimension Rhsz(*), Rhsdmz(*), Ipvtw(*)
C
      Common /Colord/  K, Ncomp, Mstar, Kdum, Mmax, M(20)
      Common /Colbas/ B(7,4), Acol(28,7), Asave(28,4)
C
C.... Compute Local Basis
C
      Fact = 1.d0
      Basm(1) = 1.d0
      Do 30 L = 1,Mmax
         Fact = Fact * H / Dfloat(L)
         Basm(L+1) = Fact
         Do 20 J = 1,K
   20       Hb(J,L) = Fact * B(J,L)
   30 Continue
C
C.... Branch According To  M O D E
C
      Go To (40, 110), Mode
C
C.... Set Right Gi-Block To Identity
C
   40 Continue
      Do 60 J = 1, Mstar
        Do 50 Ir = 1, Mstar
          Gi(Irow-1+Ir,J) = 0.d0
   50   Gi(Irow-1+Ir,Mstar+J) = 0.d0
   60 Gi(Irow-1+J,Mstar+J) = 1.d0
C
C.... Compute The Block Gi
C
      Ir = Irow
      Do 100 Icomp = 1, Ncomp
         Mj = M(Icomp)
         Ir = Ir + Mj
         Do 90 L = 1, Mj
            Id = Ir - L
            Do 80 Jcol = 1, Mstar
               Ind = Icomp
               Rsum = 0.d0
               Do 70 J = 1, K
                  Rsum = Rsum  -  Hb(J,L) * Vi(Ind,Jcol)
   70          Ind = Ind + Ncomp
               Gi(Id,Jcol) = Rsum
   80       Continue
            Jd = Id - Irow
            Do 85 Ll = 1, L
               Gi(Id,Jd+Ll) = Gi(Id,Jd+Ll) - Basm(Ll)
   85       Continue
   90    Continue
  100 Continue
      Return
C
C.... Compute The Appropriate Piece Of  Rhsz
C
  110 Continue
      Call Dgesl  (Wi, Kd, Kd, Ipvtw, Rhsdmz, 0)
      Ir = Irow
      Do 140 Jcomp = 1, Ncomp
         Mj = M(Jcomp)
         Ir = Ir + Mj
         Do 130 L = 1, Mj
            Ind = Jcomp
            Rsum = 0.d0
            Do 120 J = 1, K
               Rsum = Rsum  +  Hb(J,L) * Rhsdmz(Ind)
  120       Ind = Ind + Ncomp
            Rhsz(Ir-L) = Rsum
  130    Continue
  140 Continue
      Return
      End
C
C----------------------------------------------------------------------
C                             P A R T  4
C               Polynomial And Service Routines
C----------------------------------------------------------------------
C
      Subroutine Appsln (X, Z, Fspace, Ispace)
C
C*****************************************************************
C
C     Purpose
C
C           Set Up A Standard Call To  Approx  To Evaluate The
C           Approximate Solution  Z = Z( U(X) )  At A Point X
C           (It Has Been Computed By A Call To Colmod ).
C           The Parameters Needed For  Approx  Are Retrieved
C           From The Work Arrays  Ispace  And  Fspace .
C
C*****************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Z(*), Fspace(*), Ispace(*), A(28), Dummy(1)
      Is6 = Ispace(6)
      Is5 = Ispace(1) + 2
      Is4 = Is5 + Ispace(4) * (Ispace(1) + 1)
      I = 1
      Call Approx (I, X, Z, A, Fspace(Is6), Fspace(1), Ispace(1),
     +             Fspace(Is5), Fspace(Is4), Ispace(2), Ispace(3),
     +             Ispace(5), Ispace(8), Ispace(4), 2, Dummy, 0)
      Return
      End
      Subroutine Approx (I, X, Zval, A, Coef, Xi, N, Z, Dmz, K,
     +                   Ncomp, Mmax, M, Mstar, Mode, Dmval, Modm )
C
C**********************************************************************
C
C   Purpose
C                                    (1)       (M1-1)     (Mncomp-1)
C           Evaluate Z(U(X))=(U (X),U (X),...,U  (X),...,U  (X)      )
C                              1     1         1          Mncomp
C           At One Point X.
C
C   Variables
C     A      - Array Of Mesh Independent Rk-Basis Coefficients
C     Basm   - Array Of Mesh Dependent Monomial Coefficients
C     Xi     - The Current Mesh (Having N Subintervals)
C     Z      - The Current Solution Vector
C     Dmz    - The Array Of Mj-Th Derivatives Of The Current Solution
C     Mode   - Determines The Amount Of Initialization Needed
C            = 4  Forms Z(U(X)) Using Z, Dmz And Ha
C            = 3  As In  = 4, But Computes Local Rk-Basis
C            = 2  As In  = 3, But Determines I Such That
C                       Xi(I) .le. X .lt. Xi(I+1) (Unless X = Xi(N+1))
C            = 1  Retrieve  Z = Z(U(X(I)))  Directly
C
C**********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Zval(*), Dmval(*), Xi(*), M(*), A(7,*), Dm(7)
      Dimension Z(*), Dmz(*), Bm(4), Coef(*)
C
      Common /Colout/ Precis, Iout, Iprint
C
      Go To (10, 30, 80, 90), Mode
C
C.... Mode = 1 , Retrieve  Z( U(X) )  Directly For X = Xi(I).
C
   10 X  = Xi(I)
      Iz = (I-1) * Mstar
      Do 20 J = 1, Mstar
        Iz = Iz + 1
        Zval(J) = Z(Iz)
   20 Continue
      Return
C
C.... Mode = 2 ,  Locate I So  Xi(I) .le. X .lt. Xi(I+1)
C
   30 Continue
      If ( X .ge. Xi(1)-Precis .and. X .le. Xi(N+1)+Precis )
     +                                              Go To 40
      If (Iprint .lt. 1)  Write(Iout,1000) X, Xi(1), Xi(N+1)
      If ( X .lt. Xi(1) )  X = Xi(1)
      If ( X .gt. Xi(N+1) )  X = Xi(N+1)
   40 If ( I .gt. N .or. I .lt. 1 )  I = (N+1) / 2
      Ileft = I
      If ( X .lt. Xi(Ileft) )                       Go To 60
      Do 50 L = Ileft, N
           I = L
           If ( X .lt. Xi(L+1) )                    Go To 80
   50 Continue
      Go To 80
   60 Iright = Ileft - 1
      Do 70 L = 1, Iright
           I = Iright + 1 - L
           If ( X .ge. Xi(I) )                      Go To 80
   70 Continue
C
C.... Mode = 2 Or 3 , Compute Mesh Independent Rk-Basis.
C
   80 Continue
      S = (X - Xi(I)) / (Xi(I+1) - Xi(I))
      Call Rkbas ( S, Coef, K, Mmax, A, Dm, Modm )
C
C.... Mode = 2, 3, Or 4 , Compute Mesh Dependent Rk-Basis.
C
   90 Continue
      Bm(1) = X - Xi(I)
      Do 95 L = 2, Mmax
         Bm(L) = Bm(1) / Dfloat(L)
   95 Continue
C
C.... Evaluate  Z( U(X) ).
C
  100 Ir = 1
      Iz = (I-1) * Mstar + 1
      Idmz = (I-1) * K * Ncomp
      Do 140 Jcomp = 1, Ncomp
          Mj = M(Jcomp)
          Ir = Ir + Mj
          Iz = Iz + Mj
          Do 130 L = 1, Mj
             Ind = Idmz + Jcomp
             Zsum = 0.d0
             Do 110 J = 1, K
               Zsum = Zsum  +  A(J,L) * Dmz(Ind)
  110        Ind = Ind + Ncomp
             Do 120 Ll = 1, L
               Lb = L + 1 - Ll
  120          Zsum = Zsum * Bm(Lb)  +  Z(Iz-Ll)
  130     Zval(Ir-L) = Zsum
  140 Continue
      If ( Modm .eq. 0 )                            Return
C
C.... For Modm = 1 Evaluate  Dmval(J) = Mj-Th Derivative Of Uj.
C
      Do 150 Jcomp = 1, Ncomp
  150 Dmval(Jcomp) = 0.d0
      Idmz = Idmz + 1
      Do 170 J = 1, K
         Fact = Dm(J)
         Do 160 Jcomp = 1, Ncomp
            Dmval(Jcomp) = Dmval(Jcomp)  +  Fact * Dmz(Idmz)
            Idmz = Idmz + 1
  160    Continue
  170 Continue
      Return
C--------------------------------------------------------------------
 1000 Format(37h ****** Domain Error In Approx ******
     +       /4h X = ,D20.10, 10h   Aleft = ,D20.10,
     +       11h   Aright = ,D20.10)
      End
      Subroutine Rkbas (S, Coef, K, M, Rkb, Dm, Mode)
C
C**********************************************************************
C
C   Purpose
C           Evaluate Mesh Independent Runge-Kutta Basis For Given S
C
C   Variables
C     S      - Argument, I.e. The Relative Position For Which
C              The Basis Is To Be Evaluated ( 0. .le. S .le. 1. ).
C     Coef   - Precomputed Derivatives Of The Basis
C     K      - Number Of Collocatin Points Per Subinterval
C     M      - Maximal Order Of The Differential Equation
C     Rkb    - The Runge-Kutta Basis (0-Th To (M-1)-Th Derivatives )
C     Dm     - Basis Elements For M-Th Derivative
C
C**********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Coef(K,*), Rkb(7,*), Dm(*), T(10)
C
      If ( K .eq. 1 )                            Go To 70
      Kpm1 = K + M - 1
      Do 10 I = 1, Kpm1
   10   T(I) = S / Dfloat(I)
      Do 40 L = 1, M
         Lb = K + L + 1
         Do 30 I = 1, K
           P = Coef(1,I)
           Do 20 J = 2, K
             P = P * T(Lb-J)  + Coef(J,I)
   20      Continue
           Rkb(I,L) = P
   30    Continue
   40 Continue
      If ( Mode .eq. 0 )                         Return
      Do 60 I = 1, K
         P = Coef(1,I)
         Do 50 J = 2, K
   50       P = P * T(K+1-J) + Coef(J,I)
         Dm(I) = P
   60 Continue
      Return
   70 Rkb(1,1) = 1.0d0
      Dm(1) = 1.0d0
      Return
      End
      Subroutine Vmonde ( Rho, Coef, K )
C
C**********************************************************************
C
C   Purpose
C          Solve Vandermonde System V * X = E
C          With  V(I,J) = Rho(J)**(I-1)/(I-1)! .
C
C**********************************************************************
C
      Integer K, I,Ifac,J,Km1,Kmi
      Double Precision Rho(K), Coef(K)
C
      If ( K .eq. 1 )                             Return
      Km1 = K - 1
      Do 10 I = 1, Km1
         Kmi = K - I
         Do 10 J = 1, Kmi
           Coef(J) = (Coef(J+1) - Coef(J)) / (Rho(J+I) - Rho(J))
  10  Continue
C
      Ifac = 1
      Do 40 I = 1, Km1
         Kmi = K + 1 - I
         Do 30 J = 2, Kmi
  30        Coef(J) = Coef(J) - Rho(J+I-1) * Coef(J-1)
         Coef(Kmi) = Dfloat(Ifac) * Coef(Kmi)
         Ifac = Ifac * I
  40  Continue
      Coef(1) = Dfloat(Ifac) * Coef(1)
      Return
      End
      Subroutine Horder (I, Uhigh, Hi, Dmz, Ncomp, K)
C
C**********************************************************************
C
C   Purpose
C           Determine Highest Order (Piecewise Constant) Derivatives
C           Of The Current Collocation Solution
C
C   Variables
C     Hi     - The Stepsize, Hi = Xi(I+1) - Xi(I)
C     Dmz    - Vector Of Mj-Th Derivative Of The Solution
C     Uhigh  - The Array Of Highest Order (Piecewise Constant)
C              Derivatives Of The Approximate Solution On
C              (Xi(I),Xi(I+1)), Viz,
C                          (K+Mj-1)
C              Uhigh(J) = U   (X)    On (Xi(I),Xi(I+1))
C                          J
C
C**********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension Uhigh(*), Dmz(*)
C
      Common /Colloc/ Rho(7), Coef(49)
C
      Dn = 1.d0 / Hi**(K-1)
C
C.... Loop Over The Ncomp Solution Components
C
      Do 10 Id = 1, Ncomp
         Uhigh(Id) = 0.d0
   10 Continue
      Kin = 1
      Idmz = (I-1) * K * Ncomp + 1
      Do 30 J = 1, K
         Fact = Dn * Coef(Kin)
         Do 20 Id = 1, Ncomp
            Uhigh(Id) = Uhigh(Id)  +  Fact * Dmz(Idmz)
            Idmz = Idmz + 1
   20    Continue
         Kin = Kin + K
   30 Continue
      Return
      End
      Subroutine Dmzsol (Kd, Mstar, N, V, Z, Dmz)
C
C**********************************************************************
C
C   Purpose
C          Compute Dmz In A Blockwise Manner
C          Dmz(I) = Dmz(I)  +  V(I) * Z(I), I = 1,...,N
C
C**********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension V(Kd,*), Dmz(Kd,*), Z(*)
C
      Jz = 1
      Do 30 I = 1, N
         Do 20 J = 1, Mstar
            Fact = Z(Jz)
            Do 10 L = 1, Kd
               Dmz(L,I) = Dmz(L,I)  +  Fact * V(L,Jz)
   10       Continue
            Jz = Jz + 1
   20    Continue
   30 Continue
      Return
      End
C----------------------------------------------------------------------
C                            P A R T  5
C          We List Here A Modified (Column Oriented, Faster)
C          Version Of The Package Solveblok Of De Boor - Weiss [5].
C          We Also Give A Listing Of The Linpack
C          Routines Dgefa Und Dgesl Used By Colmod.
C----------------------------------------------------------------------
C
      Subroutine Fcblok (Bloks, Integs, Nbloks, Ipivot, Scrtch, Info)
C
C
C     Calls Subroutines  Factrb  And  Shiftb .
C
C     Fcblok  Supervises The Plu Factorization With Pivoting Of
C     Scaled Rows Of The Almost Block Diagonal Matrix Stored In The
C     Arrays  Bloks  And  Integs .
C
C     Factrb = Subprogram Which Carries Out Steps 1,...,Last Of Gauss
C            Elimination (With Pivoting) For An Individual Block.
C     Shiftb = Subprogram Which Shifts The Remaining Rows To The Top Of
C            The Next Block
C
C     Parameters
C      Bloks   An Array That Initially Contains The Almost Block Diago-
C            Nal Matrix  A  To Be Factored, And On Return Contains The
C            Computed Factorization Of  A .
C      Integs  An Integer Array Describing The Block Structure Of  A .
C      Nbloks  The Number Of Blocks In  A .
C      Ipivot  An Integer Array Of Dimension   Sum (Integs(3,N) ; N = 1,
C            ...,Nbloks) Which, On Return, Contains The Pivoting Stra-
C            Tegy Used.
C      Scrtch  Work Area Required, Of Length  Max (Integs(1,N) ; N = 1,
C            ...,Nbloks).
C      Info    Output Parameter;
C            = 0  In Case Matrix Was Found To Be Nonsingular.
C            Otherwise,
C            = N If The Pivot Element In The Nth Gauss Step Is Zero.
C
C**********************************************************************
C
      Integer Integs(3,Nbloks),Ipivot(*),Info, I,Index,Indexn,Last,
     +        Ncol,Nrow
      Double Precision Bloks(*),Scrtch(*)
      Info = 0
      Indexx = 1
      Indexn = 1
      I = 1
C
C.... Loop Over The Blocks.  I  Is Loop Index
C
   10      Index = Indexn
           Nrow = Integs(1,I)
           Ncol = Integs(2,I)
           Last = Integs(3,I)
C
C....     Carry Out Elimination On The I-Th Block Until Next Block
C....     Enters, I.e., For Columns 1,...,Last  Of I-Th Block.
C
           Call Factrb ( Bloks(Index), Ipivot(Indexx), Scrtch, Nrow,
     +                   Ncol, Last, Info)
C
C....     Check For Having Reached A Singular Block Or The Last Block
C
           If ( Info .ne. 0 )                       Go To 20
           If ( I .eq. Nbloks )                     Return
           I = I+1
           Indexn = Nrow * Ncol + Index
           Indexx = Indexx + Last
C
C....     Put The Rest Of The I-Th Block Onto The Next Block
C
           Call Shiftb ( Bloks(Index), Nrow, Ncol, Last,
     +                   Bloks(Indexn), Integs(1,I), Integs(2,I) )
      Go To 10
   20 Info = Info + Indexx - 1
      Return
      End
      Subroutine Factrb ( W, Ipivot, D, Nrow, Ncol, Last, Info)
C
C********************************************************************
C
C     Adapted From P.132 Of  Element.numer.analysis  By Conte-De Boor
C
C     Constructs A Partial Plu Factorization, Corresponding To Steps
C      1,..., Last   In Gauss Elimination, For The Matrix  W  Of
C      Order ( Nrow ,  Ncol ), Using Pivoting Of Scaled Rows.
C
C     Parameters
C       W       Contains The (Nrow,Ncol) Matrix To Be Partially Factored
C               On Input, And The Partial Factorization On Output.
C       Ipivot  An Integer Array Of Length Last Containing A Record Of
C               The Pivoting Strategy Used; Explicit Interchanges
C               Are Used For Pivoting.
C       D       A Work Array Of Length Nrow Used To Store Row Sizes
C               Temporarily.
C       Nrow    Number Of Rows Of W.
C       Ncol    Number Of Columns Of W.
C       Last    Number Of Elimination Steps To Be Carried Out.
C       Info    On Output, Zero If The Matrix Is Found To Be Non-
C               Singular, In Case A Zero Pivot Was Encountered In Row
C               N,  Info = N On Output.
C
C**********************************************************************
C
      Integer Ipivot(Nrow),Ncol,Last,Info, I,J,K,L,Kp1
      Double Precision W(Nrow,Ncol),D(Nrow), Colmax,T,S
      Double Precision Abs,Max
C
C.... Initialize  D
C
      Do 10 I = 1, Nrow
        D(I) = 0.d0
   10 Continue
      Do 20 J = 1, Ncol
        Do 20 I = 1, Nrow
          D(I) = Max( D(I) , Abs(W(I,J)))
   20 Continue
C
C.... Gauss Elimination With Pivoting Of Scaled Rows, Loop Over
C.... K = 1,.,Last
C
      K = 1
C
C.... As Pivot Row For K-Th Step, Pick Among The Rows Not Yet Used,
C.... I.e., From Rows  K ,..., Nrow , The One Whose K-Th Entry
C.... (Compared To The Row Size) Is Largest. Then, If This Row
C.... Does Not Turn Out To Be Row K, Interchange Row K With This
C.... Particular Row And Redefine Ipivot(K).
C
   30      Continue
           If ( D(K) .eq. 0.d0 )                    Go To 90
           If (K .eq. Nrow)                         Go To 80
           L = K
           Kp1 = K+1
           Colmax = Abs(W(K,K)) / D(K)
C
C....     Find The (Relatively) Largest Pivot
C
           Do 40 I = Kp1, Nrow
             If ( Abs(W(I,K)) .le. Colmax * D(I) ) Go To 40
             Colmax = Abs(W(I,K)) / D(I)
             L = I
   40      Continue
           Ipivot(K) = L
           T = W(L,K)
           S = D(L)
           If ( L .eq. K )                          Go To 50
             W(L,K) = W(K,K)
             W(K,K) = T
             D(L) = D(K)
             D(K) = S
   50      Continue
C
C....     If Pivot Element Is Too Small In Absolute Value, Declare
C....     Matrix To Be Noninvertible And Quit.
C
           If ( Abs(T)+D(K) .le. D(K) )            Go To 90
C
C....     Otherwise, Subtract The Appropriate Multiple Of The Pivot
C....     Row From Remaining Rows, I.e., The Rows (K+1),..., (Nrow)
C....     To Make K-Th Entry Zero. Save The Multiplier In Its Place.
C....     For High Performance Do This Operations Column Oriented.
C
           T = -1.0d0/T
           Do 60 I = Kp1, Nrow
   60        W(I,K) = W(I,K) * T
           Do 70 J = Kp1,Ncol
             T = W(L,J)
             If ( L .eq. K )                        Go To 62
               W(L,J) = W(K,J)
               W(K,J) = T
   62        If ( T .eq. 0.d0 )                     Go To 70
             Do 64 I = Kp1, Nrow
   64           W(I,J) = W(I,J) + W(I,K) * T
   70      Continue
           K = Kp1
C
C.... Check For Having Reached The Next Block.
C
           If ( K .le. Last )                       Go To 30
      Return
C
C.... If  Last  .eq. Nrow , Check Now That Pivot Element In Last Row
C.... Is Nonzero.
C
   80 If( Abs(W(Nrow,Nrow))+D(Nrow) .gt. D(Nrow) ) Return
C
C.... Singularity Flag Set
C
   90 Info = K
      Return
      End
      Subroutine Shiftb (Ai, Nrowi, Ncoli, Last, Ai1, Nrowi1, Ncoli1)
C
C*********************************************************************
C
C     Shifts The Rows In Current Block, Ai, Not Used As Pivot Rows, If
C     Any, I.e., Rows  (Last+1),..., (Nrowi), Onto The First Mmax  = 
C      = Nrow-Last  Rows Of The Next Block, Ai1, With Column Last+J Of
C      Ai Going To Column J , J = 1,...,Jmax = Ncoli-Last. The Remaining
C     Columns Of These Rows Of Ai1 Are Zeroed Out.
C
C                                Picture
C
C          Original Situation After         Results In A New Block I+1
C          Last = 2 Columns Have Been       Created And Ready To Be
C          Done In Factrb (Assuming No      Factored By Next Factrb
C          Interchanges Of Rows)            Call.
C                      1
C                 X  X 1x  X  X           X  X  X  X  X
C                      1
C                 0  X 1x  X  X           0  X  X  X  X
C     Block I          1                       ---------------
C     Nrowi = 4   0  0 1x  X  X           0  0 1x  X  X  0  01
C     Ncoli = 5        1                       1             1
C     Last = 2    0  0 1x  X  X           0  0 1x  X  X  0  01
C     -------------------------------          1             1   New
C                      1x  X  X  X  X          1x  X  X  X  X1  Block
C                      1                       1             1   I+1
C     Block I+1        1x  X  X  X  X          1x  X  X  X  X1
C     Nrowi1 =  5        1                       1             1
C     Ncoli1 =  5        1x  X  X  X  X          1x  X  X  X  X1
C     -------------------------------          1-------------1
C                      1
C
C*********************************************************************
C
      Integer Last, J,Jmax,Jmaxp1,M,Mmax
      Double Precision Ai(Nrowi,Ncoli),Ai1(Nrowi1,Ncoli1)
      Mmax = Nrowi - Last
      Jmax = Ncoli - Last
      If (Mmax .lt. 1 .or. Jmax .lt. 1)             Return
C
C.... Put The Remainder Of Block I Into Ai1
C
      Do 10 J = 1,Jmax
           Do 10 M = 1,Mmax
   10 Ai1(M,J) = Ai(Last+M,Last+J)
      If (Jmax .eq. Ncoli1)                         Return
C
C.... Zero Out The Upper Right Corner Of Ai1
C
      Jmaxp1 = Jmax + 1
      Do 20 J = Jmaxp1,Ncoli1
           Do 20 M = 1,Mmax
   20 Ai1(M,J) = 0.d0
      Return
      End
      Subroutine Sbblok ( Bloks, Integs, Nbloks, Ipivot, X )
C
C**********************************************************************
C
C     Calls Subroutines  Subfor  And  Subbak .
C
C     Supervises The Solution (By Forward And Backward Substitution) Of
C     The Linear System  A*X = B  For X, With The Plu Factorization Of
C     A  Already Generated In  Fcblok .  Individual Blocks Of
C     Equations Are Solved Via  Subfor  And  Subbak .
C
C    Parameters
C       Bloks, Integs, Nbloks, Ipivot    Are As On Return From Fcblok.
C       X       On Input: The Right Hand Side, In Dense Storage
C               On Output: The Solution Vector
C
C*********************************************************************
C
      Integer Integs(3,Nbloks),Ipivot(*), I,Index,Indexx,J,Last,
     +        Nbp1,Ncol,Nrow
      Double Precision Bloks(*), X(*)
C
C.... Forward Substitution Pass
C
      Index = 1
      Indexx = 1
      Do 10 I = 1, Nbloks
           Nrow = Integs(1,I)
           Last = Integs(3,I)
           Call Subfor ( Bloks(Index), Ipivot(Indexx), Nrow, Last,
     +                   X(Indexx) )
           Index = Nrow * Integs(2,I) + Index
   10 Indexx = Indexx + Last
C
C.... Back Substitution Pass
C
      Nbp1 = Nbloks + 1
      Do 20 J = 1, Nbloks
           I = Nbp1 - J
           Nrow = Integs(1,I)
           Ncol = Integs(2,I)
           Last = Integs(3,I)
           Index = Index - Nrow * Ncol
           Indexx = Indexx - Last
   20 Call Subbak ( Bloks(Index), Nrow, Ncol, Last, X(Indexx) )
      Return
      End
      Subroutine Subfor ( W, Ipivot, Nrow, Last, X )
C
C**********************************************************************
C
C     Carries Out The Forward Pass Of Substitution For The Current
C     Block, I.e., The Action On The Right Side Corresponding To The
C     Elimination Carried Out In  Factrb  For This Block.
C
C    Parameters
C       W, Ipivot, Nrow, Last  Are As On Return From Factrb.
C       X(J)  Is Expected To Contain, On Input, The Right Side Of J-Th
C             Equation For This Block, J = 1,...,Nrow.
C       X(J)  Contains, On Output, The Appropriately Modified Right
C             Side Of Equation (J) In This Block, J = 1,...,Last And
C             For J = Last+1,...,Nrow.
C
C*********************************************************************
C
      Integer Ipivot(Last), Ip,K,Kp1,Lstep
      Double Precision W(Nrow,Last), X(Nrow), T
C
      If ( Nrow .eq. 1 )                            Return
      Lstep = Min( Nrow-1 , Last )
      Do 20 K = 1, Lstep
           Kp1 = K + 1
           Ip = Ipivot(K)
           T = X(Ip)
           X(Ip) = X(K)
           X(K) = T
           If ( T .eq. 0.d0 )                       Go To 20
           Do 10 I = Kp1, Nrow
   10         X(I) = X(I) + W(I,K) * T
   20 Continue
   30 Return
      End
      Subroutine Subbak ( W, Nrow, Ncol, Last, X )
C
C*********************************************************************
C
C     Carries Out Backsubstitution For Current Block.
C
C    Parameters
C       W, Ipivot, Nrow, Ncol, Last  Are As On Return From Factrb.
C       X(1),...,X(Ncol)  Contains, On Input, The Right Side For The
C               Equations In This Block After Backsubstitution Has Been
C               Carried Up To But Not Including Equation (Last).
C               Means That X(J) Contains The Right Side Of Equation (J)
C               As Modified During Elimination, J = 1,...,Last, While
C               For J .gt. Last, X(J) Is Already A Component Of The
C               Solution Vector.
C       X(1),...,X(Ncol) Contains, On Output, The Components Of The
C               Solution Corresponding To The Present Block.
C
C*********************************************************************
C
      Integer  Last,  I,J,K,Km1,Lm1,Lp1
      Double Precision W(Nrow,Ncol),X(Ncol), T
C
      Lp1 = Last + 1
      If ( Lp1 .gt. Ncol )                          Go To 30
      Do 20 J = Lp1, Ncol
         T = - X(J)
         If ( T .eq. 0.d0 )                         Go To 20
         Do 10 I = 1, Last
   10       X(I) = X(I) + W(I,J) * T
   20 Continue
   30 If ( Last .eq. 1 )                            Go To 60
      Lm1 = Last - 1
      Do 50 Kb = 1, Lm1
        Km1 = Last - Kb
        K = Km1 + 1
        X(K) = X(K)/W(K,K)
        T = - X(K)
        If ( T .eq. 0.d0 )                          Go To 50
        Do 40 I = 1, Km1
   40     X(I) = X(I) + W(I,K) * T
   50 Continue
   60 X(1) = X(1)/W(1,1)
      Return
      End
      Subroutine Dgefa(A,Lda,N,Ipvt,Info)
      Integer Lda,N,Ipvt(*),Info
      Double Precision A(Lda,*)
C
C     Dgefa Factors A Double Precision Matrix By Gaussian Elimination.
C
C     Dgefa Is Usually Called By Dgeco, But It Can Be Called
C     Directly With A Saving In Time If  Rcond  Is Not Needed.
C     (Time For Dgeco) = (1 + 9/N)*(Time For Dgefa) .
C
C     On Entry
C
C        A       Double Precision(Lda, N)
C                The Matrix To Be Factored.
C
C        Lda     Integer
C                The Leading Dimension Of The Array  A .
C
C        N       Integer
C                The Order Of The Matrix  A .
C
C     On Return
C
C        A       An Upper Triangular Matrix And The Multipliers
C                Which Were Used To Obtain It.
C                The Factorization Can Be Written  A = L*U  Where
C                L  Is A Product Of Permutation And Unit Lower
C                Triangular Matrices And  U  Is Upper Triangular.
C
C        Ipvt    Integer(N)
C                An Integer Vector Of Pivot Indices.
C
C        Info    Integer
C                = 0  Normal Value.
C                = K  If  U(K,K) .eq. 0.0 .  This Is Not An Error
C                     Condition For This Subroutine, But It Does
C                     Indicate That Dgesl Or Dgedi Will Divide By Zero
C                     If Called.  Use  Rcond  In Dgeco For A Reliable
C                     Indication Of Singularity.
C
C     Linpack. This Version Dated 08/14/78 .
C     Cleve Moler, University Of New Mexico, Argonne National Lab.
C
C     Subroutines And Functions
C
C     Blas Daxpy,Dscal,Idamax
C
C     Internal Variables
C
      Double Precision T
      Integer Idamax,J,K,Kp1,L,Nm1
C
C
C     Gaussian Elimination With Partial Pivoting
C
      Info = 0
      Nm1 = N - 1
      If (Nm1 .lt. 1) Go To 70
      Do 60 K = 1, Nm1
         Kp1 = K + 1
C
C        Find L = Pivot Index
C
         L = Idamax(N-K+1,A(K,K),1) + K - 1
         Ipvt(K) = L
C
C        Zero Pivot Implies This Column Already Triangularized
C
         If (A(L,K) .eq. 0.0d0) Go To 40
C
C           Interchange If Necessary
C
            If (L .eq. K) Go To 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       Continue
C
C           Compute Multipliers
C
            T = -1.0d0/A(K,K)
            Call Dscal(N-K,T,A(K+1,K),1)
C
C           Row Elimination With Column Indexing
C
            Do 30 J = Kp1, N
               T = A(L,J)
               If (L .eq. K) Go To 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          Continue
               Call Daxpy(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       Continue
         Go To 50
   40    Continue
            Info = K
   50    Continue
   60 Continue
   70 Continue
      Ipvt(N) = N
      If (A(N,N) .eq. 0.0d0) Info = N
      Return
      End
      Subroutine Dgesl(A,Lda,N,Ipvt,B,Job)
      Integer Lda,N,Ipvt(*),Job
      Double Precision A(Lda,*),B(*)
C
C     Dgesl Solves The Double Precision System
C     A * X = B  Or  Trans(A) * X = B
C     Using The Factors Computed By Dgeco Or Dgefa.
C
C     On Entry
C
C        A       Double Precision(Lda, N)
C                The Output From Dgeco Or Dgefa.
C
C        Lda     Integer
C                The Leading Dimension Of The Array  A .
C
C        N       Integer
C                The Order Of The Matrix  A .
C
C        Ipvt    Integer(N)
C                The Pivot Vector From Dgeco Or Dgefa.
C
C        B       Double Precision(N)
C                The Right Hand Side Vector.
C
C        Job     Integer
C                = 0         To Solve  A*X = B ,
C                = Nonzero   To Solve  Trans(A)*X = B  Where
C                            Trans(A)  Is The Transpose.
C
C     On Return
C
C        B       The Solution Vector  X .
C
C     Error Condition
C
C        A Division By Zero Will Occur If The Input Factor Contains A
C        Zero On The Diagonal.  Technically This Indicates Singularity
C        But It Is Often Caused By Improper Arguments Or Improper
C        Setting Of Lda .  It Will Not Occur If The Subroutines Are
C        Called Correctly And If Dgeco Has Set Rcond .gt. 0.0
C        Or Dgefa Has Set Info .eq. 0 .
C
C     To Compute  Inverse(A) * C  Where  C  Is A Matrix
C     With  P  Columns
C           Call Dgeco(A,Lda,N,Ipvt,Rcond,Z)
C           If (Rcond Is Too Small) Go To ...
C           Do 10 J = 1, P
C              Call Dgesl(A,Lda,N,Ipvt,C(1,J),0)
C        10 Continue
C
C     Linpack. This Version Dated 08/14/78 .
C     Cleve Moler, University Of New Mexico, Argonne National Lab.
C
C     Subroutines And Functions
C
C     Blas Daxpy,Ddot
C
C     Internal Variables
C
      Double Precision Ddot,T
      Integer K,Kb,L,Nm1
C
      Nm1 = N - 1
      If (Job .ne. 0) Go To 50
C
C        Job = 0 , Solve  A * X = B
C        First Solve  L*Y = B
C
         If (Nm1 .lt. 1) Go To 30
         Do 20 K = 1, Nm1
            L = Ipvt(K)
            T = B(L)
            If (L .eq. K) Go To 10
               B(L) = B(K)
               B(K) = T
   10       Continue
            Call Daxpy(N-K,T,A(K+1,K),1,B(K+1),1)
   20    Continue
   30    Continue
C
C        Now Solve  U*X = Y
C
         Do 40 Kb = 1, N
            K = N + 1 - Kb
            B(K) = B(K)/A(K,K)
            T = -B(K)
            Call Daxpy(K-1,T,A(1,K),1,B(1),1)
   40    Continue
      Go To 100
   50 Continue
C
C        Job = Nonzero, Solve  Trans(A) * X = B
C        First Solve  Trans(U)*Y = B
C
         Do 60 K = 1, N
            T = Ddot(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    Continue
C
C        Now Solve Trans(L)*X = Y
C
         If (Nm1 .lt. 1) Go To 90
         Do 80 Kb = 1, Nm1
            K = N - Kb
            B(K) = B(K) + Ddot(N-K,A(K+1,K),1,B(K+1),1)
            L = Ipvt(K)
            If (L .eq. K) Go To 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       Continue
   80    Continue
   90    Continue
  100 Continue
      Return
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
