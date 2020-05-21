
!IJB Adapted to 3d grids and FD template from https://web.stanford.edu/group/SOL/software/minres/

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File minresModule.f90
!
! MINRES solves symmetric systems Ax = b or min ||Ax - b||_2,
! where the matrix A may be indefinite and/or singular.
!
! Contributors:
!     Chris Paige <chris@cs.mcgill.ca>
!     Sou-Cheng Choi <scchoi@stanford.edu>
!
!     Michael Saunders <saunders@stanford.edu>
!     Systems Optimization Laboratory (SOL)
!     Stanford University
!     Stanford, CA 94305-4026, USA
!     (650)723-1875
!
! 09 Oct 2007: F90 version constructed from the F77 version.
!              Initially used compiler option -r8, but this is nonstandard.
! 15 Oct 2007: Test on Arnorm = ||Ar|| added to recognize singular systems.
! 15 Oct 2007: Temporarily used real(8) everywhere.
! 16 Oct 2007: Use minresDataModule to define dp = selected_real_kind(15).
!              We need "use minresDataModule"
!              at the beginning of modules AND inside interfaces.
!
!              g95 compiles successfully with the following options:
!   g95 -c -g -O0 -pedantic -Wall -Wextra -fbounds-check -ftrace=full minresModule.f90
! Apr-May 2020: Adapted by Ian Bush (Oxford University) for us in SSP Poisson solver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module equation_solver_minres_module

  Use equation_solver_base_class_module, Only : equation_solver_base_class
  
  Implicit None

  Type, Public, Extends( equation_solver_base_class ) :: equation_solver_minres
   Contains
     Procedure, Public :: solve => minres
  End type equation_solver_minres

  Private
  
Contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Subroutine MINRES( method, &
       lb, ub, FD_operator, comms, halo_swapper, Msolve, b, itnlim, rtol,  precon, &
       x, istop, istop_message, itn, rnorm )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use comms_base_class_module, Only : comms_base_class
    Use halo_setter_base_module, Only : halo_setter_base_class
    Use FD_template_module     , Only : FD_template

    Class( equation_solver_minres )                       , Intent( In    ) :: method
    Integer,  Dimension( 1:3 )                            , Intent( In    ) :: lb( 1:3 )
    Integer,  Dimension( 1:3 )                            , Intent( In    ) :: ub( 1:3 )
    Class( FD_template )                                  , Intent( In    ) :: FD_operator
    Class( halo_setter_base_class )                       , Intent( InOut ) :: halo_swapper
    Class( comms_base_class       )                       , Intent( In    ) :: comms
    Real( wp ) , Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ) :: b
    Integer                                               , Intent( In    ) :: itnlim
    Real( wp )                                            , Intent( In    ) :: rtol
    Logical                                               , Intent( In    ) :: precon
    Real( wp ) , Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: x
    Integer                                               , Intent(   Out ) :: istop
    Character( Len = * )                                  , Intent(   Out ) :: istop_message
    Real( wp )                                            , Intent(   Out ) :: rnorm
    Integer                                               , Intent(   Out ) :: itn
    Interface
       Subroutine Msolve(lb,ub,x,y)                   ! Solve M*y = x
         Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
         Integer                                               , Intent( In    ) :: lb( 1:3 )
         Integer                                               , Intent( In    ) :: ub( 1:3 )
         Real( wp ) , Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ) :: x
         Real( wp ) , Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: y
       End Subroutine Msolve
    End Interface

    Logical :: checkA = .False.
    Real( wp ) :: shift = 0.0_wp
    Real(wp) :: Anorm, Acond, Arnorm, ynorm
    !-------------------------------------------------------------------
    !
    ! MINRES  is designed to solve the system of linear equations
    !
    !    Ax = b
    !
    ! or the least-squares problem
    !
    !    min ||Ax - b||_2,
    !
    ! where A is an n by n symmetric matrix and b is a given vector.
    ! The matrix A may be indefinite and/or singular.
    !
    ! 1. If A is known to be positive definite, the Conjugate Gradient
    ! Method might be preferred, since it requires the same number
    ! of iterations as MINRES but less work per iteration.
    !
    ! 2. If A is indefinite but Ax = b is known to have a solution
    ! (e.g. if A is nonsingular), SYMMLQ might be preferred,
    ! since it requires the same number of iterations as MINRES
    ! but slightly less work per iteration.
    !
    ! The matrix A is intended to be large and sparse.  It is accessed
    ! by means of a subroutine call of the form
    ! SYMMLQ development:
    !
    !    call Aprod ( n, x, y )
    !
    ! which must return the product y = Ax for any given vector x.
    !
    !
    ! More generally, MINRES is designed to solve the system
    !
    !    (A - shift*I) x = b
    ! or
    !    min ||(A - shift*I) x - b||_2,
    !
    ! where  shift  is a specified scalar value.  Again, the matrix
    ! (A - shift*I) may be indefinite and/or singular.
    ! The work per iteration is very slightly less if  shift = 0.
    !
    ! Note: If  shift  is an approximate eigenvalue of  A
    ! and  b  is an approximate eigenvector,  x  might prove to be
    ! a better approximate eigenvector, as in the methods of
    ! inverse iteration and/or Rayleigh-quotient iteration.
    ! However, we're not yet sure on that -- it may be better to use SYMMLQ.
    !
    ! A further option is that of preconditioning, which may reduce
    ! the number of iterations required.  If M = C C' is a positive
    ! definite matrix that is known to approximate  (A - shift*I)
    ! in some sense, and if systems of the form  My = x  can be
    ! solved efficiently, the parameters precon and Msolve may be
    ! used (see below).  When  precon = .true., MINRES will
    ! implicitly solve the system of equations
    !
    !    P (A - shift*I) P' xbar  =  P b,
    !
    ! i.e.             Abar xbar  =  bbar
    ! where                    P  =  C**(-1),
    !                       Abar  =  P (A - shift*I) P',
    !                       bbar  =  P b,
    !
    ! and return the solution       x  =  P' xbar.
    ! The associated residual is rbar  =  bbar - Abar xbar
    !                                  =  P (b - (A - shift*I)x)
    !                                  =  P r.
    !
    ! In the discussion below, eps refers to the machine precision.
    !
    ! Parameters
    ! ----------
    !
    ! n       input      The dimension of the matrix A.
    ! b(n)    input      The rhs vector b.
    ! x(n)    output     Returns the computed solution x.
    !
    ! Aprod   external   A subroutine defining the matrix A.
    !                       call Aprod ( n, x, y )
    !                    must return the product y = Ax
    !                    without altering the vector x.
    !
    ! Msolve  external   An optional subroutine defining a
    !                    preconditioning matrix M, which should
    !                    approximate (A - shift*I) in some sense.
    !                    M must be positive definite.
    !
    !                       call Msolve( n, x, y )
    !
    !                    must solve the linear system My = x
    !                    without altering the vector x.
    !
    !                    In general, M should be chosen so that Abar has
    !                    clustered eigenvalues.  For example,
    !                    if A is positive definite, Abar would ideally
    !                    be close to a multiple of I.
    !                    If A or A - shift*I is indefinite, Abar might
    !                    be close to a multiple of diag( I  -I ).
    !
    ! checkA  input      If checkA = .true., an extra call of Aprod will
    !                    be used to check if A is symmetric.  Also,
    !                    if precon = .true., an extra call of Msolve
    !                    will be used to check if M is symmetric.
    !
    ! precon  input      If precon = .true., preconditioning will
    !                    be invoked.  Otherwise, subroutine Msolve
    !                    will not be referenced; in this case the
    !                    actual parameter corresponding to Msolve may
    !                    be the same as that corresponding to Aprod.
    !
    ! shift   input      Should be zero if the system Ax = b is to be
    !                    solved.  Otherwise, it could be an
    !                    approximation to an eigenvalue of A, such as
    !                    the Rayleigh quotient b'Ab / (b'b)
    !                    corresponding to the vector b.
    !                    If b is sufficiently like an eigenvector
    !                    corresponding to an eigenvalue near shift,
    !                    then the computed x may have very large
    !                    components.  When normalized, x may be
    !                    closer to an eigenvector than b.
    !
    ! nout    input      A file number.
    !                    If nout > 0, a summary of the iterations
    !                    will be printed on unit nout.
    !
    ! itnlim  input      An upper limit on the number of iterations.
    !
    ! rtol    input      A user-specified tolerance.  MINRES terminates
    !                    if it appears that norm(rbar) is smaller than
    !                       rtol * norm(Abar) * norm(xbar),
    !                    where rbar is the transformed residual vector,
    !                       rbar = bbar - Abar xbar.
    !
    !                    If shift = 0 and precon = .false., MINRES
    !                    terminates if norm(b - A*x) is smaller than
    !                       rtol * norm(A) * norm(x).
    !
    ! istop   output     An integer giving the reason for termination...
    !
    !          -1        beta2 = 0 in the Lanczos iteration; i.e. the
    !                    second Lanczos vector is zero.  This means the
    !                    rhs is very special.
    !                    If there is no preconditioner, b is an
    !                    eigenvector of A.
    !                    Otherwise (if precon is true), let My = b.
    !                    If shift is zero, y is a solution of the
    !                    generalized eigenvalue problem Ay = lambda My,
    !                    with lambda = alpha1 from the Lanczos vectors.
    !
    !                    In general, (A - shift*I)x = b
    !                    has the solution         x = (1/alpha1) y
    !                    where My = b.
    !
    !           0        b = 0, so the exact solution is x = 0.
    !                    No iterations were performed.
    !
    !           1        Norm(rbar) appears to be less than
    !                    the value  rtol * norm(Abar) * norm(xbar).
    !                    The solution in  x  should be acceptable.
    !
    !           2        Norm(rbar) appears to be less than
    !                    the value  eps * norm(Abar) * norm(xbar).
    !                    This means that the residual is as small as
    !                    seems reasonable on this machine.
    !
    !           3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
    !                    which should indicate that x has essentially
    !                    converged to an eigenvector of A
    !                    corresponding to the eigenvalue shift.
    !
    !           4        Acond (see below) has exceeded 0.1/eps, so
    !                    the matrix Abar must be very ill-conditioned.
    !                    x may not contain an acceptable solution.
    !
    !           5        The iteration limit was reached before any of
    !                    the previous criteria were satisfied.
    !
    !           6        The matrix defined by Aprod does not appear
    !                    to be symmetric.
    !                    For certain vectors y = Av and r = Ay, the
    !                    products y'y and r'v differ significantly.
    !
    !           7        The matrix defined by Msolve does not appear
    !                    to be symmetric.
    !                    For vectors satisfying My = v and Mr = y, the
    !                    products y'y and r'v differ significantly.
    !
    !           8        An inner product of the form  x' M**(-1) x
    !                    was not positive, so the preconditioning matrix
    !                    M does not appear to be positive definite.
    !
    !                    If istop >= 5, the final x may not be an
    !                    acceptable solution.
    !
    ! itn     output     The number of iterations performed.
    !
    ! Anorm   output     An estimate of the norm of the matrix operator
    !                    Abar = P (A - shift*I) P',   where P = C**(-1).
    !
    ! Acond   output     An estimate of the condition of Abar above.
    !                    This will usually be a substantial
    !                    under-estimate of the true condition.
    !
    ! rnorm   output     An estimate of the norm of the final
    !                    transformed residual vector,
    !                       P (b  -  (A - shift*I) x).
    !
    ! ynorm   output     An estimate of the norm of xbar.
    !                    This is sqrt( x'Mx ).  If precon is false,
    !                    ynorm is an estimate of norm(x).
    !-------------------------------------------------------------------
    ! MINRES is an implementation of the algorithm described in
    ! the following reference:
    !
    ! C. C. Paige and M. A. Saunders (1975),
    ! Solution of sparse indefinite systems of linear equations,
    ! SIAM J. Numer. Anal. 12(4), pp. 617-629.
    !-------------------------------------------------------------------
    !
    !
    ! MINRES development:
    !    1972: First version, similar to original SYMMLQ.
    !          Later lost @#%*!
    !    Oct 1995: Tried to reconstruct MINRES from
    !              1995 version of SYMMLQ.
    ! 30 May 1999: Need to make it more like LSQR.
    !              In middle of major overhaul.
    ! 19 Jul 2003: Next attempt to reconstruct MINRES.
    !              Seems to need two vectors more than SYMMLQ.  (w1, w2)
    !              Lanczos is now at the top of the loop,
    !              so the operator Aprod is called in just one place
    !              (not counting the initial check for symmetry).
    ! 22 Jul 2003: Success at last.  Preconditioning also works.
    !              minres.f added to http://www.stanford.edu/group/SOL/.
    !
    ! 16 Oct 2007: Added a stopping rule for singular systems,
    !              as derived in Sou-Cheng Choi's PhD thesis.
    !              Note that ||Ar|| small => r is a null vector for A.
    !              Subroutine minrestest2 in minresTestModule.f90
    !              tests this option.  (NB: Not yet working.)
    !-------------------------------------------------------------------

    !     Local arrays and variables

    Real( wp ), Dimension( :, :, : ), Allocatable :: grid_with_halo

    Real( wp ), Dimension( :, :, : ), Allocatable :: r1
    Real( wp ), Dimension( :, :, : ), Allocatable :: r2
    Real( wp ), Dimension( :, :, : ), Allocatable :: v
    Real( wp ), Dimension( :, :, : ), Allocatable :: w
    Real( wp ), Dimension( :, :, : ), Allocatable :: w1
    Real( wp ), Dimension( :, :, : ), Allocatable :: w2
    Real( wp ), Dimension( :, :, : ), Allocatable :: y
    Real(wp)  :: alfa  , beta  , beta1 , cs    ,          &
         dbar  , delta , denom , diag  ,          &
         eps   , epsa  , epsln , epsr  , epsx  ,  &
         gamma , gbar  , gmax  , gmin  ,          &
         oldb  , oldeps, qrnorm, phi   , phibar,  &
         rhs1  , rhs2  , rnorml, rootl ,          &
         Arnorml,        relArnorml,              &
         s     , sn    , t     , tnorm2, ynorm2, z

    Integer :: FD_order, halo_width

    Logical   :: debug, prnt
    Integer :: error

    ! Local constants
    Real(wp),         Parameter :: zero =  0.0_wp,  one = 1.0_wp
    Real(wp),         Parameter :: ten  = 10.0_wp
    Character(len=*), Parameter :: msg(-1:8) =                  &
         (/ 'beta2 = 0.  If M = I, b and x are eigenvectors of A', & ! -1
         'beta1 = 0.  The exact solution is  x = 0           ', & !  0
         'Requested accuracy achieved, as determined by rtol ', & !  1
         'Reasonable accuracy achieved, given eps            ', & !  2
         'x has converged to an eigenvector                  ', & !  3
         'Acond has exceeded 0.1/eps                         ', & !  4
         'The iteration limit was reached                    ', & !  5
         'Aprod  does not define a symmetric matrix          ', & !  6
         'Msolve does not define a symmetric matrix          ', & !  7
         'Msolve does not define a pos-def preconditioner    ' /) !  8
    !-------------------------------------------------------------------

    Intrinsic       :: abs, dot_product, epsilon, min, max, sqrt

    solver_block: Block

      Allocate( r1, r2, v, w, w1, w2, y, Mold = b )

      ! Print heading and initialize.

      debug = .False.
      eps   = Epsilon(eps)
      istop  = 0
      itn    = 0
      Anorm  = zero
      Acond  = zero
      rnorm  = zero
      ynorm  = zero
      x      = zero

      !IJB
      ! Note order is always even, so no worries about splitting it in 2
      FD_order  = FD_operator%get_order()

!!$    halo_width = FD_order / 2
      ! Looks like a bug in get order! Returns twice what expected ...
      halo_width = FD_order
      Call halo_swapper%allocate( lb, ub, halo_width, grid_with_halo )

      ! Otherwise Arnorml may be used uninitiliased if we have a quick exit
      ! due to errors before the iteration starts
      Arnorml = Huge( Arnorml )

      !-------------------------------------------------------------------
      ! Set up y and v for the first Lanczos vector v1.
      ! y = beta1 P' v1, where P = C**(-1).
      ! v is really P' v1.
      !-------------------------------------------------------------------
      r1     = b
      y      = b
      If ( precon ) Call Msolve( lb, ub, b, y )
      beta1 = method%contract( comms, b, y )

      If (beta1 < zero) Then     ! M must be indefinite.
         istop = 8
         Exit solver_block
      End If

      If (Abs(beta1) < Tiny(beta1)) Then    ! b = 0 exactly.  Stop with x = 0.
         istop = 0
         Exit solver_block
      End If

      beta1  = Sqrt( beta1 )     ! Normalize y to get v1 later.

      !-------------------------------------------------------------------
      ! See if Msolve is symmetric.
      !-------------------------------------------------------------------
      If (checkA  .And.  precon) Then
         Call Msolve( lb, ub, y, r2 )
         s = method%contract( comms, y , y  )
         t = method%contract( comms, r1, r2 )
         z      = Abs(s - t)
         epsa   = (s + eps) * eps**0.33333
         If (z > epsa) Then
            istop = 7
            Exit solver_block
         End If
      End If

      !-------------------------------------------------------------------
      ! See if Aprod  is symmetric.  Initialize Arnorm.
      !-------------------------------------------------------------------
      If (checkA) Then
         ! NEED TO FIX ARGUMENTS - especially lb of grid
         Call halo_swapper%fill( halo_width, Lbound( grid_with_halo ), y, grid_with_halo, error )
         Call FD_operator%apply( Lbound( grid_with_halo ), Lbound( w  ), Lbound( w  ), Ubound( w  ), &
              grid_with_halo, w  )
         Call halo_swapper%fill( halo_width, Lbound( grid_with_halo ), w, grid_with_halo, error )
         Call FD_operator%apply( Lbound( grid_with_halo ), Lbound( r2 ), Lbound( r2 ), Ubound( r2 ), &
              grid_with_halo, r2  )
         s = method%contract( comms, w, w  )
         t = method%contract( comms, y, r2 )       
         z      = Abs(s - t)
         epsa   = (s + eps) * eps**0.33333
         If (z > epsa) Then
            istop = 6
            Exit solver_block
         End If
         Arnorml = Sqrt(s);
      Else
         Call halo_swapper%fill( halo_width, Lbound( grid_with_halo ), y, grid_with_halo, error )
         Call FD_operator%apply( Lbound( grid_with_halo ), Lbound( w ), Lbound( w ), Ubound( w ), &
              grid_with_halo, w  )
         Arnorml = Sqrt( method%contract( comms, w, w  ) )
      End If

      !-------------------------------------------------------------------
      ! Initialize other quantities.
      !-------------------------------------------------------------------
      oldb   = zero
      beta   = beta1
      dbar   = zero
      epsln  = zero
      qrnorm = beta1
      phibar = beta1
      rhs1   = beta1
      rhs2   = zero
      tnorm2 = zero
      ynorm2 = zero
      cs     = - one
      sn     = zero
      w  = zero
      w2 = zero
      r2 = r1

      ! Stop gfortran complaining about might be using uninitialised
      gmin =   Huge( gmin )
      gmax = - Huge( gmax )

      !===================================================================
      ! Main iteration loop.
      !===================================================================
      iteration_loop: Do
         itn = itn + 1               ! k = itn = 1 first time through

         !----------------------------------------------------------------
         ! Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
         ! The general iteration is similar to the case k = 1 with v0 = 0:
         !
         !   p1      = Operator * v1  -  beta1 * v0,
         !   alpha1  = v1'p1,
         !   q2      = p2  -  alpha1 * v1,
         !   beta2^2 = q2'q2,
         !   v2      = (1/beta2) q2.
         !
         ! Again, y = betak P vk,  where  P = C**(-1).
         ! .... more description needed.
         !----------------------------------------------------------------
         s      = one / beta            ! Normalize previous vector (in y).
         v      = s*y                   ! v = vk if P = I

         ! NEED TO FIX ARGUMENTS - especially lb of grid
         Call halo_swapper%fill( halo_width, Lbound( grid_with_halo ), v, grid_with_halo, error )
         Call FD_operator%apply( Lbound( grid_with_halo ), Lbound( y ), Lbound( y ), Ubound( y ), &
              grid_with_halo, y  )
         y      = y - shift*v           ! call daxpy ( n, (- shift), v, 1, y, 1 )
         If (itn >= 2) Then
            y   = y - (beta/oldb)*r1    ! call daxpy ( n, (- beta/oldb), r1, 1, y, 1 )
         End If

         alfa = method%contract( comms, v, y  )
         y      = y - (alfa/beta)*r2    ! call daxpy ( n, (- alfa/beta), r2, 1, y, 1 )
         r1     = r2
         r2     = y
         If ( precon ) Call Msolve( lb, ub, r2, y )

         oldb   = beta                  ! oldb = betak
         beta = method%contract( comms, r2, y  )
         If (beta < zero) Then
            istop = 6
            Exit solver_block
         End If

         beta   = Sqrt( beta )          ! beta = betak+1
         tnorm2 = tnorm2 + alfa**2 + oldb**2 + beta**2

         If (itn == 1) Then                   ! Initialize a few things.
            If (beta/beta1 <= ten*eps) Then   ! beta2 = 0 or ~ 0.
               istop = -1                     ! Terminate later.
            End If
            !tnorm2 = alfa**2
            gmax   = Abs( alfa )              ! alpha1
            gmin   = gmax                     ! alpha1
         End If

         ! Apply previous rotation Qk-1 to get
         !   [deltak epslnk+1] = [cs  sn][dbark    0   ]
         !   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

         oldeps = epsln
         delta  = cs * dbar  +  sn * alfa ! delta1 = 0         deltak
         gbar   = sn * dbar  -  cs * alfa ! gbar 1 = alfa1     gbar k
         epsln  =               sn * beta ! epsln2 = 0         epslnk+1
         dbar   =            -  cs * beta ! dbar 2 = beta2     dbar k+1

         ! Compute the next plane rotation Qk

         gamma  = Sqrt( gbar**2 + beta**2 )   ! gammak
         cs     = gbar / gamma                ! ck
         sn     = beta / gamma                ! sk
         phi    = cs * phibar                 ! phik
         phibar = sn * phibar                 ! phibark+1

         ! Update  x.

         denom = one/gamma

         w1 = w2
         w2 = w
         w = ( v - oldeps * w1 - delta * w2 ) * denom
         x = x + phi * w

         ! Go round again.

         gmax   = Max( gmax, gamma )
         gmin   = Min( gmin, gamma )
         z      = rhs1 / gamma
         ynorm2 = z**2  +  ynorm2
         rhs1   = rhs2  -  delta * z
         rhs2   =       -  epsln * z

         ! Estimate various norms and test for convergence.

         Anorm  = Sqrt( tnorm2 )
         ynorm  = Sqrt( ynorm2 )
         epsa   = Anorm * eps
         epsx   = Anorm * ynorm * eps
         epsr   = Anorm * ynorm * rtol
         diag   = gbar
         If (Abs(diag) < Tiny(diag)) diag = epsa

         qrnorm = phibar
         rnorml = rnorm
         rnorm  = qrnorm
         rootl       = Sqrt( gbar**2 +dbar**2  )  ! norm([gbar; dbar]);
         Arnorml     = rnorml*rootl               ! ||A r_{k-1} ||
         relArnorml  = rootl  /  Anorm;           ! ||Ar|| / (||A|| ||r||)     
         !relArnorml = Arnorml / Anorm;           ! ||Ar|| / ||A|| 

         ! Estimate  cond(A).
         ! In this version we look at the diagonals of  R  in the
         ! factorization of the lower Hessenberg matrix,  Q * H = R,
         ! where H is the tridiagonal matrix from Lanczos with one
         ! extra row, beta(k+1) e_k^T.

         Acond  = gmax / gmin

         ! See if any of the stopping criteria are satisfied.
         ! In rare cases, istop is already -1 from above (Abar = const*I).

         If (istop == 0) Then
            If (itn    >= itnlim    ) istop = 5
            If (Acond  >= 0.1d+0/eps) istop = 4
            If (epsx   >= beta1     ) istop = 3
            If (qrnorm <= epsx  .Or.  relArnorml <= epsx) istop = 2
            If (qrnorm <= epsr  .Or.  relArnorml <= epsr) istop = 1
         End If

         If (istop /= 0) Exit

      End Do iteration_loop
      !===================================================================
      ! End of iteration loop.
      !===================================================================

    End Block solver_block

    ! Display final status.

    Arnorm = Arnorml

    istop_message = msg(istop)

  End Subroutine MINRES

End Module equation_solver_minres_module
