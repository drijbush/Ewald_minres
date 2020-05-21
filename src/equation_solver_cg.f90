Module equation_solver_conjugate_gradient_module

  Use equation_solver_base_class_module, Only : equation_solver_base_class
  
  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Type, Public, Extends( equation_solver_base_class ) :: equation_solver_conjugate_gradient
   Contains
     Procedure :: solve => cg
  End type equation_solver_conjugate_gradient

Contains

  Subroutine cg( method, &
       lb, ub, FD_operator, comms, halo_swapper, Msolve, b, itnlim, rtol,  precon, &
       x, istop, istop_message, itn, rnorm )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use comms_base_class_module, Only : comms_base_class
    Use halo_setter_base_module, Only : halo_setter_base_class
    Use FD_template_module     , Only : FD_template

    Implicit None
    
    Class( equation_solver_conjugate_gradient )           , Intent( In    ) :: method
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

    Real( wp ), Dimension( :, :, : ), Allocatable :: r, p, w
    Real( wp ), Dimension( :, :, : ), Allocatable :: grid_with_halo

    Real( wp ) :: alpha, beta
    Real( wp ) :: r_dot_r_old, r_dot_r, p_dot_w

    Integer :: halo_width
    Integer :: FD_order
    Integer :: iteration
    Integer :: error
    
    Allocate( r, p, w, Mold = b )

    rnorm = Huge( rnorm )
    
    x = 0.0_wp
    
    !IJB
    ! Note order is always even, so no worries about splitting it in 2
    FD_order  = FD_operator%get_order()

!!$    halo_width = FD_order / 2
    ! Looks like a bug in get order! Returns twice what expected ...
    halo_width = FD_order
    Call halo_swapper%allocate( lb, ub, halo_width, grid_with_halo )

    Call halo_swapper%fill( halo_width, Lbound( grid_with_halo ), x, grid_with_halo, error )
    Call FD_operator%apply( Lbound( grid_with_halo ), Lbound( w ), Lbound( w ), Ubound( w ), &
         grid_with_halo, w  )
    r = b - w
    p = r
    r_dot_r_old = method%contract( comms, r, r )
    Do iteration = 1, itnlim
!!$    Do iteration = 1, 30
       Call halo_swapper%fill( halo_width, Lbound( grid_with_halo ), p, grid_with_halo, error )
       Call FD_operator%apply( Lbound( grid_with_halo ), Lbound( w ), Lbound( w  ), Ubound( w  ), &
            grid_with_halo, w  )
       p_dot_w = method%contract( comms, p, w )
       alpha = r_dot_r_old / p_dot_w
       x = x + alpha * p
       r = r - alpha * w
       r_dot_r = method%contract( comms, r, r )
       ! NEED BETTER CONVERGENCE CRITERION!!!
       rnorm = Sqrt( r_dot_r )
       If( rnorm < rtol ) Exit
       beta = r_dot_r / r_dot_r_old
       p = r + beta * p
       r_dot_r_old = r_dot_r
    End Do

    itn = iteration

    If( iteration <= itnlim ) Then
       istop = 1
       istop_message = "CG worked"
    Else
       istop = -1
       istop_message = "CG FAILED!!!!!"
    End If
    
  End Subroutine cg

End Module equation_solver_conjugate_gradient_module
