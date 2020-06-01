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
       lb, ub, Msolve, b, rtol,  precon, &
       x, istop, istop_message, itn, rnorm )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use halo_setter_base_module, Only : halo_setter_base_class
    Use FD_template_module     , Only : FD_template

    Implicit None
    
    Class( equation_solver_conjugate_gradient )           , Intent( InOut ) :: method
    Integer,  Dimension( 1:3 )                            , Intent( In    ) :: lb( 1:3 )
    Integer,  Dimension( 1:3 )                            , Intent( In    ) :: ub( 1:3 )
    Real( wp ) , Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ) :: b
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
    Real( wp ) :: mu
    
    Integer :: halo_width
    Integer :: FD_order
    Integer :: iteration
    Integer :: np
    Integer :: error

    Call method%comms%get_size( np )
    
    Allocate( r, p, w, Mold = b )

    rnorm = Huge( rnorm )
    
    x = 0.0_wp
    
    !IJB
    ! Note order is always even, so no worries about splitting it in 2
    FD_order  = method%FD_operator%get_order()

!!$    halo_width = FD_order / 2
    ! Looks like a bug in get order! Returns twice what expected ...
    halo_width = FD_order
    Call method%halo_swapper%allocate( lb, ub, halo_width, grid_with_halo )

    Call method%halo_swapper%fill( halo_width, Lbound( grid_with_halo ), x, grid_with_halo, error )
    Call method%FD_operator%apply( Lbound( grid_with_halo ), Lbound( w ), Lbound( w ), Ubound( w ), &
         grid_with_halo, w  )
    r = b - w
    mu = Sum( r )
    ! Need n_grid in here really!!!! So can for av value of r
    Call method%comms%reduce( mu )
    mu = mu / ( np * Size( r ) )
    r = r - mu
    p = r
    r_dot_r_old = method%contract( method%comms, r, r )
    rnorm = Sqrt( r_dot_r_old )
    If( rnorm > rtol ) Then
       Do iteration = 1, method%max_iter
          Call method%halo_swapper%fill( halo_width, Lbound( grid_with_halo ), p, grid_with_halo, error )
          Call method%FD_operator%apply( Lbound( grid_with_halo ), Lbound( w ), Lbound( w  ), Ubound( w  ), &
               grid_with_halo, w  )
          p_dot_w = method%contract( method%comms, p, w )
          alpha = r_dot_r_old / p_dot_w
          x = x + alpha * p
          r = r - alpha * w
    mu = Sum( r )
    mu = mu / Size( r )
    Call method%comms%reduce( mu )
    mu = mu / ( np * Size( r ) )
    r = r - mu
          r_dot_r = method%contract( method%comms, r, r )
          ! NEED BETTER CONVERGENCE CRITERION!!!
          rnorm = Sqrt( r_dot_r )
          If( rnorm < rtol ) Exit
          beta = r_dot_r / r_dot_r_old
          p = r + beta * p
          r_dot_r_old = r_dot_r
       End Do
    Else
       iteration = 0
    End If
       
    itn = iteration

    If( iteration <= method%max_iter ) Then
       istop = 1
       istop_message = "CG worked"
    Else
       istop = -1
       istop_message = "CG FAILED!!!!!"
    End If
    
  End Subroutine cg

End Module equation_solver_conjugate_gradient_module
