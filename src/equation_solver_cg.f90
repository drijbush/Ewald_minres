Module equation_solver_conjugate_gradient_module

  Use constants, Only : wp
  Use equation_solver_precon_base_class_module, Only : equation_solver_precon_base_class

  Implicit None

  Type, Public, Extends( equation_solver_precon_base_class ) :: equation_solver_conjugate_gradient
   Contains
     Procedure :: solve => cg
  End type equation_solver_conjugate_gradient

Contains

  Subroutine cg( method, &
       lb, ub, b, rtol, &
       x, istop, istop_message, itn, rnorm )

    Use halo_setter_base_module, Only : halo_setter_base_class
    Use FD_template_module,      Only : FD_template

    Implicit None

    Class( equation_solver_conjugate_gradient ),            Intent( InOut ) :: method
    Integer,  Dimension( 1:3 ),                             Intent( In    ) :: lb( 1:3 )
    Integer,  Dimension( 1:3 ),                             Intent( In    ) :: ub( 1:3 )
    Real( wp ),  Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ) :: b
    Real( wp ),                                             Intent( In    ) :: rtol
    Real( wp ),  Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: x
    Integer,                                                Intent(   Out ) :: istop
    Character( Len = * ),                                   Intent(   Out ) :: istop_message
    Real( wp ),                                             Intent(   Out ) :: rnorm
    Integer,                                                Intent(   Out ) :: itn

    Real( wp ), Dimension( :, :, : ), Allocatable :: r, p, w, z
    Real( wp ), Dimension( :, :, : ), Allocatable :: grid_with_halo

    Real( wp ) :: alpha, beta
    Real( wp ) :: z_dot_r_old, z_dot_r, r_dot_r, p_dot_w

    Integer :: halo_width
    Integer :: FD_order
    Integer :: iteration
    Integer :: error

    Allocate( r, p, w, z, Mold = b )

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
    If( Allocated( method%precon ) ) Then
       Call method%precon%solve( lb, ub, r, rtol, z, istop, istop_message, itn, rnorm )
    Else
       z = r
    End If
    p = z
    z_dot_r_old = method%contract( method%comms, z, r )
    r_dot_r = method%contract( method%comms, r, r )
    rnorm = Sqrt( r_dot_r )
    If( rnorm > rtol ) Then
       Do iteration = 1, method%max_iter
          Call method%halo_swapper%fill( halo_width, Lbound( grid_with_halo ), p, grid_with_halo, error )
          Call method%FD_operator%apply( Lbound( grid_with_halo ), Lbound( w ), Lbound( w  ), Ubound( w  ), &
               grid_with_halo, w  )
          p_dot_w = method%contract( method%comms, p, w )
          alpha = z_dot_r_old / p_dot_w
          x = x + alpha * p
          r = r - alpha * w
          r_dot_r = method%contract( method%comms, r, r )
          ! NEED BETTER CONVERGENCE CRITERION!!!
          rnorm = Sqrt( r_dot_r )
          If( rnorm < rtol ) Exit
          If( Allocated( method%precon ) ) Then
             Call method%precon%solve( lb, ub, r, rtol, z, istop, istop_message, itn, rnorm )
          Else
             z = r
          End If
          z_dot_r = method%contract( method%comms, z, r )
          beta = z_dot_r / z_dot_r_old
          p = z + beta * p
          z_dot_r_old = z_dot_r
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
