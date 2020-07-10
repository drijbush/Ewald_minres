Module equation_solver_weighted_jacobi_module

  Use equation_solver_precon_base_class_module, Only : equation_solver_precon_base_class

  Use constants, Only : wp

  Implicit None

  Type, Public, Extends( equation_solver_precon_base_class ) :: equation_solver_weighted_jacobi
!!$     Real( wp ) :: weight = 0.666666666666666666666666666666666666666666666666666666666666_wp
     Real( wp ) :: weight = 1.0_wp
   Contains
     Procedure :: solve => weighted_jacobi
  End type equation_solver_weighted_jacobi

Contains

  Subroutine weighted_jacobi( method, &
       lb, ub, b, rtol, &
       x, istop, istop_message, itn, rnorm )

    Use constants, Only : wp
    Use halo_setter_base_module, Only : halo_setter_base_class
    Use FD_template_module,      Only : FD_template

    Implicit None

    Class( equation_solver_weighted_jacobi ),               Intent( InOut ) :: method
    Integer,  Dimension( 1:3 ),                             Intent( In    ) :: lb( 1:3 )
    Integer,  Dimension( 1:3 ),                             Intent( In    ) :: ub( 1:3 )
    Real( wp ),  Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ) :: b
    Real( wp ),                                             Intent( In    ) :: rtol
    Real( wp ),  Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: x
    Integer,                                                Intent(   Out ) :: istop
    Character( Len = * ),                                   Intent(   Out ) :: istop_message
    Real( wp ),                                             Intent(   Out ) :: rnorm
    Integer,                                                Intent(   Out ) :: itn

    Real( wp ), Dimension( :, :, : ), Allocatable :: soln_in

    Real( wp ) :: diff, max_diff

    Integer :: halo_width
    Integer :: FD_order
    Integer :: iteration
    Integer :: error
    Integer :: i1, i2, i3
    Integer :: me

    Call method%comms%get_rank( me )

    ! No solution so norm of residual huge ...
    rnorm = Huge( rnorm )
    max_diff = rnorm

    ! Initial guess at result - assume diagonal matrix
    x = b * method%FD_operator%get_diag_inv()

    !IJB
    ! Note order is always even, so no worries about splitting it in 2
    FD_order  = method%FD_operator%get_order()

!!$    halo_width = FD_order / 2
    ! Looks like a bug in get order! Returns twice what expected ...
    halo_width = FD_order
    Call method%halo_swapper%allocate( lb, ub, halo_width, soln_in  )

    jacobi_iteration: Do iteration = 1, method%max_iter
       Call method%halo_swapper%fill( halo_width, Lbound( soln_in ), x, soln_in, error )
       Call method%FD_operator%jacobi_sweep( Lbound( b ), Lbound( soln_in ), Lbound( b ), Ubound( b ), &
            method%weight, b, soln_in, x )
       ! NEED BETTER CONV CHECK
       max_diff = - 1.0_wp
       Do i3 = Lbound( b, Dim = 3 ), Ubound( b, Dim = 3 )
          Do i2 = Lbound( b, Dim = 2 ), Ubound( b, Dim = 2 )
             Do i1 = Lbound( b, Dim = 1 ), Ubound( b, Dim = 1 )
                diff = soln_in( i1, i2, i3 ) - x( i1, i2, i3 )
                max_diff = Max( max_diff, Abs( diff ) )
             End Do
          End Do
       End Do
!!$       Call method%comms%max( max_diff )
!!$       If( me == 0 ) Then
!!$          Write( *, * ) iteration, max_diff, x( 10, 10, 10 ), b( 10, 10, 10 )
!!$       End If
       If( max_diff < rtol ) Then
          Exit
       End If
    End Do jacobi_iteration

    itn = iteration

    If( iteration <= method%max_iter ) Then
       istop = 1
       istop_message = "Weighted Jacobi worked"
    Else
       istop = -1
       istop_message = "Weighted Jacobi FAILED!!!!!"
    End If

    rnorm = max_diff

  End Subroutine weighted_jacobi

End Module equation_solver_weighted_jacobi_module
