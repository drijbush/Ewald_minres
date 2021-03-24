Module equation_solver_hypre_pfmg_module

  Use, Intrinsic :: iso_C_binding, Only : c_ptr

  Use constants                               , Only : wp
  Use equation_solver_precon_base_class_module, Only : equation_solver_precon_base_class
  Use FD_template_module                      , Only : FD_template

  Implicit None

  Private
  
  Type, Public, Extends( equation_solver_precon_base_class ) :: equation_solver_hypre_pfmg
     Integer, Dimension( 1:3 )                     :: n                     ! Number of grid point
     Real( wp )                                    :: V                     ! Volume
     Real( wp )                                    :: iter_tol = 1.0e-06_wp ! Tolerance for iterations around order 2 solver
     Logical                                       :: orthog_stencil        ! If the stencil only has elements along the axes can save some message passing
     Real( wp ), Dimension( :, :, : ), Allocatable :: stencil_2
     Real( wp )                                    :: stencil_diag_2
     Type( c_ptr )                                 :: hypre_objects
   Contains
     Procedure, Public :: solve => pfmg
     Procedure, Public :: pfmg_init
     Procedure, Public :: set_iter_tol
!!$     Final             :: pfmg_destroy
  End Type equation_solver_hypre_pfmg

  ! Experiment partial mixing in "SCF" procedure
  Real( wp ), Parameter :: mix = 1.0_wp

  ! Set to true to get report during iteration. Otherwise no printing
  Logical :: report = .True.

  ! If true the returned residual is consistent with the whole calculation.
  ! However calculation of the residual requies significant extra calculation
  ! which is not actually needed in most situations. Thus set to false to avoid this overhead.
  Logical :: consistent_residual = .False.
  
  ! Interfaces for C functions which talk to hpyre
  
  Interface

     Function ssp_hypre_struct_setup( comm, n, lb, ub, n_stencil, stencil_elements, stencil_values ) &
          Result( hypre_objects ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr, c_int, c_double
       Implicit None
       Type   ( c_ptr    )                                                     :: hypre_objects
       Integer( c_int    )                               , Intent( In ), Value :: comm
       Integer( c_int    ), Dimension( 1:3 )             , Intent( In )        :: n
       Integer( c_int    ), Dimension( 1:3 )             , Intent( In )        :: lb
       Integer( c_int    ), Dimension( 1:3 )             , Intent( In )        :: ub
       Integer( c_int    )                               , Intent( In ), Value :: n_stencil
       Integer( c_int    ), Dimension( 1:3, 1:n_stencil ), Intent( In )        :: stencil_elements
       Real   ( c_double ), Dimension(      1:n_stencil ), Intent( In )        :: stencil_values
     End Function ssp_hypre_struct_setup

     Subroutine ssp_hypre_struct_pfmg_solve( hypre_objects, n1, n2, n3, rtol, b, x, n_iter, residual, info ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr, c_int, c_double
       Implicit None
       Type   ( c_ptr    )                               , Intent( In    ), Value :: hypre_objects
       Integer( c_int    )                               , Intent( In    ), Value :: n1
       Integer( c_int    )                               , Intent( In    ), Value :: n2
       Integer( c_int    )                               , Intent( In    ), Value :: n3
       Real   ( c_double )                               , Intent( In    ), Value :: rtol
       Real   ( c_double ), Dimension( 1:n1, 1:n2, 1:n3 ), Intent( In    )        :: b
       Real   ( c_double ), Dimension( 1:n1, 1:n2, 1:n3 ), Intent( InOut )        :: x
       Integer( c_int    )                               , Intent(   Out )        :: n_iter
       Real   ( c_double )                               , Intent(   Out )        :: residual
       Integer( c_int    )                               , Intent(   Out )        :: info
     End Subroutine ssp_hypre_struct_pfmg_solve

     Subroutine ssp_hypre_struct_free( hypre_objects ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr
       Implicit None
       Type( c_ptr ), Intent( In ), Value :: hypre_objects
     End Subroutine ssp_hypre_struct_free
     
     Function ssp_hypre_semi_struct_setup( comm, n, lb, ub, n_stencil, stencil_elements, stencil_values ) &
          Result( hypre_objects ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr, c_int, c_double
       Implicit None
       Type   ( c_ptr    )                                                     :: hypre_objects
       Integer( c_int    )                               , Intent( In ), Value :: comm
       Integer( c_int    ), Dimension( 1:3 )             , Intent( In )        :: n
       Integer( c_int    ), Dimension( 1:3 )             , Intent( In )        :: lb
       Integer( c_int    ), Dimension( 1:3 )             , Intent( In )        :: ub
       Integer( c_int    )                               , Intent( In ), Value :: n_stencil
       Integer( c_int    ), Dimension( 1:3, 1:n_stencil ), Intent( In )        :: stencil_elements
       Real   ( c_double ), Dimension(      1:n_stencil ), Intent( In )        :: stencil_values
     End Function ssp_hypre_semi_struct_setup

     Subroutine ssp_hypre_semi_struct_pfmg_solve( hypre_objects, n1, n2, n3, b, x, n_iter, residual, info ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr, c_int, c_double
       Implicit None
       Type   ( c_ptr    )                               , Intent( In    ), Value :: hypre_objects
       Integer( c_int    )                               , Intent( In    ), Value :: n1
       Integer( c_int    )                               , Intent( In    ), Value :: n2
       Integer( c_int    )                               , Intent( In    ), Value :: n3
       Real   ( c_double ), Dimension( 1:n1, 1:n2, 1:n3 ), Intent( In    )        :: b
       Real   ( c_double ), Dimension( 1:n1, 1:n2, 1:n3 ), Intent( InOut )        :: x
       Integer( c_int    )                               , Intent(   Out )        :: n_iter
       Real   ( c_double )                               , Intent(   Out )        :: residual
       Integer( c_int    )                               , Intent(   Out )        :: info
     End Subroutine ssp_hypre_semi_struct_pfmg_solve

     Subroutine ssp_hypre_semi_struct_free( hypre_objects ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr
       Implicit None
       Type( c_ptr ), Intent( In ), Value :: hypre_objects
     End Subroutine ssp_hypre_semi_struct_free
     
     Function ssp_hypre_ij_setup( comm, n, lb, ub, n_stencil, stencil_elements, stencil_values ) &
          Result( hypre_objects ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr, c_int, c_double
       Implicit None
       Type   ( c_ptr    )                                                     :: hypre_objects
       Integer( c_int    )                               , Intent( In ), Value :: comm
       Integer( c_int    ), Dimension( 1:3 )             , Intent( In )        :: n
       Integer( c_int    ), Dimension( 1:3 )             , Intent( In )        :: lb
       Integer( c_int    ), Dimension( 1:3 )             , Intent( In )        :: ub
       Integer( c_int    )                               , Intent( In ), Value :: n_stencil
       Integer( c_int    ), Dimension( 1:3, 1:n_stencil ), Intent( In )        :: stencil_elements
       Real   ( c_double ), Dimension(      1:n_stencil ), Intent( In )        :: stencil_values
     End Function ssp_hypre_ij_setup

     Subroutine ssp_hypre_ij_solve( hypre_objects, n1, n2, n3, b, x, n_iter, residual, info ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr, c_int, c_double
       Implicit None
       Type   ( c_ptr    )                               , Intent( In    ), Value :: hypre_objects
       Integer( c_int    )                               , Intent( In    ), Value :: n1
       Integer( c_int    )                               , Intent( In    ), Value :: n2
       Integer( c_int    )                               , Intent( In    ), Value :: n3
       Real   ( c_double ), Dimension( 1:n1, 1:n2, 1:n3 ), Intent( In    )        :: b
       Real   ( c_double ), Dimension( 1:n1, 1:n2, 1:n3 ), Intent( InOut )        :: x
       Integer( c_int    )                               , Intent(   Out )        :: n_iter
       Real   ( c_double )                               , Intent(   Out )        :: residual
       Integer( c_int    )                               , Intent(   Out )        :: info
     End Subroutine ssp_hypre_ij_solve

     Subroutine ssp_hypre_ij_free( hypre_objects ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr
       Implicit None
       Type( c_ptr ), Intent( In ), Value :: hypre_objects
     End Subroutine ssp_hypre_ij_free
     
  End Interface

  ! HACK to get rid of compiler warnings about unused arguments
  Class( FD_template ), Allocatable :: last_operator
  
Contains

  Subroutine pfmg_init( method, comms, n, lb, ub, FD_operator, V )

    Use constants              , Only : wp
    Use comms_base_class_module, Only : comms_base_class
    Use FD_template_module     , Only : FD_template

    Use FD_Laplacian_3d_module, Only : FD_Laplacian_3d
    
    Implicit None

    Real( wp ), Parameter :: stencil_tol = 1e-15_wp

    Class( equation_solver_hypre_pfmg ), Intent( InOut ) :: method
    Class( comms_base_class       )    , Intent( In    ) :: comms
    Integer, Dimension( 1:3 )          , Intent( In    ) :: n
    Integer, Dimension( 1:3 )          , Intent( In    ) :: lb
    Integer, Dimension( 1:3 )          , Intent( In    ) :: ub
    Class( FD_template )               , Intent( In    ) :: FD_operator
    Real( wp )                         , Intent( In    ) :: V

    Real( wp ), Dimension( :, :, : ), Allocatable :: stencil

    Real( wp ), Dimension( : ), Allocatable :: stencil_values
    
    Real( wp ), Dimension( 1:3, 1:3 ) :: dGrid_vecs

    Integer, Dimension( :, : ), Allocatable :: stencil_elements

    Integer :: comm
    Integer :: n_stencil
    Integer :: i1, i2, i3
    Integer :: i

    Type( FD_Laplacian_3d ) :: FD_order_2

    !HACK to get rid of silly compiler warnings
    last_operator = FD_operator

    method%n = n
    method%V = V
    
    Call comms%get_comm( comm )

    ! get higher accuracy by the iteration described by Matt and Phil - multigrid paper, need to look up ref
    dGrid_vecs( :, 1 ) = method%FD_operator%get_dir_vec( 1 )
    dGrid_vecs( :, 2 ) = method%FD_operator%get_dir_vec( 2 )
    dGrid_vecs( :, 3 ) = method%FD_operator%get_dir_vec( 3 )
    Call FD_order_2%init( 2, dGrid_vecs )
    Call FD_order_2%get_stencil( stencil )
    method%stencil_2      = stencil
    method%stencil_diag_2 = stencil( 0, 0, 0 )
    
    ! Build up stencil list
    Allocate( stencil_values( 1:0 ) )
    ! Hack, but shouldn't waste too much memory
    Allocate( stencil_elements( 1:3, 1:Size( stencil ) ) )
    n_stencil = 0

    Do i3 = Lbound( stencil, Dim = 3 ), Ubound( stencil, Dim = 3 )
       Do i2 = Lbound( stencil, Dim = 2 ), Ubound( stencil, Dim = 2 )
          Do i1 = Lbound( stencil, Dim = 1 ), Ubound( stencil, Dim = 1 )
             If( Abs( stencil( i1, i2, i3 ) ) >= stencil_tol ) Then
                n_stencil = n_stencil + 1
                stencil_values = [ stencil_values, stencil( i1, i2, i3 ) ]
                stencil_elements( :, n_stencil ) = [ i1, i2, i3 ]
             End If
          End Do
       End Do
    End Do

    ! Work out if we have an orthogonal stencil
    method%orthog_stencil = .True.
    Do i = 1, n_stencil
       method%orthog_stencil = method%orthog_stencil .And. ( Count( stencil_elements( :, i ) /= 0 ) == 1 )
       If( .Not. method%orthog_stencil ) Exit
    End Do

    method%hypre_objects = ssp_hypre_struct_setup( comm, n, lb, ub, n_stencil, stencil_elements, stencil_values )

  End Subroutine pfmg_init

  Subroutine pfmg( method, &
       lb, ub, b, rtol, &
       x, istop, istop_message, itn, rnorm )
    
    Use constants, Only : wp, pi
    Use halo_setter_base_module, Only : halo_setter_base_class
    Use FD_template_module,      Only : FD_template

    Implicit None

    Class( equation_solver_hypre_pfmg ),                    Intent( InOut ) :: method
    Integer,  Dimension( 1:3 ),                             Intent( In    ) :: lb( 1:3 )
    Integer,  Dimension( 1:3 ),                             Intent( In    ) :: ub( 1:3 )
    Real( wp ),  Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ) :: b
    Real( wp ),                                             Intent( In    ) :: rtol
    Real( wp ),  Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: x
    Integer,                                                Intent(   Out ) :: istop
    Character( Len = * ),                                   Intent(   Out ) :: istop_message
    Real( wp ),                                             Intent(   Out ) :: rnorm
    Integer,                                                Intent(   Out ) :: itn

    Integer, Parameter :: max_iter = 50
    
    Real( wp ), Dimension( :, :, : ), Allocatable :: x_with_halo, r, e, w

    Real( wp ) :: r2, x2, x2_old, dx2, dx_max
    Real( wp ) :: tol
    Real( wp ) :: Energy, Energy_old, dE
    Real( wp ) :: lambda

    Integer, Dimension( 1:3 ) :: n

    Integer :: FD_order, halo_width
    Integer :: rank
    Integer :: lb1, lb2, lb3
    Integer :: ub1, ub2, ub3
    Integer :: error
    Integer :: outer_iter, solver_iter

    Logical :: do_io
    Logical :: do_corners
    
    n = Ubound( b ) - Lbound( b ) + 1

    Call method%comms%get_rank( rank )
    do_io = rank == 0 .And. report
    
    ! Use tol to think about varying tolerance at different stages
    tol = rtol

    ! Get the initial Order 2 solution
    x = b * ( 1.0_wp /  method%stencil_diag_2 )
    Call ssp_hypre_struct_pfmg_solve( method%hypre_objects, n( 1 ), n( 2 ), n( 3 ), tol, b, x, solver_iter, rnorm, istop )
    If( report ) Then
       x2_old = Sqrt( method%contract( method%comms, x, x  ) )
    Else
       x2_old = 0.0_wp
    End If

    ! Set up the data structures for the higher order solution
    ! gfortran 9 has a bug with Mold that doesn't propagate the bounds of the arrays
    ! properly, sigh ... So do it "manually"
!!$    Allocate( r, e, w, Mold = b )
    lb1 = Lbound( b, Dim = 1 )
    lb2 = Lbound( b, Dim = 2 )
    lb3 = Lbound( b, Dim = 3 )
    ub1 = Ubound( b, Dim = 1 )
    ub2 = Ubound( b, Dim = 2 )
    ub3 = Ubound( b, Dim = 3 )
    Allocate( r( lb1:ub1, lb2:ub2, lb3:ub3 ) )
    Allocate( e( lb1:ub1, lb2:ub2, lb3:ub3 ) )
    Allocate( w( lb1:ub1, lb2:ub2, lb3:ub3 ) )
    FD_order   = method%FD_operator%get_order()
    halo_width = FD_order
    Call method%halo_swapper%allocate( lb, ub, halo_width, x_with_halo )

    ! Initial Energy
    Energy = method%contract( method%comms, x, b  ) * method%V / Product( method%n )
    ! And scale with the appropriate constants
    Energy = - 0.5_wp * Energy / ( 4.0_wp * pi )
    Energy_old = Energy

    If( do_io ) Then
       Write( *, * ) 'Initial energy: ', Energy
       Write( *, * ) 'Mix = ', mix
       Write( *, * )
       Write( *, '( a7, 2a5 )' ) '  FD   ', 'Outer', 'FD 2 '
       Write( *, '( a7, 2a5, t20, 4a13, 2a16 )' ) &
            'order ', 'iter ', 'iter ', &
            'FD 2 norm    ', 'Full norm    ', 'dSoln Norm   ', 'dSoln Max    ', 'Energy       ', 'dE           '
    End If

    ! If orthog stencil can save some message passing
    If( method%orthog_stencil ) Then
       do_corners = method%halo_swapper%do_corners()
       Call method%halo_swapper%set_corners( .False. )
    End If
    
    iteration_loop: Do outer_iter = 1, max_iter

       ! Get the residual for the higher order operator
       Call method%halo_swapper%fill( halo_width, Lbound( x_with_halo ), x, x_with_halo, error )
       Call method%FD_operator%apply( Lbound( x_with_halo ), Lbound( w ), Lbound( w ), Ubound( w ), &
            x_with_halo, w  )
       r = b - w
       ! The residual should not contain any null space component as b can not contain any, and the
       ! action of the matrix on the solution projects out the null space, by definition
       ! Hence we can just call the solver
       e = r * ( 1.0_wp /  method%stencil_diag_2 ) 
       Call ssp_hypre_struct_pfmg_solve( method%hypre_objects, n( 1 ), n( 2 ), n( 3 ), tol, r, e, solver_iter, rnorm, istop )
       ! Update the solution - written this way to allow mixing c.f. SCF solvers
       x = x + mix * e
       ! New energy
       ! Note we are essentially do the quadrature by hand as a little painful to interface to the quadrature type from here
       Energy = method%contract( method%comms, x, b  ) * method%V / Product( method%n )
       ! And scale with the appropriate constants
       Energy = - 0.5_wp * Energy / ( 4.0_wp * pi )
       dE = Energy - Energy_old
       ! If reqd tell the world things
       If( report ) Then
          r2 = Sqrt( method%contract( method%comms, r, r  ) )
          x2 = Sqrt( method%contract( method%comms, x, x  ) )
          dx2 = Abs( x2 - x2_old )
          dx_max = Maxval( Abs( e ) )
          Call method%comms%max( dx_max )
          If( do_io ) Then
             Write( *, '( 2x, 3( i4, 1x ), t20, 4( g12.6, 1x ), 2( g16.10, 1x ) )' ) &
                  FD_order * 2, outer_iter, solver_iter, rnorm, r2, dx2, dx_max, Energy, Abs( Energy - Energy_old )
          End If
          x2_old = x2
       End If
       If( Abs( dE ) < method%iter_tol ) Then
          Exit iteration_loop
       End If
       ! After the solve e may now be contaminated with components of the Null space, so project null space out
       ! so allowing us to go around again
       Call project_out_null( method%n, method%comms, e, lambda )
       ! And around we go
       Energy_old = Energy
    End Do iteration_loop

    If( consistent_residual ) Then
       ! This is not strictly needed, but without it the reported residual is not consistent with the above calculation
       ! However the solution is consistent
       Call method%halo_swapper%fill( halo_width, Lbound( x_with_halo ), x, x_with_halo, error )
       Call method%FD_operator%apply( Lbound( x_with_halo ), Lbound( w ), Lbound( w ), Ubound( w ), &
            x_with_halo, w  )
       r = b - w
    End If
    r2 = Sqrt( method%contract( method%comms, r, r  ) )
    rnorm = r2

    ! Reset the method as we found it
    If( method%orthog_stencil ) Then
       Call method%halo_swapper%set_corners( do_corners )
    End If

    istop_message = "Who knows?"

    itn = outer_iter

  Contains

    Subroutine project_out_null( n_tot, comms, v, lambda_out )

      ! Project out the Null space from the vector v
      ! We know as the matrix is a circulant matrix that the null space is the
      ! vector ( 1, 1, 1, ..., 1, 1 ) represented by n below

      ! Use Gram-Schmidt for the projection
      
      Use constants              , Only : wp
      Use comms_base_class_module, Only : comms_base_class

      Integer   , Dimension( : )      , Intent( In    )           :: n_tot
      Class( comms_base_class )       , Intent( In    )           :: comms
      Real( wp ), Dimension( :, :, : ), Intent( InOut )           :: v
      Real( wp ),                       Intent(   Out ), Optional :: lambda_out

      Real( wp ), Parameter :: proj_tol = 1.0e-14_wp
      
      Real( wp ) :: lambda
      Real( wp ) :: v_dot_n, n_dot_n

      ! Calculate v.n - the dot product of the vector with the null space
      v_dot_n = kahan_Sum( v )
      Call comms%reduce( v_dot_n )

      ! Calculate n.n - the dot product of the null space with itself,
      ! which is simply the total size of the grid
      n_dot_n = Product( n_tot )

      ! Project out the null space
      lambda = - v_dot_n / n_dot_n
      If( Abs( lambda ) > proj_tol ) Then
         v = v + lambda
      End If

      If( Present( lambda_out ) ) Then
         lambda_out = lambda
      End If
      
    End Subroutine project_out_null

    Pure Function kahan_sum( v ) Result( r )

      Use constants              , Only : wp

      Real( wp ) :: r

      Real( wp ), Dimension( :, :, : ), Intent( In ) :: v

      Real( wp ) :: c, y, t

      Integer :: i1, i2, i3
      
      r = 0.0_wp
      !$omp parallel default( none ) shared( v, r ) private( c, y, t, i1, i2, i3 )
      c = 0.0_wp
      !$omp do collapse( 3 ) reduction( +:r )
      Do i3 = Lbound( v, Dim = 3 ), Ubound( v, Dim = 3 )
         Do i2 = Lbound( v, Dim = 2 ), Ubound( v, Dim = 2 )
            Do i1 = Lbound( v, Dim = 1 ), Ubound( v, Dim = 1 )
               y = v( i1, i2, i3 ) - c
               t = r + y
               c = ( t - r ) - y
               r = t
            End Do
         End Do
      End Do
      !$omp end do
      !$omp end parallel
      
    End Function kahan_sum

  End Subroutine pfmg

!!$  Subroutine pfmg_destroy( method )
!!$
!!$    Implicit None
!!$
!!$    Type( equation_solver_hypre_pfmg ), Intent( InOut ) :: method
!!$
!!$    Call  ssp_hypre_struct_free( method%hypre_objects )
!!$
!!$  End Subroutine pfmg_destroy
  
  Subroutine set_iter_tol( method, tol )

    Implicit None

    Class( equation_solver_hypre_pfmg ), Intent( InOut ) :: method
    Real( wp )                         , Intent( In    ) :: tol

    method%iter_tol = tol

  End Subroutine set_iter_tol
  
End Module equation_solver_hypre_pfmg_module
