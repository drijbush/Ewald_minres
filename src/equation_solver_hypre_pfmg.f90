Module equation_solver_hypre_pfmg_module

  Use, Intrinsic :: iso_C_binding, Only : c_ptr

  Use constants, Only : wp
  Use equation_solver_precon_base_class_module, Only : equation_solver_precon_base_class

  Implicit None

  Private
  
  Type, Public, Extends( equation_solver_precon_base_class ) :: equation_solver_hypre_pfmg
     Type( c_ptr ) :: hypre_objects
   Contains
     Procedure, Public :: solve => pfmg
     Procedure, Public :: pfmg_init
!!$     Final             :: pfmg_destroy
  End Type equation_solver_hypre_pfmg

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

     Subroutine ssp_hypre_struct_pfmg_solve( hypre_objects, n1, n2, n3, b, x, n_iter, residual, info ) Bind( C )
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
     End Subroutine ssp_hypre_struct_pfmg_solve

     Subroutine ssp_hypre_struct_free( hypre_objects ) Bind( C )
       Use, Intrinsic :: iso_C_binding, Only : c_ptr
       Implicit None
       Type( c_ptr ), Intent( In ), Value :: hypre_objects
     End Subroutine ssp_hypre_struct_free
     
  End Interface

Contains

  Subroutine pfmg_init( method, comms, n, lb, ub, FD_operator )

    Use comms_base_class_module,           Only : comms_base_class
    Use FD_template_module,                Only : FD_template

    ! HACK
    Use FD_Laplacian_3d_module, Only : FD_Laplacian_3d
    
    Implicit None

    Real( wp ), Parameter :: stencil_tol = 1e-15_wp

    Class( equation_solver_hypre_pfmg ), Intent( InOut ) :: method
    Class( comms_base_class       )    , Intent( In    ) :: comms
    Integer, Dimension( 1:3 )          , Intent( In    ) :: n
    Integer, Dimension( 1:3 )          , Intent( In    ) :: lb
    Integer, Dimension( 1:3 )          , Intent( In    ) :: ub
    Class( FD_template ),                Intent( In    ) :: FD_operator

    Real( wp ), Dimension( :, :, : ), Allocatable :: stencil

    Real( wp ), Dimension( : ), Allocatable :: stencil_values
    
    Integer, Dimension( :, : ), Allocatable :: stencil_elements

    Integer :: comm
    Integer :: n_stencil
    Integer :: i1, i2, i3

    Call comms%get_comm( comm )
    
    ! HACK - should fix FD library
    Select Type( FD_operator )
    Class is ( FD_Laplacian_3d )
       Call FD_operator%get_stencil( stencil )
    Class Default
       Error Stop "Unknown class"
    End Select

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
    Write( *, * ) stencil_values( 1:10 )

    method%hypre_objects = ssp_hypre_struct_setup( comm, n, lb, ub, n_stencil, stencil_elements, stencil_values )

  End Subroutine pfmg_init

  Subroutine pfmg( method, &
       lb, ub, b, rtol, &
       x, istop, istop_message, itn, rnorm )

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

    Integer, Dimension( 1:3 ) :: n

    n = Ubound( b ) - Lbound( b ) + 1

    x = 0.1_wp
    Call ssp_hypre_struct_pfmg_solve( method%hypre_objects, n( 1 ), n( 2 ), n( 3 ), b, x, itn, rnorm, istop )

    istop_message = "Who knows?"

  End Subroutine pfmg

  Subroutine pfmg_destroy( method )

    Implicit None

    Type( equation_solver_hypre_pfmg ), Intent( InOut ) :: method

    Call  ssp_hypre_struct_free( method%hypre_objects )

  End Subroutine pfmg_destroy
  
End Module equation_solver_hypre_pfmg_module
