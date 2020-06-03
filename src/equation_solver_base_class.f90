Module equation_solver_base_class_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Use comms_base_class_module, Only : comms_base_class
  Use halo_setter_base_module, Only : halo_setter_base_class
  Use FD_template_module     , Only : FD_template
  
  Implicit None
       
  Type, Public, Abstract :: equation_solver_base_class
     Integer                                         , Public :: max_iter = 1000         
     Class( comms_base_class           ), Allocatable, Public :: comms
     Class( FD_template                ), Allocatable, Public :: FD_operator
     Class( halo_setter_base_class     ), Allocatable, Public :: halo_swapper
     ! Bug in gcc 7.4 gives ICE with allocatable
     ! Therefore split off in own module
!!$     Class( equation_solver_base_class ), Pointer    , Private :: preconditioner
   Contains
     Generic                      , Public  :: init => base_init
     Procedure( solver ), Deferred, Public  :: solve
     Procedure          , NoPass  , Public  :: contract
     Procedure                    , Private :: base_init
  End type equation_solver_base_class

  Abstract Interface

     Subroutine solver( method, &
          lb, ub, b, rtol, &
          x, istop, istop_message, itn, rnorm )
       Use halo_setter_base_module, Only : halo_setter_base_class
       Use FD_template_module     , Only : FD_template
       Import :: wp
       Import :: equation_solver_base_class
       Class( equation_solver_base_class )                   , Intent( InOut ) :: method
       Integer,  Dimension( 1:3 )                            , Intent( In    ) :: lb( 1:3 )
       Integer,  Dimension( 1:3 )                            , Intent( In    ) :: ub( 1:3 )
       Real( wp ) , Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ) :: b
       Real( wp )                                            , Intent( In    ) :: rtol
       Real( wp ) , Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: x
       Integer                                               , Intent(   Out ) :: istop
       Character( Len = * )                                  , Intent(   Out ) :: istop_message
       Real( wp )                                            , Intent(   Out ) :: rnorm
       Integer                                               , Intent(   Out ) :: itn
    End Subroutine solver
     
  End Interface

  Private
  
Contains

  Subroutine base_init( method, max_iter, comms, FD_operator, halo_swapper )

    Use comms_base_class_module          , Only : comms_base_class
    Use halo_setter_base_module          , Only : halo_setter_base_class
    Use FD_template_module               , Only : FD_template
    
    Implicit None

    Class( equation_solver_base_class ), Intent( InOut )           :: method
    Integer                            , Intent( In    ), Optional :: max_iter
    Class( comms_base_class       )    , Intent( In    ), Optional :: comms
    Class( FD_template )               , Intent( In    ), Optional :: FD_operator
    Class( halo_setter_base_class )    , Intent( In    ), Optional :: halo_swapper

    If( Present( max_iter ) ) Then
       method%max_iter = max_iter
    End If

    If( Present( comms ) ) Then
       method%comms = comms
    End If

    If( Present( FD_operator ) ) Then
       method%FD_operator = FD_operator
    End If
    
    If( Present( halo_swapper ) ) Then
       method%halo_swapper = halo_swapper
    End If
    
  End Subroutine base_init
  
  Function contract( comms, x, y ) Result( d )

    Use comms_base_class_module, Only : comms_base_class

    Implicit None

    Real( wp ) :: d

    Class( comms_base_class )       , Intent( In ) :: comms
    Real( wp ), Dimension( :, :, : ), Intent( In ) :: x
    Real( wp ), Dimension( :, :, : ), Intent( In ) :: y

    Real( wp ) :: c, yk, t

    Integer :: i1, i2, i3

    ! Sum using Kahan summation for accuracy
    d = 0.0_wp
    !$omp parallel default( none ) shared( grid, d ) private( c, yk, t, i1, i2, i3 )
    c = 0.0_wp
    !$omp do collapse( 3 ) reduction( +:d )
    Do i3 = Lbound( x, Dim = 3 ), Ubound( x, Dim = 3 )
       Do i2 = Lbound( x, Dim = 2 ), Ubound( x, Dim = 2 )
          Do i1 = Lbound( x, Dim = 1 ), Ubound( x, Dim = 1 )
             yk = x( i1, i2, i3 ) * y( i1, i2, i3 ) - c
             t  = d + yk
             c  = ( t - d ) - yk
             d  = t
          End Do
       End Do
    End Do
    !$omp end do
    !$omp end parallel    ! Use Kahan for accuracy

    Call comms%reduce( d )

  End Function contract
   
End Module equation_solver_base_class_module
