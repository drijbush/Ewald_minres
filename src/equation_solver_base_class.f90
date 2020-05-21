Module equation_solver_base_class_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Type, Public, Abstract :: equation_solver_base_class
   Contains
     Procedure( solver ), Deferred, Public :: solve
     Procedure          , NoPass  , Public :: contract
  End type equation_solver_base_class

  Abstract Interface

     Subroutine solver( method, &
          lb, ub, FD_operator, comms, halo_swapper, Msolve, b, itnlim, rtol,  precon, &
          x, istop, istop_message, itn, rnorm )
       Use comms_base_class_module, Only : comms_base_class
       Use halo_setter_base_module, Only : halo_setter_base_class
       Use FD_template_module     , Only : FD_template
       Import :: wp
       Import :: equation_solver_base_class
       Class( equation_solver_base_class )                   , Intent( In    ) :: method
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
    End Subroutine solver
     
  End Interface

Contains
  
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