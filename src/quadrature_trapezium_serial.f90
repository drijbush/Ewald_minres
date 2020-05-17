Module quadrature_trapezium_serial_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Use quadrature_base_module, Only : quadrature_base_class

  Implicit None

  Type, Public, Extends( quadrature_base_class ) :: quadrature_trapezium_serial
     Private
   Contains
     Procedure, NoPass, Public :: integrate => trapezium_serial
  End type quadrature_trapezium_serial

Contains

  Function trapezium_serial( l, n_grid, grid ) Result( r )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use lattice_module, Only : lattice

    Real( wp ) :: r

    Type( lattice )                 , Intent( In ) :: l
    Integer   , Dimension( :       ), Intent( In ) :: n_grid
    Real( wp ), Dimension( :, :, : ), Intent( In ) :: grid

    Real( wp ) :: dV
    Real( wp ) :: c, y, t

    Integer :: i1, i2, i3

    dV =  l%get_volume() / Product( n_grid )

    ! Integrate using Kahan summation for accuracy
    r = 0.0_wp
    !$omp parallel default( none ) shared( grid, r ) private( c, y, t, i1, i2, i3 )
    c = 0.0_wp
    !$omp do collapse( 3 ) reduction( +:r )
    Do i3 = Lbound( grid, Dim = 3 ), Ubound( grid, Dim = 3 )
       Do i2 = Lbound( grid, Dim = 2 ), Ubound( grid, Dim = 2 )
          Do i1 = Lbound( grid, Dim = 1 ), Ubound( grid, Dim = 1 )
             y = grid( i1, i2, i3 ) - c
             t = r + y
             c = ( t - r ) - y
             r = t
          End Do
       End Do
    End Do
    !$omp end do
    !$omp end parallel
    r = r * dV    
    
  End Function trapezium_serial
  
End Module quadrature_trapezium_serial_module
