Module quadrature_trapezium_parallel_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Use quadrature_trapezium_serial_module, Only : quadrature_trapezium_serial

  Implicit None

  Type, Public, Extends( quadrature_trapezium_serial ) :: quadrature_trapezium_parallel
     Private
   Contains
     Procedure, NoPass, Public :: integrate => trapezium_parallel
  End type quadrature_trapezium_parallel

Contains

  Function trapezium_parallel( l, n_grid, grid ) Result( r )
    
    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use lattice_module, Only : lattice

    Real( wp ) :: r

    Type( lattice )                 , Intent( In ) :: l
    Integer   , Dimension( :       ), Intent( In ) :: n_grid
    Real( wp ), Dimension( :, :, : ), Intent( In ) :: grid

    Type( quadrature_trapezium_serial ) :: g
    
    Real( wp ) :: r_loc

    r_loc = g%integrate( l, n_grid, grid )

!!$    Call global_sum( r_loc )
    r = r_loc
    
  End Function trapezium_parallel

End Module quadrature_trapezium_parallel_module
