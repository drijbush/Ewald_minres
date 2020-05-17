Module quadrature_base_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Type, Public, Abstract :: quadrature_base_class
     Private
   Contains
     Procedure( integrate ), Deferred, NoPass, Public :: integrate
  End Type quadrature_base_class

  Abstract Interface
     
     Function integrate( l, n_grid, grid ) Result( r )
       Use lattice_module, Only : lattice
       Import :: wp
       Real( wp ) :: r
       Type( lattice )                 , Intent( In ) :: l
       Integer   , Dimension( :       ), Intent( In ) :: n_grid
       Real( wp ), Dimension( :, :, : ), Intent( In ) :: grid
     End Function integrate
     
  End Interface

  Private
  
End Module quadrature_base_module
