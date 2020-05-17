Module quadrature_base_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Type, Public, Abstract :: quadrature_base_class
     Private
   Contains
     Procedure( integrate ), Deferred, Public :: integrate
  End Type quadrature_base_class

  Abstract Interface
     
     Function integrate( q, comms, l, n_grid, grid ) Result( r )
       Use comms_base_class_module, Only : comms_base_class
       Use lattice_module         , Only : lattice
       Import :: wp
       Import :: quadrature_base_class
       Real( wp ) :: r
       Class( quadrature_base_class   ), Intent( In ) :: q
       Class( comms_base_class        ), Intent( In ) :: comms
       Type( lattice )                 , Intent( In ) :: l
       Integer   , Dimension( :       ), Intent( In ) :: n_grid
       Real( wp ), Dimension( :, :, : ), Intent( In ) :: grid
     End Function integrate
     
  End Interface

  Private
  
End Module quadrature_base_module
