Module quadrature_trapezium_rule_module

  Use constants,              Only : wp
  Use quadrature_base_module, Only : quadrature_base_class

  Implicit None

  Type, Public, Extends( quadrature_base_class ) :: quadrature_trapezium_rule
     Private
   Contains
     Procedure, Public :: integrate => trapezium_rule
  End type quadrature_trapezium_rule

Contains

  Function trapezium_rule( q, comms, l, n_grid, grid ) Result( r )

    Use comms_base_class_module, Only : comms_base_class
    Use lattice_module,          Only : lattice

    Real( wp ) :: r

    Class( quadrature_trapezium_rule   ), Intent( In ) :: q
    Class( comms_base_class            ), Intent( In ) :: comms
    Type( lattice ),                      Intent( In ) :: l
    Integer,    Dimension( :       ),     Intent( In ) :: n_grid
    Real( wp ), Dimension( :, :, : ),     Intent( In ) :: grid

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

    Call comms%reduce( r )

  End Function trapezium_rule

End Module quadrature_trapezium_rule_module
