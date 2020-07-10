Module grid_io_module

  Implicit None

  Public :: grid_io_save

  Private

Contains

  Subroutine grid_io_save( unit, filename, l, grid )

    Use constants, Only : wp
    Use lattice_module, Only : lattice

    Integer,                             Intent( In ) :: unit
    Character( Len = * ),                Intent( In ) :: filename
    Type( lattice ),                     Intent( In ) :: l
    Real( wp ), Dimension( 0:, 0:, 0: ), Intent( In ) :: grid

    Real( wp ), Dimension( :,  : ), Allocatable  :: lattice_vectors

    Real( wp ), Dimension( 1:3 ) :: f_point
    Real( wp ), Dimension( 1:3 ) :: r_point

    Integer, Dimension( 1:3 ) :: n_grid
    Integer, Dimension( 1:3 ) :: i_grid

    Integer :: i
    Integer :: i1, i2, i3

    n_grid = Ubound( grid ) + 1

    lattice_vectors = l%get_direct_vectors()

    Open( unit, file = filename )
    Write( unit, * ) 'Sum over grid           ', Sum( grid )
    Write( unit, * ) 'Average per unit volume ', Sum( grid ) / l%get_volume()
    Write( unit, '( 3( i0, 1x ), t60, a )' ) n_grid, '# grid size'
    Do i = 1, Ubound( lattice_vectors,  Dim = 2 )
       Write( unit, '( 3( f11.6, 1x ), t60, a, i0  )' ) lattice_vectors( :, i ), '# lattice vector ', i
    End Do
    Do i3 = 0, n_grid( 3 ) - 1
       Do i2 = 0, n_grid( 2 ) - 1
          Do i1 = 0, n_grid( 1 ) - 1
             i_grid = [ i1, i2, i3 ]
             f_point = Real( i_grid, wp ) / n_grid
             Call l%to_direct( f_point, r_point )
             Write( unit, '( 3( f6.2, 1x ), 5x, f30.16 )' ) r_point, grid( i1, i2, i3 )
          End Do
       End Do
    End Do
    Close( unit )

  End Subroutine grid_io_save


End Module grid_io_module
