Module charge_grid_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Use lattice_module, Only : lattice

  Implicit None

  Public :: charge_grid_calculate, charge_grid_find_range

  Private

  Real( wp ), Parameter, Private :: pi = 3.141592653589793238462643383279502884197_wp

Contains

    Subroutine charge_grid_calculate( l, alpha, q, r, range_gauss, q_grid )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use lattice_module, Only : lattice

    Implicit none
    
    Type( lattice )                    , Intent( In    ) :: l
    Real( wp ),                          Intent( In    ) :: alpha
    Real( wp ), Dimension( 1: ),         Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1: ),     Intent( In    ) :: r
    Integer   , Dimension( 1:3        ), Intent( In    ) :: range_gauss
    Real( wp ), Dimension( 0:, 0:, 0: ), Intent(   Out ) :: q_grid

    Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp

    Real( wp ), Dimension( 1:3 ) :: ri
    Real( wp ), Dimension( 1:3 ) :: fi
    Real( wp ), Dimension( 1:3 ) :: f_point
    Real( wp ), Dimension( 1:3 ) :: r_point
    Real( wp ), Dimension( 1:3 ) :: grid_vec

    Real( wp ) :: q_norm
    Real( wp ) :: qi_norm
    Real( wp ) :: q_val
    
    Integer, Dimension( 1:3 ) :: n_grid
    Integer, Dimension( 1:3 ) :: i_atom_centre
    Integer, Dimension( 1:3 ) :: i_atom_grid
    Integer, Dimension( 1:3 ) :: i_point
    Integer, Dimension( 1:3 ) :: i_grid

    Integer :: n
    Integer :: i1, i2, i3
    Integer :: i

    n = Size( q )
    n_grid = Ubound( q_grid ) + 1

    q_norm = ( ( ( alpha * alpha ) / pi ) ** 1.5_wp )
    
    !$omp parallel default( none ) shared( n, l, alpha, r, q, q_norm, n_grid, q_grid, range_gauss ) &
    !$omp                          private( i, i1, i2, i3, qi_norm, ri, fi, i_atom_centre,   &
    !$omp                                   i_atom_grid, i_point, f_point, r_point, grid_vec, q_val, i_grid )
    !$omp do collapse( 3 )
    Do i3 = 0, n_grid( 3 ) - 1
       Do i2 = 0, n_grid( 2 ) - 1
          Do i1 = 0, n_grid( 1 ) - 1
             q_grid( i1, i2, i3 ) = 0.0_wp
          End Do
       End Do
    End Do
    !$omp end do
    ! Loop over atoms
    !$omp do reduction( +:q_grid )
    Do i = 1, n
       ! Loop over points associated with atoms
       ! Find point nearest to the atom, and call this the centre for the atom grid
       ! Assumes atom in fractional 0 < ri < 1
       ri = r( :, i )
       qi_norm = q_norm * q( i )
       Call l%to_fractional( ri, fi )
       i_atom_centre = Nint( fi * n_grid )
       Do i3 = - range_gauss( 3 ), range_gauss( 3 )
          Do i2 = - range_gauss( 2 ), range_gauss( 2 )
             Do i1 = - range_gauss( 1 ), range_gauss( 1 )
                i_atom_grid = [ i1, i2, i3 ]
                ! The indices of the point in space
                i_point = i_atom_centre + i_atom_grid
                ! Transform to fractional coordinates
                f_point = Real( i_point, wp ) / n_grid
                ! And fractional to real
                Call l%to_direct( f_point, r_point )
                ! Calculate the contribution to the total charge at the point
                ! R_POINT due to the charge distribution I
                grid_vec = r_point - ri
                q_val = qi_norm * Exp( - alpha * alpha * Dot_product( grid_vec, grid_vec ) )
                ! Reflect grid indices into reference Cell
                i_grid = Modulo( i_point, n_grid )
                ! And add in
                q_grid( i_grid( 1 ), i_grid( 2 ), i_grid( 3 ) ) = &
                     q_grid( i_grid( 1 ), i_grid( 2 ), i_grid( 3 ) ) + q_val
             End Do
          End Do
       End Do
    End Do
    !$omp end do
    !$omp end parallel
    
  End Subroutine charge_grid_calculate

  Subroutine charge_grid_find_range( l, alpha, n_grid, range_gauss )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

    Use lattice_module, Only : lattice

    Implicit None

    Type( lattice )             , Intent( In    ) :: l
    Real( wp ),                   Intent( In    ) :: alpha
    Integer   , Dimension( 1:3 ), Intent( In    ) :: n_grid
    Integer   , Dimension( 1:3 ), Intent(   Out ) :: range_gauss

    Real( wp ), Dimension( 1:3 ) :: G
    Real( wp ), Dimension( 1:3 ) :: grid_vec

    Real( wp ) :: q_norm
    Real( wp ) :: q_val

    Integer, Dimension( 1:3 ) :: n_vec

    Integer :: i_dir

    n_vec = [ 1, 0, 0 ] ! Will use cshift to move onto next vector
    q_norm = ( ( ( alpha * alpha ) / pi ) ** 1.5_wp )
    Do i_dir = 1, 3
       Call l%get_dir_vec( n_vec, G )
       range_gauss( i_dir ) = 1
       Do
          grid_vec = ( range_gauss( i_dir ) * G ) / ( n_grid( i_dir ) )
          q_val = q_norm * Exp( - alpha * alpha * Dot_product( grid_vec, grid_vec ) )
          If( q_val < 1e-15_wp ) Exit
          range_gauss( i_dir ) = range_gauss( i_dir ) + 1
       End Do
       n_vec = Cshift( n_vec, -1 )
    End Do

  End Subroutine charge_grid_find_range

End Module charge_grid_module
