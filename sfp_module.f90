Module sfp_module

  Implicit None

  Public :: sfp_long_range
  
  Private

Contains

  Subroutine sfp_long_range(  l, q, r, alpha, ew_func, q_grid, pot_grid, recip_E, t_grid, t_recip )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64
    
    Use lattice_module        , Only : lattice
    Use charge_grid_module    , Only : charge_grid_calculate, charge_grid_find_range

    Implicit None

    Type( lattice )                    , Intent( In    ) :: l
    Real( wp )   , Dimension( :     )     , Intent( In    ) :: q
    Real( wp )   , Dimension( :, :  )     , Intent( In    ) :: r
    Real( wp )                            , Intent( In    ) :: alpha
    Complex( wp ), Dimension( 0: )        , Intent( In    ) :: ew_func
    Real( wp )   , Dimension( 0:, 0:, 0: ), Intent(   Out ) :: q_grid
    Real( wp )   , Dimension( 0:, 0:, 0: ), Intent(   Out ) :: pot_grid
    Real( wp )                            , Intent( Out   ) :: recip_E
    Real( wp )                            , Intent( Out   ) :: t_grid
    Real( wp )                            , Intent( Out   ) :: t_recip

    Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp

    Complex( wp ) :: pot
    
    Real( wp ), Dimension( 1:3 ) :: f_point
    Real( wp ), Dimension( 1:3 ) :: r_point
    Real( wp ), Dimension( 1:3 ) :: G
    Real( wp ), Dimension( 1:3 ) :: ri

    Real( wp ) :: potr

    Integer, Dimension( 1:3 ) :: range_gauss
    Integer, Dimension( 1:3 ) :: i_grid
    
    Integer, Dimension( 1:3 ) :: n_grid
    
    Integer :: iG
    Integer :: i1, i2, i3
    
    Integer( li ) :: start, finish, rate

    n_grid = Ubound( q_grid ) + 1
    
    ! Grid the charge
    Call system_clock( start, rate )
    ! First find range of the gaussian along each of the axes of the grid
    Call charge_grid_find_range( l, alpha, n_grid, range_gauss )
    ! Now grid the charge
    Call charge_grid_calculate( l, alpha, q, r, range_gauss, q_grid )
    Call system_clock( finish, rate )
    t_grid = Real( finish - start, wp ) / rate
    
    ! Calculate the potential at each point - CAN BE SLOW!!!!
    Call system_clock( start, rate )
    !$omp parallel default( none ) shared( l, n_grid, ew_func, pot_grid ) &
    !$omp                          private( i1, i2, i3, i_grid, f_point, r_point, &
    !$omp                                   ri, potr, ig, G, pot )
    !$omp do collapse( 3 )
    Do i3 = 0, n_grid( 3 ) - 1
       Do i2 = 0, n_grid( 2 ) - 1
          Do i1 = 0, n_grid( 1 ) - 1
             i_grid = [ i1, i2, i3 ]
             f_point = Real( i_grid, wp ) / n_grid
             Call l%to_direct( f_point, r_point )
             ri = r_point
             potr = 0.0_wp
             Do iG = 1, Ubound( ew_func, Dim = 1 )
                Call l%get_nth_rec_vec( ig + 1, G )
                G = G * 2.0_wp * pi
                pot = Exp( Cmplx( 0.0_wp, Dot_product( G, ri ), wp ) ) * ew_func( ig )
                potr = potr + Real( pot + Conjg( pot ), wp )
             End Do
             potr = potr * 4.0_wp * pi / l%get_volume()
             pot_grid( i1, i2, i3 ) = potr
          End Do
       End Do
    End Do
    !$omp end do
    !$omp end parallel
    !HACK
    ! Where * V / 4 * pi come from ?
    pot_grid = pot_grid * l%get_volume() / ( 4.0_wp * pi )
    Call system_clock( finish, rate )
    t_recip = Real( finish - start, wp ) / rate

    ! And hence the energy
    ! Minus Sign on charge grid as we integrate against the screening charge, which
    ! has opposite sign to the actual charge (as it is screening, Duh!
    ! Minus sign on whole comes here just so can add up all contribs at end rather than confuse myself about signs
    recip_E = - 0.5_wp * Sum( - q_grid * pot_grid ) * ( l%get_volume() / Product( n_grid ) )

  End Subroutine sfp_long_range
  
End Module sfp_module
