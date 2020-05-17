Module fast_fourier_poisson_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

  Implicit None

  Public :: ffp_long_range
  Public :: ffp_sic
  
  Private
  
  Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp

Contains

  Subroutine ffp_long_range( l, q, r, alpha, FD_order, q_halo, r_halo, recip_E, q_grid, pot_grid, ei, f, t_grid, t_recip, error )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

    Use lattice_module        , Only : lattice
    Use charge_grid_module    , Only : charge_grid_calculate, charge_grid_find_range, charge_grid_forces
    Use FFT_module            , Only : fft_fft3d
    Use halo_serial_module     , Only : halo_serial_setter
    Use halo_setter_base_module, Only : halo_setter_base_class

    Implicit None

    Type( lattice )                    , Intent( In    ) :: l
    Real( wp ), Dimension( 1:     )    , Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1:  )   , Intent( In    ) :: r
    Real( wp )                         , Intent( In    ) :: alpha
    Integer                            , Intent( In    ) :: FD_order
    Real( wp ), Dimension( 1:     )    , Intent( In    ) :: q_halo
    Real( wp ), Dimension( 1:, 1:  )   , Intent( In    ) :: r_halo
    Real( wp )                         , Intent(   Out ) :: recip_E
    Real( wp ), Dimension( 0:, 0:, 0: ), Intent(   Out ) :: q_grid
    Real( wp ), Dimension( 0:, 0:, 0: ), Intent(   Out ) :: pot_grid
    Real( wp ), Dimension( 1: )        , Intent(   Out ) :: ei
    Real( wp ), Dimension( 1:, 1: )    , Intent(   Out ) :: f
    Real( wp )                         , Intent(   Out ) :: t_grid
    Real( wp )                         , Intent(   Out ) :: t_recip
    Integer                            , Intent(   Out ) :: error

    Class( halo_setter_base_class      ), Allocatable :: pot_swapper
    
    Complex( wp ), Dimension( :, :, : ), Allocatable :: sfac
    Complex( wp ), Dimension( :, :, : ), Allocatable :: pot_k

    Real( wp ), Dimension( :, : ), Allocatable :: r_full
    
    Real( wp ), Dimension( 1:3 ) :: G

    Real( wp ) :: G_len_sq
    Real( wp ) :: gauss_term
    
    Integer, Dimension( 1:3 ) :: n_grid
    Integer, Dimension( 1:3 ) :: range_gauss
    Integer, Dimension( 1:3 ) :: iG
    
    Integer :: iG1, iG2, iG3

    Integer( li ) :: start, finish, rate
    
    Allocate( r_full( 1:3, 1:Size( q ) + Size( q_halo ) ) )
    r_full( :, 1:Size( q )    ) = r
    r_full( :, Size( q ) + 1: ) = r_halo
    
    error = 0

    n_grid = Ubound( q_grid ) + 1

    ! Grid the charge
    Call system_clock( start, rate )
    ! First find range of the gaussian along each of the axes of the grid
    Call charge_grid_find_range( l, alpha, n_grid, range_gauss )
    ! Now grid the charge
    Call charge_grid_calculate( l, alpha, [ q, q_halo ], r_full, range_gauss, &
         Lbound( q_grid ), Ubound( q_grid ), q_grid, error )
    Call system_clock( finish, rate )
    t_grid = Real( finish - start, wp ) / rate

    ! Now calculate the long range term by Fourier transform methods

    ! First calculate the structure factor
    Call system_clock( start, rate )
    Allocate( sfac( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
    sfac = q_grid
    Call fft_fft3d( +1, sfac )
    
    ! Now calculate the potential in fourier space
    Allocate( pot_k( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
    Do iG3 = 0, n_grid( 3 ) - 1
       Do iG2 = 0, n_grid( 2 ) - 1
          Do iG1 = 0, n_grid( 1 ) - 1
             iG = [ iG1, iG2, iG3 ]
             If( All( iG == 0 ) ) Then
                pot_k( 0, 0, 0 ) = 0.0_wp
             Else
                Where( iG > n_grid / 2 )
                   iG = iG - n_grid
                End Where
                Call l%get_rec_vec( iG, G )
                G_len_sq = Dot_product( G, G )
                gauss_term = 1.0_wp / G_len_sq
                pot_k( iG1, iG2, iG3 ) = sfac( iG1, iG2, iG3 ) * gauss_term
             End If
          End Do
       End Do
    End Do
    Call fft_fft3d( -1, pot_k )
    pot_grid = Real( pot_k, wp ) / ( pi * Product( n_grid ) )
    Call system_clock( finish, rate )
    t_recip = Real( finish - start, wp ) / rate
    
    ! Calculate from the potential the long range energy
    ! Two minus signs as a) this is the SCREENED charge
    !                    b) I like to just add up the energies, having to subtract it is confusing
    recip_E = - 0.5_wp * Sum( - q_grid * pot_grid ) * ( l%get_volume() / Product( n_grid ) )
    
    ! Calculate the forces and energy per site
    ! For now while implementing halos
   Allocate( halo_serial_setter :: pot_swapper )
   Call pot_swapper%init( error )
   Call charge_grid_forces( l, alpha, q, r, range_gauss, pot_swapper, Lbound( pot_grid ), Ubound( pot_grid ), &
         pot_grid, ei, f )
    
  End Subroutine ffp_long_range

  Pure Function ffp_sic( q, alpha ) Result( sic )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Real( wp ), Dimension( : ), Intent( In    ) :: q
    Real( wp )                , Intent( In    ) :: alpha

    Real( wp ) :: sic

    sic = - alpha * Sum( q * q ) / Sqrt( 2.0_wp * pi )

  End Function ffp_sic

End Module fast_fourier_poisson_module
