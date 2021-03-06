Program test

  !$ Use omp_lib
  
  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

  Use lattice_module                      , Only : lattice
  Use charge_grid_module                  , Only : charge_grid_find_range
  Use fast_fourier_poisson_module         , Only : ffp_long_range, ffp_sic
  Use symetrically_screened_poisson_module, Only : ssp_long_range, ssp_sic
  Use sfp_module                          , Only : sfp_long_range, sfp_sic
  Use real_space_module                   , Only : real_space_energy
  Use grid_io_module                      , Only : grid_io_save
  
  Implicit None

  Logical, Parameter :: do_sfp = .False.
  
  Type( lattice ) :: l

  Complex( wp ), Dimension( : ), Allocatable :: ew_func

  Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp

  Real( wp ), Parameter :: r4pie0 = 138935.48350000000_wp

  Real( wp ) :: gauss_tol = 1e-15_wp
  
  Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_ssp
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_ffp
  
  Real( wp ), Dimension( :, : ), Allocatable :: a
  Real( wp ), Dimension( :, : ), Allocatable :: r
  Real( wp ), Dimension( :, : ), Allocatable :: force_ffp
  Real( wp ), Dimension( :, : ), Allocatable :: force_ssp

  Real( wp ), Dimension( : ), Allocatable :: q
  Real( wp ), Dimension( : ), Allocatable :: ei_ffp
  Real( wp ), Dimension( : ), Allocatable :: ei_ssp

  Real( wp ), Dimension( 1:3 ) :: t
  Real( wp ), Dimension( 1:3 ) :: dG

  Real( wp ) :: alpha
  Real( wp ) :: sic
  Real( wp ) :: recip_E, real_E, tot_E
  Real( wp ) :: recip_E_ffp, sic_ffp, real_E_ffp, tot_E_ffp
  Real( wp ) :: recip_E_ssp, sic_ssp, real_E_ssp, tot_E_ssp
  Real( wp ) :: recip_E_sfp, sic_sfp, real_E_sfp, tot_E_sfp
  Real( wp ) :: xshift
  Real( wp ) :: rms_delta_pot, rms_delta_force, q_error
  Real( wp ) :: t_grid, t_real, t_recip
  
  Integer, Dimension( 1:3 ) :: n_grid

  Integer :: range_gauss
  Integer :: FD_order
  Integer :: n, level
  Integer :: nd
  Integer :: i, j
  Integer :: max_G_shells = 2
  Integer :: error
  
  Integer( li ) :: start, finish, rate

  Logical :: fexist
  Logical :: l_bad_charge
  
  Write( *, * ) 'Ewald param ?'
  Read ( *, * ) alpha
!!$  Write( *, * ) 'Grid Dims ?'
!!$  Read ( *, * ) n_grid
  Write( *, * ) 'Number of ppoints for gaussians?'
  Read ( *, * ) range_gauss
  Write( *, * ) 'FD_order?'
  Read ( *, * ) FD_order
  Write( *, * ) 'xshift?'
  Read ( *, * ) xshift

  Write( *, * ) 'Param, dims, order ', alpha, range_gauss, FD_order, xshift
  
  !$ Write( *, * ) 'Running on ', omp_get_max_threads(), ' threads'
  
  nd = 3

  Allocate( a( 1:3, 1:nd ) )

  Open( 10, File = 'CONFIG' )
  Call read_header( 10, n, level, a )
  Write( *, * ) 'n = ', n
  Allocate( q( 1:n ) )
  Allocate( r( 1:3, 1:n ) )
  Call read_config( level, 10, q, r )
  
  Call l%initialise( nd, a, alpha )
  
  Write( *, '( a, 3( t25, 3( f10.5, 1x ) / ) )' ) &
       'Lattice vectors: ', l%get_direct_vectors()
  
  Call get_n_grid( l, alpha, range_gauss, gauss_tol, n_grid )
  Write( *, * ) 'N_grid = ', n_grid
  
!!$  Call l%print
  Do i = 1, 3
     dG( i ) = Sqrt( Dot_product( a( :, i ), a( :, i ) ) ) / n_grid( i )
  End Do

  Write( *, * ) 'r(1) before shift and reference = ', r( :, 1 )
  Do i = 1, n
     t = r( :, i )
     t( 1 ) = t( 1 ) + xshift
     Call l%to_reference( t, r( :, i ) ) 
  End Do
  Write( *, * ) 'r(1) after  shift and reference = ', r( :, 1 )

  ! How big is my gaussian?
!!$  Call  charge_grid_find_range( l, alpha, n_grid, range_gauss )
!!$  Write( *, * ) 'Gauss range: ', range_gauss
  
  ! Generate the useful fourier space function
  Call system_clock( start, rate )
  Allocate( ew_func( 0:l%get_n_rec_vecs() - 1 ) )
  Call generate_ew_func( l, q, r, alpha, ew_func )
  Call system_clock( finish, rate )
  Write( *, * ) 'ewfunc time ', Real( finish - start, wp ) / rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  TRADITIONAL EWALD SECTION
  Call trad_ewald( l, q, r, alpha, ew_func, recip_E, sic, real_E, tot_E, t_recip, t_real )
  Write( *, * ) 'Trad Ewald long  range time: ', t_recip
  Write( *, * ) 'Trad Ewald short range time: ', t_real

  ! END TRADITIONAL EWALD SECTION
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!
  !
  ! FAST FOURIER POISSON (FFP) - grid charge, FFT, get pot in k space, FFT back, energy
  !

  Allocate( q_grid( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  ! Calculate the long range term by fourier
  Allocate( pot_grid_ffp( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  Allocate( ei_ffp( 1:n ) )
  Allocate( force_ffp( 1:3, 1:n ) )
  Call ffp_long_range( l, q, r, alpha, FD_order, recip_E_ffp, q_grid, pot_grid_ffp, ei_ffp, force_ffp, t_grid, t_recip, error )
  Write( *, * ) 'Nett force ', Sum( force_ffp( 1, : ) ), Sum( force_ffp( 2, : ) ), Sum( force_ffp( 3, : ) )
  Open( 11, file = 'forces_ffp.dat' )
  Write( 11, * ) n, '     #number of particles'
  Write( 11, * ) Sum( force_ffp( 1, : ) ), Sum( force_ffp( 2, : ) ), Sum( force_ffp( 3, : ) ), '     #Nett force'
  Do i = 1, n
     Write( 11, * ) i, force_ffp( :, i )
  End Do
  Close( 11 )
  Write( *, * ) 'FFP grid  time: ', t_grid
  Write( *, * ) 'FFP solve time: ', t_recip
  l_bad_charge = .False.
  Select Case( error )
  Case( -1 )
     Write( *, * ) 'WARNING: Sum of charge on grid greater than tolerance but "fixed" (for some value of fixed)'
     l_bad_charge = .True.
  Case( -2 )
     Write( *, * ) 'WARNING: Sum of charge on grid greater than tolerance AND could not be fixed'
     l_bad_charge = .True.
     Stop
  End Select

  ! Check how well the charge grid adds up to zero
  q_error = Sum( q_grid )

  ! Save the FFP potential
  Call grid_io_save( 11, 'pot_grid_FFP.dat', l, pot_grid_ffp )

  Write( *, * ) 'FFP: Sum of charge over grid: ', Sum( q_grid )
  Write( *, * ) 'FFP: Sum over pot grid      : ', Sum( pot_grid_ffp )

  ! SIC 
  sic_ffp = ffp_sic( q, alpha )

  ! Calculate short range
  Call system_clock( start, rate )
  !$omp parallel default( none ) shared( l, q, r, alpha, max_G_shells, real_E_ffp )
  Call real_space_energy( l, q, r, alpha / Sqrt( 2.0_wp ), max_G_shells, real_E_ffp )
  !$omp end parallel
  Call system_clock( finish, rate )
  Write( *, * ) 'FFP real time: ', Real( finish - start, wp ) / rate

  ! And hence total energy 
  tot_E_ffp = real_E_ffp + sic_ffp + recip_E_ffp

  ! END FFP SECTION
  !
!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!
  !
  ! SLOW FOURIER POISSON SOLVER (SFP) i.e. no FFT
  ! Can be slow so optional

  If( do_sfp ) Then
     ! Fourier space energy
     Allocate( pot_grid( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
     Call sfp_long_range(  l, q, r, alpha, ew_func, recip_E_sfp, q_grid, pot_grid, t_grid, t_recip, error )
     l_bad_charge = .False.
     Select Case( error )
     Case( -1 )
        Write( *, * ) 'WARNING: Sum of charge on grid greater than tolerance but "fixed" (for some value of fixed)'
        l_bad_charge = .True.
     Case( -2 )
        Write( *, * ) 'WARNING: Sum of charge on grid greater than tolerance AND could not be fixed'
        l_bad_charge = .True.
        Stop
     End Select

     Write( *, * ) 'SFP grid  time: ', t_grid
     Write( *, * ) 'SFP solve time: ', t_recip
     ! Save the q grid to file
     Call grid_io_save( 11, 'q_grid.dat', l, q_grid )
     
     Write( *, * ) 'SFP: Sum of charge over grid: ', Sum( q_grid )
     Write( *, * ) 'SFP: Sum over pot grid      : ', Sum( pot_grid )

     ! Save the pot grid to file
     Call grid_io_save( 11, 'pot_grid.dat', l, pot_grid )

     ! SFP self interaction correction
     sic_sfp = sfp_sic( q, alpha )
     
     ! SFP Real Space - same as FFP
     real_E_sfp = real_E_ffp
  
     ! And hence the total energy
     tot_E_sfp = real_E_sfp + sic_sfp + recip_E_sfp
  End If

  ! END SFP SECTION
  !
!!!!!!!!!!!!!

!!!!!!!!!!
  !
  ! SYMETRICALLY SCREENED POISSON (SSP) - purely real space methods used, finite difference to solve eqns
  !
  
  ! Calculate the long range term by finite difference methods
  Allocate( pot_grid_ssp( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  Allocate( ei_ssp( 1:n ) )
  Allocate( force_ssp( 1:3, 1:n ) )
  Call ssp_long_range( l, q, r, alpha, FD_order, recip_E_ssp, q_grid, pot_grid_ssp, ei_ssp, force_ssp, t_grid, t_recip, error )
  Write( *, * ) 'Nett force ', Sum( force_ssp( 1, : ) ), Sum( force_ssp( 2, : ) ), Sum( force_ssp( 3, : ) )
  Open( 11, file = 'forces_ssp.dat' )
  Write( 11, * ) n, '     #number of particles'
  Write( 11, * ) Sum( force_ssp( 1, : ) ), Sum( force_ssp( 2, : ) ), Sum( force_ssp( 3, : ) ), '     #Nett force'
  Do i = 1, n
     Write( 11, * ) i, force_ssp( :, i )
  End Do
  Close( 11 )
  Write( *, * ) 'SSP grid  time: ', t_grid
  Write( *, * ) 'SSP solve time: ', t_recip
  l_bad_charge = .False.
  Select Case( error )
  Case( -1 )
     Write( *, * ) 'WARNING: Sum of charge on grid greater than tolerance but "fixed" (for some value of fixed)'
     l_bad_charge = .True.
  Case( -2 )
     Write( *, * ) 'WARNING: Sum of charge on grid greater than tolerance AND could not be fixed'
     l_bad_charge = .True.
     Stop
  End Select

  ! Save the SSP potential
  Call grid_io_save( 11, 'pot_grid_SSP.dat', l, pot_grid_ssp )

  Write( *, * ) 'SSP: Sum of charge over grid: ', Sum( q_grid )
  Write( *, * ) 'SSP: Sum over pot grid      : ', Sum( pot_grid_ssp )

  ! SIC ( Really same as ffp but let's be methodical ) 
  sic_ssp = ssp_sic( q, alpha )

  ! long range same as for ffp
  real_E_ssp = real_E_ffp

  ! And hence total energy 
  tot_E_ssp = real_E_ssp + sic_ssp + recip_E_ssp

  ! END SSP SECTION
  !
!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! REPORT RESULTS
  
  Write( *, '( t39, a, t149, a )' ) 'Internal', 'DL_POLY'
  Write( *, '( t14, a, t39, a, t64, a, t89, a, t124, a, t149, a, t174, a, t199, a  )' ) &
       'Reciprocal', 'SIC', 'Real', 'Total', 'Reciprocal', 'SIC', 'Real', 'Total' 
  Write( *, '( a4, 2( 4( g24.14, 1x ), :, 10x ) )' ) 'Ew  ', Real( recip_E, wp ), sic, real_E, tot_E, &
       Real( recip_E, wp ) * r4pie0, sic * r4pie0, real_E * r4pie0, tot_E * r4pie0
  If( do_sfp ) Then
     Write( *, '( a4, 2( 4( g24.14, 1x ), :, 10x ) )' ) &
          'SFP ', recip_E_sfp, sic_sfp, real_E_sfp, tot_E_sfp, &
          recip_E_sfp * r4pie0, sic_sfp * r4pie0, real_E_sfp * r4pie0, tot_E_sfp * r4pie0
  End If
  Write( *, '( a4, 2( 4( g24.14, 1x ), :, 10x ) )' ) &
       'FFP ', recip_E_ffp, sic_ffp, real_E_ffp, tot_E_ffp, &
       recip_E_ffp * r4pie0, sic_ffp * r4pie0, real_E_ffp * r4pie0, tot_E_ffp * r4pie0
  Write( *, '( a4, 2( 4( g24.14, 1x ), :, 10x ) )' ) &
       'SSP ', recip_E_ssp, sic_ssp, real_E_ssp, tot_E_ssp, &
       recip_E_ssp * r4pie0, sic_ssp * r4pie0, real_E_ssp * r4pie0, tot_E_ssp * r4pie0
  
  Write( *, * )
  If( do_sfp ) Then
     Write( *, '( "Difference in energy SFP - EW: ", g30.16 )' ) tot_E_sfp - tot_E
  End If
  Write( *, '( "Difference in energy FFP - EW:  ", g30.16 )' ) tot_E_ffp - tot_E
  Write( *, '( "Difference in energy SSP - EW:  ", g30.16 )' ) tot_E_ssp - tot_E
  Write( *, '( "Difference in energy FFP - SSP: ", g30.16 )' ) tot_E_ffp - tot_E_ssp

  rms_delta_pot   = Sum( ( pot_grid_ffp - pot_grid_ssp ) ** 2 )
  rms_delta_pot   = Sqrt( rms_delta_pot ) / Size( pot_grid_ffp )
  rms_delta_force = Sum( ( force_ffp - force_ssp ) ** 2 )
  rms_delta_force = Sqrt( rms_delta_force ) / Size( force_ffp )
  Write( *, '( a, g20.12 )' ) 'RMS (FFP-SSP) in potential: ', rms_delta_pot
  Write( *, '( a, g20.12 )' ) 'RMS (FFP-SSP) in force    : ', rms_delta_force

  ! Add a line to the summary file
  Inquire( file = 'summary.dat', exist = fexist )
  Open( 11, file = 'summary.dat', position = 'append' )
  If( .Not. fexist ) Then
     Write( 11, '( a7, 2x, a5, 2x, a2, 1x, a6, 1x, a2, 1x, a9, &
          & t40, 3( a14, 11x ), t120, 3( a14, 11x ), t200, a14, t230, a, t260, a )' ) &
          'alpha', 'dg', 'Od', 'xshift', 'Rg', 'Q error', &
          'tot E'  , 'tot E SSP'  , 'Delta tot E', &
          'recip E FFP', 'recip E SSP', 'Delta recip E', &
          'RMS delta pot', 'RMS Delta force', 'Bad Charge'
  End If
  Write( 11, '( f8.6, 1x, f6.4, 1x, i2, 1x, f6.4, 1x, i2, 1x, g9.2, &
       & t40, 3( g24.14, 1x ), t120, 3( g24.14, 1x ), t200, g24.14, t230, g24.14, t260, l1 )' ) &
       alpha, dg( 1 ), FD_order, xshift, range_gauss, q_error, &
       tot_E, tot_E_ssp, Abs( tot_E - tot_E_ssp ), &
       Real( recip_E_ffp, wp ), recip_E_ssp, Abs( Real( recip_E_ffp, wp ) - recip_E_ssp ), &
       rms_delta_pot, rms_delta_force, l_bad_charge
  Close( 11 )

  ! END REPORTING RESULTS SECTION
  !
!!!!!!!!!!!

Contains

  Subroutine read_header( unit, n, level, vecs )

    Integer                      , Intent( In    ) :: unit
    Integer                      , Intent(   Out ) :: n
    Integer                      , Intent(   Out ) :: level    
    Real( wp ), Dimension( :, : ), Intent(   Out ) :: vecs

    Integer :: junk
    
    Read( unit, * )
    Read( unit, * ) level, junk, n
    Read( unit, * ) vecs

  End Subroutine read_header

  Subroutine read_config( level, unit, q, r )

    Integer                      , Intent( In    ) :: level
    Integer                      , Intent( In    ) :: unit
    Real( wp ), Dimension( :    ), Intent(   Out ) :: q
    Real( wp ), Dimension( :, : ), Intent(   Out ) :: r

    Integer :: n
    Integer :: i

    Character( Len = 3 ) :: what
    
    n = Size( q )

    Do i = 1, n
       Read( unit, * ) what
       If( what == 'Na+' ) Then
          q( i ) = +1.0_wp
       Else If( what == 'Cl-' ) Then
          q( i ) = -1.0_wp
       Else
          Write( *, * ) 'what is a ', what, ' ?'
          Stop 'Error in CONFIG reading'
       End If
       Read( unit, * ) r( :, i )
       If( level >= 1 ) Then 
          Read( unit, * )
       End If
       If( level >= 2 ) Then 
          Read( unit, * )
       End If
    End Do
    
  End Subroutine read_config

  Subroutine generate_ew_func( l, q, r, alpha, ew_func )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Type( lattice )               , Intent( In    ) :: l
    Real( wp ), Dimension( :     ), Intent( In    ) :: q
    Real( wp ), Dimension( :, :  ), Intent( In    ) :: r
    Real( wp ),                     Intent( In    ) :: alpha
    Complex( wp ), Dimension( 0: ), Intent(   Out ) :: ew_func

    Complex( wp ) :: s_fac

    Real( wp ), Dimension( 1:3 ) :: G

    Real( wp ) :: G_len_sq, G_fac
    
    Integer :: iG
    
    ew_func( 0 ) = 0.0_wp
    !$omp parallel default( none ) shared( n, ew_func, l, alpha, q, r ) &
    !$omp                          private( s_fac, ig, G, G_len_sq, G_fac )
    !$omp do
    Do iG = 1, Ubound( ew_func, Dim = 1 )
       Call l%get_nth_rec_vec( ig + 1, G )
       G = G * 2.0_wp * pi
       ! Calculate the structure factor for this G vector
       s_fac = 0.0_wp
       Do j = 1, n
          s_fac = s_fac + &
               q( j ) * Exp( Cmplx( 0.0_wp, - Dot_product( G, r( :, j ) ), wp ) )
       End Do
       G_len_sq = Dot_product( G, G )
       G_fac = Exp( - G_len_sq / ( 4.0_wp * alpha * alpha ) ) / G_len_sq
       ! Note difference from Darden paper - here in the ewfunc we have adsorbed 2*pi into the lattice vectors
       ! Hence we multiply though by 4 * pi * pi and then another factor of 2 as this is for the potential
       ! not the energy. This formulation is the same as in Appendix B of Kittel
       ew_func( iG ) = ( 4.0_wp * pi / l%get_volume() ) * s_fac * g_fac
    End Do
    !$omp end do
    !$omp end parallel

  End Subroutine generate_ew_func

  Subroutine trad_ewald( l, q, r, alpha, ew_func, recip_E, sic, real_E, tot_E, t_recip, t_real )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

    Use lattice_module, Only : lattice

    Type( lattice )               , Intent( In    ) :: l
    Real( wp ), Dimension( :     ), Intent( In    ) :: q
    Real( wp ), Dimension( :, :  ), Intent( In    ) :: r
    Real( wp ),                     Intent( In    ) :: alpha
    Complex( wp ), Dimension( 0: ), Intent( In    ) :: ew_func
    Real( wp )                    , Intent( Out   ) :: recip_E
    Real( wp )                    , Intent( Out   ) :: SIC
    Real( wp )                    , Intent( Out   ) :: real_E
    Real( wp )                    , Intent( Out   ) :: tot_E
    Real( wp )                    , Intent( Out   ) :: t_recip
    Real( wp )                    , Intent( Out   ) :: t_real

    Complex( wp ) :: potg
    Complex( wp ) :: cc, cy, ct

    Real( wp ), Dimension( 1:3 ) :: G
    Real( wp ), Dimension( 1:3 ) :: ri
    
    Real( wp ) :: pot
    
    Integer :: max_G_shells = 2
    Integer :: n
    Integer :: i, iG
    
    Integer( li ) :: start, finish, rate

    n = Size( q )
    
    Call system_clock( start, rate )
    recip_E = 0.0_wp
    !$omp parallel default( none ) shared( l, r, recip_E, q, n, ew_func ) &
    !$omp                          private( i, iG, pot, G, ri, cc, cy, ct, potg )
    cc = 0.0_wp
    !$omp do reduction( +:recip_E )
    Do i = 1, n
       pot = 0.0_wp
       ri = r( :, i )
       Do iG = 1, Ubound( ew_func, Dim = 1 )
          Call l%get_nth_rec_vec( ig + 1, G )
          G = G * 2.0_wp * pi
          potG = Exp( Cmplx( 0.0_wp, Dot_product(   G, ri ), wp ) ) *        ew_func( ig )
          pot = pot + Real( potG + Conjg( potG ), wp )
       End Do
       cy = q( i ) * pot - cc
       ct = recip_E + cy
       cc = ( ct - recip_E ) - cy
       recip_E = Real( ct, wp )
    End Do
    !$omp end do
    !$omp end parallel
    recip_E = recip_E * 0.5_wp
    Call system_clock( finish, rate )
    t_recip = Real( finish - start, wp ) / rate

    ! Self interaction correction
    sic = - Sum( q * q ) * alpha / Sqrt( pi )

    Call system_clock( start, rate )
    !$omp parallel default( none ) shared( l, q, r, alpha, max_G_shells, real_E )
    Call real_space_energy( l, q, r, alpha, max_G_shells, real_E )
    !$omp end parallel
    Call system_clock( finish, rate )
    t_real = Real( finish - start, wp ) / rate

    tot_E = real_E + Real( recip_E, wp ) + sic

  End Subroutine trad_ewald

  Subroutine get_n_grid( l, alpha, range_gauss, gauss_tol, n_grid )

    Type( lattice ),              Intent( In    ) :: l
    Real( wp ),                   Intent( In    ) :: alpha
    Integer   ,                   Intent( In    ) :: range_gauss
    Real( wp ),                   Intent( In    ) :: gauss_tol
    Integer   , Dimension( 1:3 ), Intent(   Out ) :: n_grid

    Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp

    Real( wp ), Dimension( 1:3, 1:3 ) :: dir_vecs
    
    Real( wp ), Dimension( 1:3 ) :: dr

    Real( wp ) :: Gsq
    
    Real( wp ) :: max_r

    Integer :: i

    ! Find range after which the NORMALISED gaussian decays to below the tolerance
    max_r = pi * Sqrt( pi ) * gauss_tol / ( alpha * alpha * alpha )
    max_r = - Log( max_r ) / ( alpha * alpha )
    max_r = Sqrt( max_r )

    ! Now work out what resolution that require given the number of grid points
    ! we are using to represent the gaussian
    dr = max_r / range_gauss

    ! And so how many grid points along each lattice vector, rounding up
    ! to ensure accuracy
    dir_vecs = l%get_direct_vectors()

    Do i = 1, 3
       Gsq = Dot_product( dir_vecs( :, i ), dir_vecs( :, i ) )
       n_grid( i ) = Ceiling( Sqrt( Gsq ) / dr( i ) )
    End Do

  End Subroutine get_n_grid
  
End Program test


