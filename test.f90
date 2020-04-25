Program test

  !$ Use omp_lib
  
  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

  Use lattice_module, Only : lattice
  Use FD_Laplacian_3d_module, Only : FD_Laplacian_3D
  Use minresmodule, Only : minres
  Use fft_module, Only : fft_fft3d
  
  Implicit None

  Type( lattice ) :: l
  Type( FD_Laplacian_3d ) :: FD

  Complex( wp ), Dimension( :, :, : ), Allocatable :: struc_fac
  Complex( wp ), Dimension( :, :, : ), Allocatable :: pot_k

  Complex( wp ), Dimension( : ), Allocatable :: ew_func

  Complex( wp ) :: s_fac
  Complex( wp ) :: recip_E, pot, potG
  Complex( wp ) :: cc, cy, ct
  
  Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp

  Real( wp ), Parameter :: r4pie0 = 138935.48350000000_wp
  
  Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_fd
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_ffp
  
  Real( wp ), Dimension( :, : ), Allocatable :: a
  Real( wp ), Dimension( :, : ), Allocatable :: r

  Real( wp ), Dimension( : ), Allocatable :: q

  Real( wp ), Dimension( 1:3, 1:3 ) :: dGrid_vecs

  Real( wp ), Dimension( 1:3 ) :: t
  Real( wp ), Dimension( 1:3 ) :: G, dG
  Real( wp ), Dimension( 1:3 ) :: ri, fi
  Real( wp ), Dimension( 1:3 ) :: Grid_vec
  Real( wp ), Dimension( 1:3 ) :: f_point, r_point

  Real( wp ) :: alpha, q_norm, qi_norm
  Real( wp ) :: q_val, potr
  Real( wp ) :: G_len_sq, G_fac, term
  Real( wp ) :: sic
  Real( wp ) :: real_E, tot_E
  Real( wp ) :: recip_E_ffp, sic_ffp, real_E_ffp, tot_E_ffp
  Real( wp ) :: recip_E_ffp_fd, tot_E_ffp_fd
  Real( wp ) :: recip_E_sfp, tot_E_sfp
  Real( wp ) :: Anorm, Arnorm, Acond, rnorm, ynorm, rtol
  Real( wp ) :: xshift
  Real( wp ) :: rms_delta_pot, q_error
  
  Integer, Dimension( 1:3 ) :: n_grid, i_atom_centre, i_atom_grid, i_point, i_grid
  Integer, Dimension( 1:3 ) :: n_vec
  Integer, Dimension( 1:3 ) :: range_gauss
  Integer, Dimension( 1:3 ) :: iG_vec

  Integer :: FD_order
  Integer :: n, level
  Integer :: nd
  Integer :: i, j, iG, i_dir
  Integer :: max_G_shells = 2
  Integer :: istop
  Integer :: itn
  Integer :: i1, i2, i3
  
  Integer( li ) :: start, finish, rate

  Logical :: fexist
  
  Interface
     Subroutine dummy_Msolve(lb,ub,x,y)                   ! Solve M*y = x
       Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
       Integer,  Intent(in)    :: lb( 1:3 ), ub( 1:3 )
       Real(wp), Intent(in)    :: x( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) )
       Real(wp), Intent(out)   :: y( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) )
     End Subroutine Dummy_Msolve
  End Interface
  
  Write( *, * ) 'Ewald param ?'
  Read ( *, * ) alpha
  Write( *, * ) 'Grid Dims ?'
  Read ( *, * ) n_grid
  Write( *, * ) 'FD_order?'
  Read ( *, * ) FD_order
  Write( *, * ) 'xshift?'
  Read ( *, * ) xshift

  Write( *, * ) 'Param, dims, order ', alpha, n_grid, FD_order, xshift

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
!!$  Call l%print

  Write( *, * ) 'r(1) before shift and reference = ', r( :, 1 )
  Do i = 1, n
     t = r( :, i )
     t( 1 ) = t( 1 ) + xshift
     Call l%to_reference( t, r( :, i ) ) 
  End Do
  Write( *, * ) 'r(1) after shift and reference = ', r( :, 1 )

  Call system_clock( start, rate )
  Allocate( ew_func( 0:l%get_n_rec_vecs() - 1 ) )
  Write( *, * ) '!!!!!!!!!!!!!!!!!!!!!!!!', l%get_rec_vec_max_index(), l%get_n_rec_vecs()
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
  Call system_clock( finish, rate )
  Write( *, * ) 'ewfunc time ', Real( finish - start, wp ) / rate

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
        pot = pot + potG + Conjg( potG )
     End Do
     cy = q( i ) * pot - cc
     ct = recip_E + cy
     cc = ( ct - recip_E ) - cy
     recip_E = ct
  End Do
  !$omp end do
  !$omp end parallel
  recip_E = recip_E * 0.5_wp
  sic = - Sum( q * q ) * alpha / Sqrt( pi )
  Call system_clock( finish, rate )
  Write( *, * ) 'recip_E time ', Real( finish - start, wp ) / rate

  Call system_clock( start, rate )
  !$omp parallel default( none ) shared( l, q, r, alpha, max_G_shells, real_E )
  Call real_space_energy( l, q, r, alpha, max_G_shells, real_E )
  !$omp end parallel
  Call system_clock( finish, rate )
  Call system_clock( finish, rate )
  Write( *, * ) 'real_E time ', Real( finish - start, wp ) / rate

  tot_E = real_E + Real( recip_E, wp ) + sic
        
  Call system_clock( start, rate )
  Allocate( q_grid( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  ! Find range of gaussian
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
  Write( *, * ) 'range_gauss = ', range_gauss
  Call grid_charge( l, alpha, q, r, range_gauss, q_grid )
  Call system_clock( finish, rate )
  Write( *, * ) 'q grid time ', Real( finish - start, wp ) / rate, ( Real( finish - start, wp ) / rate ) / Product( n_grid )
  q_error = Sum( q_grid )
  Write( *, * ) 'Sum of charge over grid ', Sum( q_grid )
  Open ( 11, file = 'q_grid.dat' )
  Write( 11, * ) 'Sum of charge over grid ', Sum( q_grid ), Sum( q_grid ) / l%get_volume()
  Do i3 = 0, n_grid( 3 ) - 1
     Do i2 = 0, n_grid( 2 ) - 1
        Do i1 = 0, n_grid( 1 ) - 1
           i_grid = [ i1, i2, i3 ]
           f_point = Real( i_grid, wp ) / n_grid
           Call l%to_direct( f_point, r_point )
           Write( 11, '( 3( f6.2, 1x ), 5x, f15.12 )' ) r_point, q_grid( i1, i2, i3 )
        End Do
     End Do
  End Do
  Close( 11 )
  
  Allocate( pot_grid( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
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
  Write( *, * ) 'pot grid time ', Real( finish - start, wp ) / rate, ( Real( finish - start, wp ) / rate ) / Product( n_grid )
  Write( *, * ) 'Sum over pot grid ', Sum( pot_grid ), Sum( pot_grid ) / l%get_volume()
  Open ( 11, file = 'pot_grid.dat' )
  Write( 11, * ) 'Sum over pot grid ', Sum( pot_grid ), Sum( pot_grid ) / l%get_volume()
  Do i3 = 0, n_grid( 3 ) - 1
     Do i2 = 0, n_grid( 2 ) - 1
        Do i1 = 0, n_grid( 1 ) - 1
           i_grid = [ i1, i2, i3 ]
           f_point = Real( i_grid, wp ) / n_grid
           Call l%to_direct( f_point, r_point )
           Write( 11, '( 3( f6.2, 1x ), 5x, f15.12, 1x, f15.8 )' ) &
                r_point, pot_grid( i1, i2, i3 ), pot_grid( i1, i2, i3 ) * r4pie0
        End Do
     End Do
  End Do
  Close( 11 )

  ! Minus Sign on charge grid as we integrate against the screening charge, which
  ! has opposite sign to the actual charge (as it is screening, Duh!
  ! Minus sign on whole comes here just so can add up all contribs at end rather than confuse myself about signs
  recip_E_sfp = - 0.5_wp * Sum( - q_grid * pot_grid ) * ( l%get_volume() / Product( n_grid ) )

  sic_ffp = - alpha * Sum( q * q ) / Sqrt( 2.0_wp * pi )
  
  Call system_clock( start, rate )
  !$omp parallel default( none ) shared( l, q, r, alpha, max_G_shells, real_E_ffp )
  Call real_space_energy( l, q, r, alpha / Sqrt( 2.0_wp ), max_G_shells, real_E_ffp )
  !$omp end parallel
  Call system_clock( finish, rate )
  Write( *, * ) 'real_E_ffp time ', Real( finish - start, wp ) / rate

  tot_E_sfp = real_E_ffp + sic_ffp + recip_E_sfp

  ! Now the big one
  dGrid_vecs = l%get_direct_vectors()
  Do i = 1, 3
     dGrid_vecs( :, i ) = dGrid_vecs( :, i ) / n_grid( i )
     dG = Sqrt( Dot_product( dGrid_vecs( :, i ), dGrid_vecs( :, i ) ) )
  End Do
  Call FD%init( FD_order, dGrid_vecs )

  Allocate( pot_grid_fd( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  rtol = 1.0e-12_wp
  Call system_clock( start, rate )
  Call minres( Lbound( q_grid ), Ubound( q_grid ), FD, dummy_Msolve, - 4.0_wp * pi * q_grid, 0.0_wp, .True., .False., &
       pot_grid_fd, 1000, 99, rtol,                      &
       istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm )
  pot_grid_fd = pot_grid_fd - Sum( pot_grid_fd ) / Size( pot_grid_fd )
  Call system_clock( finish, rate )
  Write( *, * ) 'MINRES time ', Real( finish - start, wp ) / rate

  Write( *, * ) 'Iterative solver summary:'
  Write( *, * ) 'alpha                = ', alpha
  Write( *, * ) 'Grid resolution      = ', dG
  Write( *, * ) 'Grid size            = ', n_grid
  Write( *, * ) 'Order                = ', FD_order
  Write( *, * ) 'iterations           = ', itn
  Write( *, * ) 'norm of the residual = ', rnorm
  
  Open ( 11, file = 'pot_grid_fd.dat' )
  Write( 11, * ) 'Sum over pot grid ', Sum( pot_grid_fd ), Sum( pot_grid_fd ) / l%get_volume()
  Do i3 = 0, n_grid( 3 ) - 1
     Do i2 = 0, n_grid( 2 ) - 1
        Do i1 = 0, n_grid( 1 ) - 1
           i_grid = [ i1, i2, i3 ]
           f_point = Real( i_grid, wp ) / n_grid
           Call l%to_direct( f_point, r_point )
           Write( 11, '( 3( f6.2, 1x ), 5x, f15.12, 1x, f15.8 )' ) &
                r_point, pot_grid_fd( i1, i2, i3 ), pot_grid_fd( i1, i2, i3 ) * r4pie0
        End Do
     End Do
  End Do
  Close( 11 )

  recip_E_ffp_fd = - 0.5_wp * Sum( - q_grid * pot_grid_fd ) * ( l%get_volume() / Product( n_grid ) )
  tot_E_ffp_fd = real_E_ffp + sic_ffp + recip_E_ffp_fd
  
  Write( *, '( t39, a, t149, a )' ) 'Internal', 'DL_POLY'
  Write( *, '( t14, a, t39, a, t64, a, t89, a, t124, a, t149, a, t174, a, t199, a  )' ) &
       'Reciprocal', 'SIC', 'Real', 'Total', 'Reciprocal', 'SIC', 'Real', 'Total' 
  Write( *, '( a4, 2( 4( g24.14, 1x ), :, 10x ) )' ) 'Ew  ', Real( recip_E, wp ), sic, real_E, tot_E, &
       Real( recip_E, wp ) * r4pie0, sic * r4pie0, real_E * r4pie0, tot_E * r4pie0
  Write( *, '( a4, 2( 4( g24.14, 1x ), :, 10x ) )' ) &
       'SFP ', recip_E_sfp, sic_ffp, real_E_ffp, tot_E_sfp, &
       recip_E_sfp * r4pie0, sic_ffp * r4pie0, real_E_ffp * r4pie0, tot_E_sfp * r4pie0
  Write( *, '( a4, 2( 4( g24.14, 1x ), :, 10x ) )' ) &
       'SSP ', recip_E_ffp_fd, sic_ffp, real_E_ffp, tot_E_ffp_fd, &
       recip_E_ffp_fd * r4pie0, sic_ffp * r4pie0, real_E_ffp * r4pie0, tot_E_ffp_fd * r4pie0
  
  Write( *, * )
  Write( *, '( "Difference in energy SFP - EW: ", g30.16 )' ) tot_E_sfp    - tot_E
  Write( *, '( "Difference in energy SSP - EW: ", g30.16 )' ) tot_E_ffp_fd - tot_E

  rms_delta_pot = Sum( ( pot_grid - pot_grid_fd ) ** 2 )
  rms_delta_pot = Sqrt( rms_delta_pot ) / Size( pot_grid )
  Inquire( file = 'summary.dat', exist = fexist )
  Open( 11, file = 'summary.dat', position = 'append' )
  If( .Not. fexist ) Then
     Write( 11, '( a7, 2x, a5, 2x, a2, 1x, a6, 1x, a2, 1x, a9, &
          & t40, 3( a14, 11x ), t120, 3( a14, 11x ), t200, a14 )' ) &
          'alpha', 'dg', 'Od', 'xshift', 'Rg', 'Q error', &
          'tot E'  , 'tot E SSP'  , 'Delta tot E', &
          'recip E FFP', 'recip E SSP', 'Delta recip E', &
          'RMS delta pot'
  End If
  Write( 11, '( f8.6, 1x, f6.4, 1x, i2, 1x, f6.4, 1x, i2, 1x, g9.2, &
       & t40, 3( g24.14, 1x ), t120, 3( g24.14, 1x ), t200, g24.14 )' ) &
       alpha, dg( 1 ), FD_order, xshift, range_gauss( 1 ), q_error, &
       tot_E, tot_E_ffp_fd, Abs( tot_E - tot_E_ffp_fd ), &
       Real( recip_E_sfp, wp ), recip_E_ffp_fd, Abs( Real( recip_E_sfp, wp ) - recip_E_ffp_fd ), &
       rms_delta_pot
  Close( 11 )

  Allocate( struc_fac( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  Allocate( pot_k( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  Allocate( pot_grid_ffp( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  struc_fac = q_grid
  Call fft_fft3d( -1, struc_fac )
  recip_E_ffp = 0.0_wp
  Do i3 = 0, n_grid( 3 ) - 1
     Do i2 = 0, n_grid( 2 ) - 1
        Do i1 = 0, n_grid( 1 ) - 1
           i_point = [ i1, i2, i3 ]
           Where( i_point < n_grid / 2 )
              iG_vec = i_point
           Else Where
              iG_vec = i_point - n_grid
           End Where
           Call l%get_rec_vec( iG_vec, G )
           G = G * 2.0_wp * pi
           G_len_sq = Dot_product( G, G )
           If( G_len_sq > 1e-12_wp ) Then
              G_fac = Exp( - G_len_sq / ( 4.0_wp * alpha * alpha ) ) / G_len_sq
              pot_k( i1, i2, i3 ) = struc_fac( i1, i2, i3 ) * G_fac * ( 4.0_wp * pi / l%get_volume() )
              recip_E_ffp = recip_E_ffp + pot_k( i1, i2, i3 ) * Conjg( struc_fac( i1, i2, i3 ) )
           Else
              pot_k( i1, i2, i3 ) = 0.0_wp
           End If
        End Do
     End Do
  End Do
  Write( *, * ) recip_E_ffp 
  Write( *, * ) recip_E_ffp / ( 2.0_wp * pi * l%get_volume() )
  ! Hacky extra value of pi, why?
  Write( *, * ) recip_E_ffp / ( 2.0_wp * pi * Product( n_grid ) )
  Call fft_fft3d( 1, pot_k )
!!$  pot_grid_ffp = Real( pot_k, wp ) / ( 2.0_wp * l%get_volume() )
!!$  pot_grid_ffp = Real( pot_k, wp ) / ( 2.0_wp * Product( n_grid ) )
  pot_grid_ffp = Real( pot_k, wp )
  Write( *, * )  pot_grid( 1, 1, 1 ), pot_grid_ffp( 1, 1, 1 )
  recip_E_ffp = - 0.5_wp * Sum( - q_grid * pot_grid_ffp ) * ( l%get_volume() / Product( n_grid ) )
  tot_E_ffp = real_E_ffp + sic_ffp + recip_E_ffp
  Write( *, '( a4, 2( 4( g24.14, 1x ), :, 10x ) )' ) &
       'FFP ', recip_E_ffp, sic_ffp, real_E_ffp, tot_E_ffp, &
       recip_E_ffp * r4pie0, sic_ffp * r4pie0, real_E_ffp * r4pie0, tot_E_ffp * r4pie0

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

  Subroutine real_space_energy( l, q, r, alpha, n_G_shells, E )

    Type( lattice )              , Intent( In    ) :: l
    Real( wp ), Dimension( :    ), Intent( In    ) :: q
    Real( wp ), Dimension( :, : ), Intent( In    ) :: r
    Real( wp ),                    Intent( In    ) :: alpha
    Integer   ,                    Intent( In    ) :: n_G_shells
    Real( wp ),                    Intent(   Out ) :: E

    Real( wp ), Dimension( 1:3 ) :: G
    Real( wp ), Dimension( 1:3 ) :: ri, rjG
    Real( wp ), Dimension( 1:3 ) :: rj, rij

    Real( wp ) :: rijG, rijG_sq
    Real( wp ) :: c, y, t

    Integer :: n
    Integer :: side, max_G_vector
    Integer :: i, j, iG

    n = Size( q )
    
    ! Implicit barrier at end of single makes sure shared E is inited by all threds
    !$omp single
    E = 0.0_wp
    !$omp end single
    
    ! Assume cubic cell this is the number of G vectors we need go out to reach the
    ! number of shells required excepting the reference cell
    ! CHECK THIS!! Pretty sure it is right but worth going over again

    Select Case( n_G_shells )

    Case( 1: )
       ! Include lattice sum over the number of shells of G vectors
       side =  2 * n_G_shells + 1
       max_G_vector = n_G_shells * side * side + ( ( side + 1 ) * side ) / 2 - 1
       c = 0.0_wp
       !$omp do reduction( +:E )
       Do i = 1, n
          ri = r( :, i )
          Do j = 1, n
             If( i /= j ) Then
                ! Zero G vector
                rjG = r( :, j )
                rijG_sq = Dot_product( ri - rjG, ri - rjG )
                rijG = Sqrt( rijG_sq )
                y = q( i ) * q( j ) * Erfc( alpha * rijG ) / rijG - c
                t = E + y
                c = ( t - E ) - y
                E = t
             End If
             Do iG = 1, max_G_vector
                !+ve G vector
                Call l%get_nth_dir_vec( ig + 1, G )
                rjG = r( :, j ) + G
                rijG_sq = Dot_product( ri - rjG, ri - rjG )
                rijG = Sqrt( rijG_sq )
                y = q( i ) * q( j ) * Erfc( alpha * rijG ) / rijG - c
                t = E + y
                c = ( t - E ) - y
                E = t
                !-ve G vector
                rjG = r( :, j ) - G
                rijG_sq = Dot_product( ri - rjG, ri - rjG )
                rijG = Sqrt( rijG_sq )
                y = q( i ) * q( j ) * Erfc( alpha * rijG ) / rijG - c
                t = E + y
                c = ( t - E ) - y
                E = t
             End Do
          End Do
       End Do
       !$omp end do

    Case( 0 )
       ! Just minimum image convention
       c = 0.0_wp
       !$omp do reduction( +:E )
       Do i = 1, n
          ri = r( :, i )
          Do j = 1, n
             If( i /= j ) Then
                rj = r( :, j )
                ! Apply minimum image convention
                Call l%minimum_image( ri - rj, rij )
                rijG_sq = Dot_product( rij, rij ) 
                rijG = Sqrt( rijG_sq )
                ! Should implement rcut at some point
!!$                If( rijG <= 11.0_wp ) Then
                   y = q( i ) * q( j ) * Erfc( alpha * rijG ) / rijG - c
                   t = E + y
                   c = ( t - E ) - y
                   E = t
!!$                End If
             End If
          End Do
       End Do
       !$omp end do

    Case Default
       Stop "Negative number of shells in real space lattice summation"
       
    End Select

    ! Make sure scaling is just applied by one thread
    !$omp single
    E = E * 0.5_wp
    !$omp end single

  End Subroutine real_space_energy

  Subroutine grid_charge( l, alpha, q, r, range_gauss, q_grid )

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
    
  End Subroutine grid_charge

End Program test

  
Subroutine dummy_msolve(lb,ub,x,y)                   ! Solve M*y = x
  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
  Integer,  Intent(in)    :: lb( 1:3 ), ub( 1:3 )
  Real(wp), Intent(in)    :: x( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) )
  Real(wp), Intent(out)   :: y( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) )
End Subroutine dummy_msolve

