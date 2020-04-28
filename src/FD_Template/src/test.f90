Program testit
  !!----------------------------------------------------
  !! Program to test the execution of the finite difference
  !! library by calculating the charge from a screened
  !! Gaussian according to
  !! $$ \rho _i^s(\mathbf r)=-q_i\left(\frac{\alpha }{\sqrt{\pi }}\right)^3e^{-\alpha^2 |\mathbf r-\mathbf r_i|^2} $$
  !!
  !! Written by I.J. Bush
  !!----------------------------------------------------
  
  Use numbers_module        , Only : wp
  Use FD_Laplacian_3d_module, Only : FD_Laplacian_3D

  Implicit None

  Type( FD_Laplacian_3d ) :: FD

  Real( wp ), Dimension( :, :, : ), Allocatable :: grid
  Real( wp ), Dimension( :, :, : ), Allocatable :: laplacian
  Real( wp ), Dimension( :, :, : ), Allocatable :: fd_laplacian

  Real( wp ), Dimension( 1:3, 1:3 ) :: grid_vecs

  Real( wp ), Dimension( 1:3 ) :: r, ri

  Real( wp ) :: alpha, alpha_sq
  Real( wp ) :: x, y, z
  Real( wp ) :: norm
  Real( wp ) :: arg
  Real( wp ) :: gauss, gauss_x2, gauss_y2, gauss_z2
  Real( wp ) :: max_error, rel_error, rel_value, max_rel_error

  Integer, Dimension( 1:3 ) :: grid_size, grid_with_halo

  Integer :: order
  Integer :: i3, i2, i1
  Integer :: start, finish, rate

  Write( *, * ) 'Grid vecs?'
  Read ( *, * ) grid_vecs

  Write( *, * ) 'ri ?'
  Read ( *, * ) ri

  Write( *, * ) 'alpha ?'
  Read ( *, * ) alpha
  alpha_sq = alpha * alpha

  Write( *, * ) 'order ?'
  Read ( *, * ) order

  Write( *, * ) 'grid_size ?'
  Read ( *, * ) grid_size

  grid_with_halo = grid_size + order / 2

  Allocate( grid( &
    & -grid_with_halo( 1 ):grid_with_halo( 1 ), &
    & -grid_with_halo( 2 ):grid_with_halo( 2 ), &
    & -grid_with_halo( 3 ):grid_with_halo( 3 ) ) )
  Allocate( laplacian( &
    & -grid_with_halo( 1 ):grid_with_halo( 1 ), &
    & -grid_with_halo( 2 ):grid_with_halo( 2 ), &
    & -grid_with_halo( 3 ):grid_with_halo( 3 ) ) )
  Allocate( fd_laplacian( &
    & -grid_size( 1 ):grid_size( 1 ), &
    & -grid_size( 2 ):grid_size( 2 ), &
    & -grid_size( 3 ):grid_size( 3 ) ) )

  Call system_clock( start, rate )
  norm = ( alpha / Sqrt( 3.1415926535897932384626433832795_wp ) ) ** 3
  Do i3 = -grid_with_halo( 3 ), grid_with_halo( 3 )
    Do i2 = -grid_with_halo( 2 ), grid_with_halo( 2 )
      Do i1 = -grid_with_halo( 1 ), grid_with_halo( 1 )
        r = i1 * grid_vecs( :, 1 ) + i2 * grid_vecs( :, 2 ) + i3 * grid_vecs( :, 3 ) - ri
        arg = alpha_sq * Dot_Product( r, r )
        gauss = norm * Exp( - arg )
        grid( i1, i2, i3 ) = gauss
        x = r( 1 )
        y = r( 2 )
        z = r( 3 )
        gauss_x2 = gauss * ( - 2.0_wp * alpha_sq + 4.0 * alpha_sq * alpha_sq * x * x )
        gauss_y2 = gauss * ( - 2.0_wp * alpha_sq + 4.0 * alpha_sq * alpha_sq * y * y )
        gauss_z2 = gauss * ( - 2.0_wp * alpha_sq + 4.0 * alpha_sq * alpha_sq * z * z )
        laplacian( i1, i2, i3 ) = gauss_x2 + gauss_y2 + gauss_z2
      End Do
    End Do
  End Do
  Call system_clock( finish, rate )
  Write( *, '( a, t30, f8.3 )' ) 'Time for griding charge', Real( finish - start ) / rate

  Call FD%init( order, grid_vecs )
  Call system_clock( start, rate )
  Call FD%apply( -grid_with_halo, -grid_size, -grid_size, grid_size, grid, fd_laplacian )
  Call system_clock( finish, rate )
  Write( *, '( a, t30, f8.3 )' ) 'Time for FD laplacian', Real( finish - start ) / rate

  max_error =  Maxval( Abs( fd_laplacian - &
    laplacian( &
    & -grid_size( 1 ):grid_size( 1 ), &
    & -grid_size( 2 ):grid_size( 2 ), &
    & -grid_size( 3 ):grid_size( 3 ) ) ) )
  Write( *, '( a, t40, ": ", g22.16 )' ) 'Maximum error in laplacian', max_error
  max_rel_error = -1.0_wp
  rel_error = max_rel_error
  Do i3 = - grid_size( 3 ), grid_size( 3 )
     Do i2 = - grid_size( 2 ), grid_size( 2 )
        Do i1 = - grid_size( 1 ), grid_size( 1 )
           If( Abs( laplacian( i3, i2, i1 ) ) > 1.0e-12_wp ) Then
              rel_error = Abs( fd_laplacian( i1, i2, i3 ) - laplacian( i1, i2, i3 ) )
              rel_error = rel_error / Abs( laplacian( i1, i2, i3 ) )
              If( rel_error > max_rel_error ) Then
                 max_rel_error = rel_error
                 rel_value = laplacian( i1, i2, i3 )
              End If
           End If
        End Do
     End Do
  End Do
  Write( *, '( a, t40, ": ", g22.16, a, g22.16 )' ) &
       'Maximum relative error in laplacian', max_rel_error, ' where laplacian = ', rel_value
  Write( *, '( a, t40, ": ", g22.16 )' ) 'Maximum value of laplacian ', &
       Maxval( Abs( laplacian( &
    & -grid_size( 1 ):grid_size( 1 ), &
    & -grid_size( 2 ):grid_size( 2 ), &
    & -grid_size( 3 ):grid_size( 3 ) ) ) )

End Program testit
