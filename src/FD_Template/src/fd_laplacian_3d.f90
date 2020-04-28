Module FD_Laplacian_3d_module

  Use numbers_module,     Only : wp
  Use FD_template_module, Only : FD_template

  Implicit None

  ! Default for orthogonality: If the dot product is less than this vecs are deemed orthogonal
  Real( wp ), Parameter :: orthog_tol_default = 1.0e-14_wp 
                                                 
  ! Default number of reals in cache
  Integer, Parameter :: n_cache_default = ( 2 ** 18 ) / ( 8 ) 

  Integer, Parameter, Private :: XX = 1
  Integer, Parameter, Private :: XY = 2
  Integer, Parameter, Private :: XZ = 3
  Integer, Parameter, Private :: YY = 4
  Integer, Parameter, Private :: YZ = 5
  Integer, Parameter, Private :: ZZ = 6

  Type, Extends( FD_template ), Public :: FD_Laplacian_3d
     Integer                     , Private :: n_cache    = n_cache_default
     Real( wp )                  , Private :: orthog_tol = orthog_tol_default 
     Integer   , Dimension( 1:3 ), Private :: nc_block_apply
     Integer   , Dimension( 1:3 ), Private :: nc_block_jacobi
     Real( wp ), Dimension( 1:6 ), Private :: deriv_weights
     Real( wp )                  , Private :: diag_inv
     Logical                     , Private :: need_XY     ! Indicate if do  not need
     Logical                     , Private :: need_XZ     ! these compinents due to
     Logical                     , Private :: need_YZ     ! grid orthogonality
   Contains
     Procedure, Public :: init             ! With a bit of thought could probably
     Procedure, Public :: reset_vecs       ! move these two into the lower level template

     Procedure, Public :: apply            ! Should have deferred versions of
     Procedure, Public :: jacobi_sweep     ! these two in the lower level template
  End type FD_Laplacian_3d

  Private

Contains

  Subroutine init( FD, order, vecs, n_cache, orthog_tol )
    !!----------------------------------------------------
    !! Initialise Laplacian differentiator and calculate
    !! optimal blocking in cache.
    !!
    !! Written by I.J. Bush
    !!----------------------------------------------------
    Class( FD_Laplacian_3d )     , Intent(   Out ) :: FD
    Integer                      , Intent( In    ) :: order
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs
    Integer   , Optional         , Intent( In    ) :: n_cache
    Real( wp ), Optional         , Intent( In    ) :: orthog_tol

    ! Set the order and vectors associated witht the basic template
    Call FD%FD_init( 2, order / 2, vecs )

    ! If given set the estimated size of cache in terms of the number of reals
    ! it can gold
    If( Present( n_cache ) ) Then
       FD%n_cache = n_cache
    End If

    ! If present set the orthogonality tolerance
    If( Present( orthog_tol ) ) Then
       FD%orthog_tol = orthog_tol
    End If

    ! Calculate some cache blocking parameters
    FD%nc_block_apply  = get_cache_block_facs( FD%n_cache, 1, 1, 2 * FD%get_order() )
    FD%nc_block_jacobi = get_cache_block_facs( FD%n_cache, 1, 2, 2 * FD%get_order() )
    
    ! Set the weights for the derivatives on the different directions
    Call FD%reset_vecs( vecs )

  Contains

    Pure Function get_cache_block_facs( n_cache, n_grids, n_fd, accuracy ) Result( nb )

      ! Calculate some cache blocking parameters

      Integer, Dimension( 1:3 ) :: nb

      Integer, Intent( In ) :: n_cache
      Integer, Intent( In ) :: n_grids
      Integer, Intent( In ) :: n_fd
      Integer, Intent( In ) :: accuracy

      nb = 1

      Do While( usage( n_grids, n_fd, nb, accuracy ) < n_cache )
         nb = nb + 1
      End Do

      nb( 1 ) = nb( 1 ) + 1

      ! Use a strictly greater than relationship juts to allow
      ! a little more leeway - the estimate cahce usage is probably slightly low ...
      If( usage( n_grids, n_fd, nb, accuracy ) > n_cache ) Then
         nb( 1 ) = nb( 1 ) - 1

      Else
         nb( 2 ) = nb( 2 ) + 1

         If( usage( n_grids, n_fd, nb, accuracy ) > n_cache ) Then
            nb( 2 ) = nb( 2 ) - 1
         End If
      End If

    End Function get_cache_block_facs

  End Subroutine init

  Subroutine reset_vecs( FD, vecs )
    !!-----------------------------------------------------------
    !! Calculate weights due to grid offset for non-orthogonal grids
    !! This is separated from the init as this is the onlt part that needs
    !! to be called if the vectors describing the grid change - hence the name
    !! reset_vecs
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------
    Class( FD_Laplacian_3d )     , Intent( InOut ) :: FD
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs

    Real( wp ), Dimension( : ), Allocatable :: w2

    Integer :: this
    Integer :: i, j

    Call FD%set_dir_vecs( vecs )

    this = 0
    Do i = 1, 3
       Do j = i, 3
          this = this + 1
          FD%deriv_weights( this ) = Dot_Product( FD%get_inv_vec( i ), FD%get_inv_vec( j ) )
          If( i /= j ) FD%deriv_weights( this ) = 2.0_wp * FD%deriv_weights( this )
       End Do
    End Do

    w2 = FD%get_weight( 2 )
    ! Note allocation on assignment means indexing starts at 1
    FD%diag_inv = 1.0_wp / ( w2( 1 ) * ( FD%deriv_weights( XX ) + FD%deriv_weights( YY ) + FD%deriv_weights( ZZ ) ) )
    ! Work out if we need the off diagonal derivatives
    FD%need_XY = Abs( FD%deriv_weights( XY ) ) > FD%orthog_tol
    FD%need_XZ = Abs( FD%deriv_weights( XZ ) ) > FD%orthog_tol
    FD%need_YZ = Abs( FD%deriv_weights( YZ ) ) > FD%orthog_tol
    
  End Subroutine reset_vecs

  Subroutine apply( FD, grid_lb, lap_lb, start, final, grid, laplacian )
    !!-----------------------------------------------------------
    !! Calculate the resulting derivative by applying sequentially
    !! to the precalculated cache-blocks
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------
    Class( FD_Laplacian_3d )                                      , Intent( In    ) :: FD
    Integer, Dimension( 1:3 )                                     , Intent( In    ) :: grid_lb !! lower bounds grid
    Integer, Dimension( 1:3 )                                     , Intent( In    ) :: lap_lb  !! lower bounds laplacian
    Integer, Dimension( 1:3 )                                     , Intent( In    ) :: start   !! start point for calculation
    Integer, Dimension( 1:3 )                                     , Intent( In    ) :: final   !! final point for calculation
    Real( wp ), Dimension( grid_lb(1):, grid_lb(2):, grid_lb(3): ), Intent( In    ) :: grid    !! Thing to be differentiated
    Real( wp ), Dimension(  lap_lb(1):,  lap_lb(2):,  lap_lb(3): ), Intent(   Out ) :: laplacian

    Real( wp ), Dimension( : ), Allocatable :: w1, w2

    Integer :: order

    Integer :: i_block_3, i_block_2, i_block_1

    order = FD%get_order()

    w1 = FD%get_weight( 1 )
    w2 = FD%get_weight( 2 )

    !$omp do collapse( 3 )
    Do i_block_3 = start(3), final(3), FD%nc_block_apply( 3 )
       Do i_block_2 = start(2), final(2), FD%nc_block_apply( 2 )
          Do i_block_1 = start(1), final(1), FD%nc_block_apply( 1 )
             Call apply_block( [ i_block_1, i_block_2, i_block_3 ], &
                  grid_lb, lap_lb, FD%nc_block_apply, final, &
                  order, w1, w2, FD%deriv_weights, FD%need_XY, FD%need_XZ, FD%need_YZ, &
                  grid, laplacian )
          End Do
       End Do
    End Do
    !$omp end do

  End Subroutine apply

  Subroutine apply_block( s, lg, ll, nb, f, order, w1, w2, deriv_weights, need_XY, need_XZ, need_YZ, &
       grid, laplacian )

    !!-----------------------------------------------------------
    !! Apply the FD laplacian operator to part of the grid. An effort has been made
    !! to make sure the required data is in cache
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------

    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: s             !! Where to start
    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: lg            !! Lower bound of grid array
    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: ll            !! Lower bound of laplacian array
    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: nb            !! Cache blocking factors
    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: f             !! Where to finish
    Integer                         ,                      Intent( In    ) :: order         !! Order of the FD approximation
    Real( wp ), Dimension( -order: ),                      Intent( In    ) :: w1            !! FD Weights for first derivs
    Real( wp ), Dimension( -order: ),                      Intent( In    ) :: w2            !! FD Weights for second derivs
    Real( wp ), Dimension( 1:6     ),                      Intent( In    ) :: deriv_weights !! See below
    Logical                                              , Intent( In    ) :: need_XY       !! Need the XY deriv
    Logical                                              , Intent( In    ) :: need_XZ       !! Need the XZ deriv
    Logical                                              , Intent( In    ) :: need_YZ       !! Need the YZ deriv
    Real( wp ), Dimension( lg( 1 ):, lg( 2 ):, lg( 3 ): ), Intent( In    ) :: grid          !! The source
    Real( wp ), Dimension( ll( 1 ):, ll( 2 ):, ll( 3 ): ), Intent( InOut ) :: laplacian     !! The result

    ! Deriv_weights: As we do NOT assume the grid is orthogonal our FD laplacian is of the form
    ! d_xx * del_xx + d_xy * del_xy + d_xz * del_xz + d_yy * del_yy + d_yz * del_yz + d_zz * del_zz
    ! as we must finite difference along the directions of the grid. deriv_weights are the d_xx, d_xy
    ! etc. coefficients in this expression

    Real( wp ) :: st1, st2, st3
    Real( wp ) :: fac1, fac2, fac3
    Real( wp ) :: st12, st13, st23
    Real( wp ) :: fac12, fac13, fac23

    Integer :: i3, i2, i1
    Integer :: it, it1, it2, it3

    ! Order the loops for the various terms so that the inner loop is stride 1

    ! First do the xx, yy, zz terms. Assume the weights of these are always non-zero
    Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3  ) )
       Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
          Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )

             ! FD x^2, y^2, z^2 at grid point
             laplacian( i1, i2, i3 ) =                           w2( 0 ) * deriv_weights( XX ) * grid( i1, i2, i3 )
             laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + w2( 0 ) * deriv_weights( YY ) * grid( i1, i2, i3 )
             laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + w2( 0 ) * deriv_weights( ZZ ) * grid( i1, i2, i3 )

             ! FD x^2
             Do it = 1, order

                fac1 = w2( it ) * deriv_weights( XX )
                st1 = ( grid( i1 + it, i2     , i3      ) + grid( i1 - it, i2     , i3       ) )
                laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac1 * st1

             End Do

          End Do
       End Do
    End Do

    ! FD y^2, z^2
    Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
       Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
          Do it = 1, order
             Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )

                fac2 = w2( it ) * deriv_weights( YY )
                st2 = ( grid( i1     , i2 + it, i3      ) + grid( i1     , i2 - it, i3       ) )
                laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac2 * st2

                fac3 = w2( it ) * deriv_weights( ZZ )
                st3 = ( grid( i1     , i2     , i3 + it ) + grid( i1     , i2     , i3 - it  ) )
                laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac3 * st3

             End Do
          End Do
       End Do
    End Do

    ! xy
    If( need_XY ) Then
       Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
          Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
             Do it2 = 1, order
                Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
                   ! first derivs have zero weight at the grid point for central differences
                   Do it1 = 1, order
                      fac12 = w1( it1 ) * w1( it2 ) * deriv_weights( XY )
                      st12 = grid( i1 + it1, i2 + it2, i3 ) - grid( i1 - it1, i2 + it2, i3 ) - &
                           ( grid( i1 + it1, i2 - it2, i3 ) - grid( i1 - it1, i2 - it2, i3 ) )
                      laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac12 * st12
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If

    ! xz
    If( need_XZ ) Then
       Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
          Do it3 = 1, order
             Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
                Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
                   ! first derivs have zero weight at the grid point for central differences
                   Do it1 = 1, order
                      fac13 = w1( it1 ) * w1( it3 ) * deriv_weights( XZ )
                      st13 = grid( i1 + it1, i2, i3 + it3 ) - grid( i1 - it1, i2, i3 + it3 ) - &
                           ( grid( i1 + it1, i2, i3 - it3 ) - grid( i1 - it1, i2, i3 - it3 ) )
                      laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac13 * st13
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If

    ! yz
    If( need_YZ ) Then
       Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
          Do it3 = 1, order
             Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
                ! first derivs have zero weight at the grid point for central differences
                Do it2 = 1, order
                   Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
                      fac23 = w1( it2 ) * w1( it3 ) * deriv_weights( YZ )
                      st23 = grid( i1, i2 + it2, i3 + it3 ) - grid( i1, i2 - it2, i3 + it3 )  - &
                           ( grid( i1, i2 + it2, i3 - it3 ) - grid( i1, i2 - it2, i3 - it3 ) )
                      laplacian( i1, i2, i3 ) = laplacian( i1, i2, i3 ) + fac23 * st23
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If

  End Subroutine apply_block

  Subroutine jacobi_sweep( FD, grid_lb, lap_lb, start, final, jac_weight, grid, soln_in, soln_out )
    !!-----------------------------------------------------------
    !! Use the finite difference approximation of the laplacian operator
    !! to perform a sweep of the (weighted) Jacobi method
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------
    Class( FD_Laplacian_3d )                                      , Intent( In    ) :: FD
    Integer, Dimension( 1:3 )                                     , Intent( In    ) :: grid_lb    !! lower bounds grid
    Integer, Dimension( 1:3 )                                     , Intent( In    ) :: lap_lb     !! lower bounds laplacian
    Integer, Dimension( 1:3 )                                     , Intent( In    ) :: start      !! start point for calculation
    Integer, Dimension( 1:3 )                                     , Intent( In    ) :: final      !! final point for calculation
    Real( wp )                                                    , Intent( In    ) :: jac_weight !! The weight for the weighted jacobi
    Real( wp ), Dimension( grid_lb(1):, grid_lb(2):, grid_lb(3): ), Intent( In    ) :: grid       !! Thing to be differentiated
    Real( wp ), Dimension(  lap_lb(1):,  lap_lb(2):,  lap_lb(3): ), Intent( In    ) :: soln_in    !! Solution at start of sweep 
    Real( wp ), Dimension(  lap_lb(1):,  lap_lb(2):,  lap_lb(3): ), Intent(   Out ) :: soln_out   !! Solution at end of sweep

    Real( wp ), Dimension( : ), Allocatable :: w1, w2

    Integer :: order

    Integer :: i_block_3, i_block_2, i_block_1

    order = FD%get_order()

    w1 = FD%get_weight( 1 )
    w2 = FD%get_weight( 2 )

    !$omp do collapse( 3 )
    Do i_block_3 = start(3), final(3), FD%nc_block_jacobi( 3 )
       Do i_block_2 = start(2), final(2), FD%nc_block_jacobi( 2 )
          Do i_block_1 = start(1), final(1), FD%nc_block_jacobi( 1 )
             Call jacobi_sweep_block( [ i_block_1, i_block_2, i_block_3 ], &
                  grid_lb, lap_lb, FD%nc_block_apply, final, &
                  order, FD%diag_inv, w1, w2, FD%deriv_weights, &
                  FD%need_XY, FD%need_XZ, FD%need_YZ, &
                  jac_weight, grid, &
                  soln_in, soln_out )
          End Do
       End Do
    End Do
    !$omp end do

  End Subroutine jacobi_sweep

  Subroutine jacobi_sweep_block( s, lg, ll, nb, f, order, diag_inv, w1, w2, deriv_weights, &
       need_XY, need_XZ, need_YZ, jac_weight, grid, soln_in, soln_out )

    !!-----------------------------------------------------------
    !! Apply the FD laplacian operator to perform a Jacobi sweep on part of the grid
    !! An effort has been made to make sure all the data fits in cache
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------

    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: s             !! Where to start
    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: lg            !! Lower bound of grid array
    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: ll            !! Lower bound of laplacian array
    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: nb            !! Cache blocking factors
    Integer, Dimension( 1:3 )       ,                      Intent( In    ) :: f             !! Where to finish
    Integer                         ,                      Intent( In    ) :: order         !! Order of the FD approximation
    Real( wp )                      ,                      Intent( In    ) :: diag_inv      !! Inverse of the diagonal element of the matrix
    Real( wp ), Dimension( -order: ),                      Intent( In    ) :: w1            !! FD Weights for first derivs
    Real( wp ), Dimension( -order: ),                      Intent( In    ) :: w2            !! FD Weights for second derivs
    Real( wp ), Dimension( 1:6     ),                      Intent( In    ) :: deriv_weights !! See below
    Logical                                              , Intent( In    ) :: need_XY       !! Need the XY deriv
    Logical                                              , Intent( In    ) :: need_XZ       !! Need the XZ deriv
    Logical                                              , Intent( In    ) :: need_YZ       !! Need the YZ deriv
    Real( wp )                                           , Intent( In    ) :: jac_weight    !! The weight of the NEW solution 
    Real( wp ), Dimension( lg( 1 ):, lg( 2 ):, lg( 3 ): ), Intent( In    ) :: grid          !! The source
    Real( wp ), Dimension( ll( 1 ):, ll( 2 ):, ll( 3 ): ), Intent( In    ) :: soln_in       !! The solution on input
    Real( wp ), Dimension( ll( 1 ):, ll( 2 ):, ll( 3 ): ), Intent( InOut ) :: soln_out      !! The updated solution

    ! Deriv_weights: As we do NOT assume the grid is orthogonal our FD laplacian is of the form
    ! d_xx * del_xx + d_xy * del_xy + d_xz * del_xz + d_yy * del_yy + d_yz * del_yz + d_zz * del_zz
    ! as we must finite difference along the directions of the grid. deriv_weights are the d_xx, d_xy
    ! etc. coefficients in this expression

    Real( wp ) :: new_soln
    Real( wp ) :: st1, st2, st3
    Real( wp ) :: fac1, fac2, fac3
    Real( wp ) :: st12, st13, st23
    Real( wp ) :: fac12, fac13, fac23

    Integer :: i3, i2, i1
    Integer :: it, it1, it2, it3

    ! Order the loops for the various terms so that the inner loop is stride 1

    ! First do the xx, yy, zz terms. Assume the weights of these are always non-zero
    Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3  ) )
       Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
          Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )

             new_soln = grid( i1, i2, i3 )

             ! FD x^2
             Do it = 1, order

                fac1 = w2( it ) * deriv_weights( XX )
                st1 = ( soln_in( i1 + it, i2     , i3      ) + soln_in( i1 - it, i2     , i3       ) )
                new_soln = new_soln - fac1 * st1

             End Do

             soln_out( i1, i2, i3 ) = new_soln

          End Do
       End Do
    End Do

    ! FD y^2, z^2
    Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
       Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
          Do it = 1, order
             Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )

                fac2 = w2( it ) * deriv_weights( YY )
                st2 = ( soln_in( i1     , i2 + it, i3      ) + soln_in( i1     , i2 - it, i3       ) )
                soln_out( i1, i2, i3 ) = soln_out( i1, i2, i3 ) - fac2 * st2

                fac3 = w2( it ) * deriv_weights( ZZ )
                st3 = ( soln_in( i1     , i2     , i3 + it ) + soln_in( i1     , i2     , i3 - it  ) )
                soln_out( i1, i2, i3 ) = soln_out( i1, i2, i3 ) - fac3 * st3

             End Do
          End Do
       End Do
    End Do

    ! xy
    If( need_XY ) Then
       Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
          Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
             ! first derivs have zero weight at the diag point for central differences
             Do it2 = 1, order
                Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
                   Do it1 = 1, order
                      fac12 = w1( it1 ) * w1( it2 ) * deriv_weights( XY )
                      st12 = soln_in( i1 + it1, i2 + it2, i3 ) - soln_in( i1 - it1, i2 + it2, i3 ) - &
                           ( soln_in( i1 + it1, i2 - it2, i3 ) - soln_in( i1 - it1, i2 - it2, i3 ) )
                      soln_out( i1, i2, i3 ) = soln_out( i1, i2, i3 ) - fac12 * st12
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If

    ! xz
    If( need_XZ ) Then
       Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
          Do it3 = 1, order
             Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
                Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
                   ! first derivs have zero weight at the diag point for central differences
                   Do it1 = 1, order
                      fac13 = w1( it1 ) * w1( it3 ) * deriv_weights( XZ )
                      st13 = soln_in( i1 + it1, i2, i3 + it3 ) - soln_in( i1 - it1, i2, i3 + it3 ) - &
                           ( soln_in( i1 + it1, i2, i3 - it3 ) - soln_in( i1 - it1, i2, i3 - it3 ) )
                      soln_out( i1, i2, i3 ) = soln_out( i1, i2, i3 ) - fac13 * st13
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If

    ! yz
    If( need_YZ ) Then
       Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3 ) )
          Do it3 = 1, order
             Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
                ! first derivs have zero weight at the diag point for central differences
                Do it2 = 1, order
                   Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
                      fac23 = w1( it2 ) * w1( it3 ) * deriv_weights( YZ )
                      st23 = soln_in( i1, i2 + it2, i3 + it3 ) - soln_in( i1, i2 - it2, i3 + it3 )  - &
                           ( soln_in( i1, i2 + it2, i3 - it3 ) - soln_in( i1, i2 - it2, i3 - it3 ) )
                      soln_out( i1, i2, i3 ) = soln_out( i1, i2, i3 ) - fac23 * st23
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If

    ! Scale by diag element inverse and weight as required
    Do i3 = s( 3 ), Min( s( 3 ) + nb( 3 ) - 1, f( 3  ) )
       Do i2 = s( 2 ), Min( s( 2 ) + nb( 2 ) - 1, f( 2 ) )
          Do i1 = s( 1 ), Min( s( 1 ) + nb( 1 ) - 1, f( 1 ) )
             soln_out( i1, i2, i3 ) = soln_out( i1, i2, i3 ) * diag_inv
             soln_out( i1, i2, i3 ) = jac_weight * soln_out( i1, i2, i3 ) + ( 1.0_wp - jac_weight ) * soln_in( i1, i2, i3 )
          End Do
       End Do
    End Do

  End Subroutine jacobi_sweep_block

  Pure Function usage( n_grids, n_fd, nc_block, accuracy ) Result( reals )
    !!-----------------------------------------------------------
    !! Estimate memory usage of an algorithm that needs N_GRIDS grids
    !! without halos and N_FD grids with halos
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------
    Integer :: reals

    Integer,                   Intent( In ) :: n_grids
    Integer,                   Intent( In ) :: n_fd
    Integer, Dimension( 1:3 ), Intent( In ) :: nc_block
    Integer,                   Intent( In ) :: accuracy

    Integer :: grid
    Integer :: fd

    grid = n_grids * Product( nc_block )
    fd   = n_fd    * Product( nc_block + 1 + accuracy )

    reals = grid + fd

  End Function usage

End Module FD_Laplacian_3d_module
