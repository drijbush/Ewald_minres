Module symetrically_screened_poisson_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

  Implicit None

  Public :: ssp_long_range
  Public :: ssp_sic
  
  Private
  
  Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp

Contains

  Subroutine ssp_long_range( l, q, r, alpha, FD_order, q_halo, r_halo,      &
       recip_E, q_grid, pot_grid, comms, fd_swapper, pot_swapper, grid_integrator, &
       ei, f, t_grid, t_recip, error )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

    Use lattice_module         , Only : lattice
    Use charge_grid_module     , Only : charge_grid_calculate, charge_grid_find_range, charge_grid_forces
    Use minresmodule           , Only : minres
    Use FD_Laplacian_3d_module , Only : FD_Laplacian_3D
    Use comms_base_class_module, Only : comms_base_class
    Use halo_setter_base_module, Only : halo_setter_base_class
    Use quadrature_base_module , Only : quadrature_base_class
    
    Implicit None

    Real( wp ), Dimension( :, :, : ), Allocatable :: rhs
    
    Type( lattice )                    , Intent( In    ) :: l
    Real( wp ), Dimension( 1:     )    , Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1: )    , Intent( In    ) :: r
    Real( wp )                         , Intent( In    ) :: alpha
    Integer                            , Intent( In    ) :: FD_order
    Real( wp ), Dimension( 1:     )    , Intent( In    ) :: q_halo
    Real( wp ), Dimension( 1:, 1: )    , Intent( In    ) :: r_halo
    Real( wp )                         , Intent(   Out ) :: recip_E
    Real( wp ), Dimension( 0:, 0:, 0: ), Intent(   Out ) :: q_grid
    Real( wp ), Dimension( 0:, 0:, 0: ), Intent(   Out ) :: pot_grid
    Class( comms_base_class        )   , Intent( InOut ) :: comms
    Class( halo_setter_base_class  )   , Intent( InOut ) :: fd_swapper
    Class( halo_setter_base_class  )   , Intent( InOut ) :: pot_swapper
    Class( quadrature_base_class   )   , Intent( InOut ) :: grid_integrator
    Real( wp ), Dimension( 1: )        , Intent(   Out ) :: ei
    Real( wp ), Dimension( 1:, 1: )    , Intent(   Out ) :: f
    Real( wp )                         , Intent(   Out ) :: t_grid
    Real( wp )                         , Intent(   Out ) :: t_recip
    Integer                            , Intent(   Out ) :: error

    ! Standardize the potential so it sums to zero over the cell
    ! Not required, but useful for comparison of accuracy with Fourier methods.
    ! This should be set to .False. for production
    Logical, Parameter :: standardise = .True.

    Type( FD_Laplacian_3d    ) :: FD

    Real( wp ), Dimension( :, : ), Allocatable :: r_full
    
    Real( wp ), Dimension( 1:3, 1:3 ) :: dGrid_vecs

    Real( wp ), Dimension( 1:3 ) :: dG

    Integer, Dimension( 1:3 ) :: n_grid
    Integer, Dimension( 1:3 ) :: range_gauss

    Real( wp ) :: Anorm, Arnorm, Acond, rnorm, ynorm, rtol

    Integer :: istop, itn
    Integer :: i
    
    Integer( li ) :: start, finish, rate

    Character( Len = 132 ) :: istop_message

    Allocate( r_full( 1:3, 1:Size( q ) + Size( q_halo ) ) )
    r_full( :, 1:Size( q )    ) = r
    r_full( :, Size( q ) + 1: ) = r_halo

    error = 0

    n_grid = Ubound( q_grid ) + 1

    ! Grid the charge
    Call system_Clock( start, rate )
    ! First find range of the gaussian along each of the axes of the grid
    Call charge_grid_find_range( l, alpha, n_grid, range_gauss )
    ! Now grid the charge
    Call charge_grid_calculate( l, alpha, [ q, q_halo ], r_full, range_gauss, &
         Lbound( q_grid ), Ubound( q_grid ), comms, grid_integrator, q_grid, error )
    Call system_Clock( finish, rate )
    t_grid = Real( finish - start, wp ) / rate

    ! Now calculate the long range potential by Finite difference

    ! Initialise the FD template
    dGrid_vecs = l%get_direct_vectors()
    Do i = 1, 3
       dGrid_vecs( :, i ) = dGrid_vecs( :, i ) / n_grid( i )
       dG( i )            = Sqrt( Dot_Product( dGrid_vecs( :, i ), dGrid_vecs( :, i ) ) )
    End Do
    Call FD%init( FD_order, dGrid_vecs )

    ! And solve  Possion equation on the grid by FDs
    rtol = 1.0e-12_wp
    Call system_Clock( start, rate )
    rhs = - 4.0_wp * pi * q_grid
    Call minres( Lbound( q_grid ), Ubound( q_grid ), FD, comms, fd_swapper, dummy_Msolve, rhs, 0.0_wp, .True., .False., &
         pot_grid, 1000, 99, rtol,                      &
         istop, istop_message, itn, Anorm, Acond, rnorm, Arnorm, ynorm )
    If( standardise ) Then
       ! Standardise to potential averages to zero over grid
       ! In real calculation don't need to do this!
       pot_grid = pot_grid - grid_integrator%integrate( comms, l, n_grid, pot_grid ) / l%get_volume()
    End If
    Call system_Clock( finish, rate )
    t_recip = Real( finish - start, wp ) / rate

    ! Summarise the iterative solver
    Write( *, * ) 'Iterative solver summary:'
    Write( *, * ) 'alpha                = ', alpha
    Write( *, * ) 'Grid resolution      = ', dG
    Write( *, * ) 'Grid size            = ', n_grid
    Write( *, * ) 'Order                = ', FD_order
    Write( *, * ) 'iterations           = ', itn
    Write( *, * ) 'istop                = ', istop, Trim( istop_message )
    Write( *, * ) 'norm of the residual = ', rnorm

    ! Calculate from the potential the long range energy
    ! Two minus signs as a) this is the SCREENED charge
    !                    b) I like to just add up the energies, having to subtract it is confusing
    recip_E = - 0.5_wp * grid_integrator%integrate( comms, l, n_grid, - q_grid * pot_grid )
    
    ! Calculate the forces and energy per site
    ! Initalise the halo swapper
    Call pot_swapper%init( error )
    Call charge_grid_forces( l, alpha, q, r, range_gauss, pot_swapper, Lbound( pot_grid ), Ubound( pot_grid ), &
         pot_grid, ei, f )
    
  End Subroutine ssp_long_range

  Pure Function ssp_sic( q, alpha ) Result( sic )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Real( wp ), Dimension( : ), Intent( In    ) :: q
    Real( wp )                , Intent( In    ) :: alpha

    Real( wp ) :: sic

    sic = - alpha * Sum( q * q ) / Sqrt( 2.0_wp * pi )

  End Function ssp_sic

  Subroutine dummy_msolve( lb, ub, x, y )                   ! Solve M*y = x
    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
    Integer,    Intent( In    )   :: lb( 1:3 ), ub( 1:3 )
    Real( wp ), Intent( In    )   :: x( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) )
    Real( wp ), Intent(   out )   :: y( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) )
    
    ! Shut up compiler
    y = x

  End Subroutine dummy_msolve


End Module symetrically_screened_poisson_module
