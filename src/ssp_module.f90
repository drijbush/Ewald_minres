Module symetrically_screened_poisson_module

  Use constants, Only : wp, li, pi
  Implicit None

  Public :: ssp_long_range
  Public :: ssp_sic

  Private

Contains

  Subroutine ssp_long_range( l, q, r, alpha, q_halo, r_halo, range_gauss, n_grid, lb, rtol, &
       recip_E, q_grid, pot_grid, solver, pot_swapper, grid_integrator, &
       ei, f, t_grid, t_pot_solve, t_forces, itn, istop, istop_message, rnorm, error, q_grid_old, pot_grid_old )

    Use lattice_module,                           Only : lattice
    Use charge_grid_module,                       Only : charge_grid_calculate, charge_grid_forces
    Use halo_setter_base_module,                  Only : halo_setter_base_class
    Use quadrature_base_module,                   Only : quadrature_base_class
    Use equation_solver_precon_base_class_module, Only : equation_solver_precon_base_class

    Implicit None

    Real( wp ), Dimension( :, :, : ), Allocatable :: rhs

    Type( lattice ),                     Intent( In    ) :: l
    Real( wp ), Dimension( 1:     ),     Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1: ),     Intent( In    ) :: r
    Real( wp ),                          Intent( In    ) :: alpha
    Real( wp ), Dimension( 1:     ),     Intent( In    ) :: q_halo
    Real( wp ), Dimension( 1:, 1: ),     Intent( In    ) :: r_halo
    Integer,                             Intent( In    ) :: range_gauss
    Integer,    Dimension( 1:3 ),        Intent( In    ) :: n_grid
    Integer,    Dimension( 1:3 ),        Intent( In    ) :: lb
    Real( wp ),                          Intent(   Out ) :: recip_E
    Real( wp ),                          Intent( In    ) :: rtol
    Real( wp ), Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: q_grid
    Real( wp ), Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: pot_grid
    Class( equation_solver_precon_base_class ), Intent( InOut ) :: solver ! Need to check and fix these InOuts
    Class( halo_setter_base_class  ),    Intent( InOut ) :: pot_swapper
    Class( quadrature_base_class   ),    Intent( InOut ) :: grid_integrator
    Real( wp ), Dimension( 1: ),         Intent(   Out ) :: ei
    Real( wp ), Dimension( 1:, 1: ),     Intent(   Out ) :: f
    Real( wp ),                          Intent(   Out ) :: t_grid
    Real( wp ),                          Intent(   Out ) :: t_pot_solve
    Real( wp ),                          Intent(   Out ) :: t_forces
    Integer,                             Intent(   Out ) :: itn
    Integer,                             Intent(   Out ) :: istop
    Character( Len = * ),                Intent(   Out ) :: istop_message
    Real( wp ),                          Intent(   Out ) :: rnorm
    Integer,                             Intent(   Out ) :: error
    Real( wp ), Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ), Optional :: q_grid_old
    Real( wp ), Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ), Optional :: pot_grid_old

    ! Standardize the potential so it sums to zero over the cell
    ! Not required, but useful for comparison of accuracy with Fourier methods.
    ! This should be set to .False. for production
    Logical, Parameter :: standardise = .True.

    Real( wp ), Dimension( :, : ), Allocatable :: r_full

    Integer( li ) :: start, finish, rate

    Allocate( r_full( 1:3, 1:Size( q ) + Size( q_halo ) ) )
    r_full( :, 1:Size( q )    ) = r
    r_full( :, Size( q ) + 1: ) = r_halo

    error = 0

    ! Grid the charge
    Call System_clock( start, rate )
    Call charge_grid_calculate( l, alpha, [ q, q_halo ], r_full, range_gauss, &
         n_grid, Lbound( q_grid ), Ubound( q_grid ), solver%comms, grid_integrator, q_grid, error )
    Call System_clock( finish, rate )
    t_grid = Real( finish - start, wp ) / rate

    ! Now calculate the long range potential by Finite difference

    ! And solve  Possion equation on the grid by FDs
    Call System_clock( start, rate )
    ! Solve on delta if the old grid is present
    If( Present( q_grid_old ) .And. Present( pot_grid_old ) ) Then
       rhs = - 4.0_wp * pi * ( q_grid - q_grid_old )
    Else
       rhs = - 4.0_wp * pi * q_grid
    End If
    Call solver%solve( Lbound( q_grid ), Ubound( q_grid ), rhs, rtol, pot_grid, istop, istop_message, itn, rnorm )
    ! If delta solve add back in old potential
    If( Present( q_grid_old ) .And. Present( pot_grid_old ) ) Then
       pot_grid = pot_grid + pot_grid_old
    End If
    If( standardise ) Then
       ! Standardise to potential averages to zero over grid
       ! In real calculation don't need to do this!
       pot_grid = pot_grid - grid_integrator%integrate( solver%comms, l, n_grid, pot_grid ) / l%get_volume()
    End If
    Call System_clock( finish, rate )
    t_pot_solve = Real( finish - start, wp ) / rate

    ! Calculate from the potential the long range energy
    ! Two minus signs as a) this is the SCREENED charge
    !                    b) I like to just add up the energies, having to subtract it is confusing
    recip_E = - 0.5_wp * grid_integrator%integrate( solver%comms, l, n_grid, - q_grid * pot_grid )

    ! Calculate the forces and energy per site
    ! Initalise the halo swapper
    Call System_clock( start, rate )
    Call charge_grid_forces( l, alpha, q, r, range_gauss, n_grid, pot_swapper, Lbound( pot_grid ), Ubound( pot_grid ), &
         pot_grid, ei, f )
    Call System_clock( finish, rate )
    t_forces = Real( finish - start, wp ) / rate

  End Subroutine ssp_long_range

  Pure Function ssp_sic( q, alpha ) Result( sic )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Real( wp ), Dimension( : ), Intent( In    ) :: q
    Real( wp ),                 Intent( In    ) :: alpha

    Real( wp ) :: sic

    sic = - alpha * Sum( q * q ) / Sqrt( 2.0_wp * pi )

  End Function ssp_sic

  Subroutine dummy_msolve( lb, ub, x, y )                   ! Solve M*y = x

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Integer,                                                Intent( In    ) :: lb( 1:3 )
    Integer,                                                Intent( In    ) :: ub( 1:3 )
    Real( wp ),  Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent( In    ) :: x
    Real( wp ),  Dimension( lb( 1 ):, lb( 2 ):, lb( 3 ): ), Intent(   Out ) :: y

    ! Shut up compiler
    if (lb(1) /= 0) continue
    if (ub(1) /= 0) continue

    y = x

  End Subroutine dummy_msolve

End Module symetrically_screened_poisson_module
