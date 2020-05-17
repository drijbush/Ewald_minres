Module charge_grid_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Use lattice_module, Only : lattice

  Implicit None

  Public :: charge_grid_calculate, charge_grid_find_range, charge_grid_forces

  Private

  Real( wp ), Parameter, Private :: pi = 3.141592653589793238462643383279502884197_wp

Contains

  Subroutine charge_grid_calculate( l, alpha, q, r, range_gauss, lb, ub, q_grid, error )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    !$ Use omp_lib

    Use lattice_module, Only : lattice
    Use quadrature_base_module , Only : quadrature_base_class
    Use quadrature_trapezium_serial_module, Only : quadrature_trapezium_serial

    Implicit none

    Type( lattice )                    , Intent( In    ) :: l
    Real( wp ),                          Intent( In    ) :: alpha
    Real( wp ), Dimension( 1: ),         Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1: ),     Intent( In    ) :: r
    Integer   , Dimension( 1:3        ), Intent( In    ) :: range_gauss
    Integer   , Dimension( 1:3        ), Intent( In    ) :: lb( 1:3 )
    Integer   , Dimension( 1:3        ), Intent( In    ) :: ub( 1:3 )
    Real( wp ), Dimension( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) ), Intent(   Out ) :: q_grid
    Integer                            , Intent(   Out ) :: error

    Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp

    ! If set adjust the charge so it adds to zero. If we have zero charge there is
    ! gauranteed to be a solution of the linear equations, as the RHS vector
    ! has no component of the null space of the LHS matrix
    ! If we have to do this this really tells us that our representation of the charge density
    ! is not good enough, so print a warning
    Logical   , Parameter :: stabilise_q     = .True.
!!$    Logical   , Parameter :: stabilise_q     = .False.

    Real( wp ), Parameter :: stabilise_q_tol = 1e-10_wp

    Class( quadrature_base_class   ), Allocatable :: grid
    
    Real( wp ), Dimension( :, :, :, : ), Allocatable :: q_grid_red_hack

    Real( wp ), Dimension( 1:3 ) :: ri
    Real( wp ), Dimension( 1:3 ) :: fi
    Real( wp ), Dimension( 1:3 ) :: f_point
    Real( wp ), Dimension( 1:3 ) :: r_point
    Real( wp ), Dimension( 1:3 ) :: grid_vec

    Real( wp ) :: q_norm
    Real( wp ) :: qi_norm
    Real( wp ) :: q_val
    Real( wp ) :: q_tot, q_av

    Integer, Dimension( 1:3 ) :: domain_lo, n_domain, domain_hi
    
    Integer, Dimension( 1:3 ) :: n_grid
    Integer, Dimension( 1:3 ) :: i_atom_centre
    Integer, Dimension( 1:3 ) :: i_point
    Integer, Dimension( 1:3 ) :: i_grid
    Integer, Dimension( 1:3 ) :: i_g_lo, i_g_hi

    Integer :: n
    Integer :: n_th, iam
    Integer :: i1, i2, i3
    Integer :: i, i_th

    error = 0

    n = Size( q )
    n_grid = Ubound( q_grid ) + 1

    !HACK - test in serial before parallel
    domain_lo = [ 0, 0, 0 ]
    n_domain = n_grid
    
    domain_hi = domain_lo + n_domain - 1
    
    q_norm = ( ( ( alpha * alpha ) / pi ) ** 1.5_wp )

    ! What follows is a complete HACK to avoid gfortran stupidly putting
    ! arrays it is reducing on the stack and hence seg faulting once the
    ! array gets beyond a certain size - beyond belief, going to whinge about this one!
    !
    ! The proper fix is to //ise the loops over the grid points assoicated with a given
    ! atom. This is a bit complicated currently due to the Modulo meaning if we have
    ! a wide gaussian we can get a race condition. But we have to fix exactly this sort of
    ! thing when going to MPI, hence fix this properly then. At the moment
    ! just want correctness, hence simple HACK as grids really aren't that large

    n_th = 1
    ! Yes, omp_max_threads is a HACK as well ....
    !$ n_th = omp_get_max_threads()
    Allocate( q_grid_red_hack( Lbound( q_grid, Dim = 1 ):Ubound( q_grid, Dim = 1 ), &
         Lbound( q_grid, Dim = 2 ):Ubound( q_grid, Dim = 2 ), &
         Lbound( q_grid, Dim = 3 ):Ubound( q_grid, Dim = 3 ), 0:n_th - 1 ) )

    ! Should have default none - unfortunately bug in gfortran 7 stops this working
    !$omp parallel shared( n, l, alpha, r, q, q_norm, n_grid, q_grid, range_gauss, q_grid_red_hack, n_th, &
    !$omp                                  domain_lo, domain_hi ) &
    !$omp                          private( i, i1, i2, i3, qi_norm, ri, fi, i_atom_centre,   &
    !$omp                                   i_point, f_point, r_point, grid_vec, q_val, i_grid, iam, i_th,   &
    !$omp                                   i_g_lo, i_g_hi )
    !$omp do collapse( 3 )
    Do i3 = Lbound( q_grid, Dim = 3 ), Ubound( q_grid, Dim = 3 )
       Do i2 = Lbound( q_grid, Dim = 2 ), Ubound( q_grid, Dim = 2 )
          Do i1 = Lbound( q_grid, Dim = 1 ), Ubound( q_grid, Dim = 1 )
             q_grid( i1, i2, i3 ) = 0.0_wp
          End Do
       End Do
    End Do
    iam = 0
    !$ iam = omp_get_thread_num()
    Do i3 = Lbound( q_grid, Dim = 3 ), Ubound( q_grid, Dim = 3 )
       Do i2 = Lbound( q_grid, Dim = 2 ), Ubound( q_grid, Dim = 2 )
          Do i1 = Lbound( q_grid, Dim = 1 ), Ubound( q_grid, Dim = 1 )
             q_grid_red_hack( i1, i2, i3, iam ) = 0.0_wp
          End Do
       End Do
    End Do
    ! No need to sync at this point as in the next loop nest we only write to the
    ! bit we have just initalised
    ! Loop over atoms
    !$omp do 
    Do i = 1, n
       ! Loop over points associated with this atom which are in this domain
       ! Find point nearest to the atom, and call this the centre for the atom grid
       ! Assumes atom in fractional 0 < ri < 1
       ri = r( :, i )
       qi_norm = q_norm * q( i )
       Call l%to_fractional( ri, fi )
       i_atom_centre = Nint( fi * n_grid )
       i_g_lo = Max( i_atom_centre - range_gauss, domain_lo )
       i_g_hi = Min( i_atom_centre + range_gauss, domain_hi )
       Do i3 = i_g_lo( 3 ), i_g_hi( 3 )
          Do i2 = i_g_lo( 2 ), i_g_hi( 2 )
             Do i1 = i_g_lo( 1 ), i_g_hi( 1 )
                ! The indices of the point in space
                i_point = [ i1, i2, i3 ]
                ! Transform to fractional coordinates
                f_point = Real( i_point, wp ) / n_grid
                ! And fractional to real
                Call l%to_direct( f_point, r_point )
                ! Calculate the contribution to the total charge at the point
                ! R_POINT due to the charge distribution I
                grid_vec = r_point - ri
                q_val = qi_norm * Exp( - alpha * alpha * Dot_product( grid_vec, grid_vec ) )
                ! And add in
                i_grid = i_point
                q_grid_red_hack( i_grid( 1 ), i_grid( 2 ), i_grid( 3 ), iam ) = &
                     q_grid_red_hack( i_grid( 1 ), i_grid( 2 ), i_grid( 3 ), iam ) + q_val
             End Do
          End Do
       End Do
    End Do
    !$omp end do
    ! Synced here so OK to add up
    ! Now do the reduction manually - not needed once hacky way is fixed
    Do i_th = 0, n_th - 1
       !$omp do collapse( 3 )
       Do i3 = Lbound( q_grid, Dim = 3 ), Ubound( q_grid, Dim = 3 )
          Do i2 = Lbound( q_grid, Dim = 2 ), Ubound( q_grid, Dim = 2 )
             Do i1 = Lbound( q_grid, Dim = 1 ), Ubound( q_grid, Dim = 1 )
                q_grid( i1, i2, i3 ) = q_grid( i1, i2, i3 ) + q_grid_red_hack( i1, i2, i3, i_th )
             End Do
          End Do
       End Do
       !$omp end do
    End Do
    !$omp end parallel

    ! If required carefully make sure the charge on the grid adds to zero
    If( stabilise_q ) Then
       Allocate( quadrature_trapezium_serial :: grid )
       q_tot = grid%integrate( l, n_grid, q_grid ) * Product( n_grid ) / l%get_volume()
       If( Abs( q_tot ) > stabilise_q_tol ) Then
          error = -1
       End If
       q_av = q_tot / Product( n_grid )
       q_grid = q_grid - q_av
       q_tot = grid%integrate( l, n_grid, q_grid ) * Product( n_grid ) / l%get_volume()
       If( Abs( q_tot ) > stabilise_q_tol ) Then
          error = -2
       End If
    End If

  End Subroutine charge_grid_calculate

  Subroutine charge_grid_forces( l, alpha, q, r, range_gauss, pot_swapper, lb, ub, pot_grid, ei, f )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use lattice_module          , Only : lattice
    Use halo_setter_base_module , Only : halo_setter_base_class

    Implicit none
    
    Type( lattice )                    , Intent( In    ) :: l
    Real( wp ),                          Intent( In    ) :: alpha
    Real( wp ), Dimension( 1:         ), Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1:     ), Intent( In    ) :: r
    Integer   , Dimension( 1:3        ), Intent( In    ) :: range_gauss
    Class( halo_setter_base_class )    , Intent( InOut ) :: pot_swapper
    Integer   , Dimension( 1:3        ), Intent( In    ) :: lb( 1:3 )
    Integer   , Dimension( 1:3        ), Intent( In    ) :: ub( 1:3 )
    Real( wp ), Dimension( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) ), Intent( In    ) :: pot_grid
    Real( wp ), Dimension( 1:         ), Intent(   Out ) :: ei
    Real( wp ), Dimension( 1:, 1:     ), Intent(   Out ) :: f

    Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp

    Real( wp ), Dimension( :, :, : ), Allocatable :: pot_with_halo

    Real( wp ), Dimension( 1:3, 1:3 ) :: stress
    
    Real( wp ), Dimension( 1:3 ) :: ri
    Real( wp ), Dimension( 1:3 ) :: fi
    Real( wp ), Dimension( 1:3 ) :: f_point
    Real( wp ), Dimension( 1:3 ) :: r_point
    Real( wp ), Dimension( 1:3 ) :: grid_vec

    Real( wp ) :: q_norm
    Real( wp ) :: qi_norm
    Real( wp ) :: g_val
    Real( wp ) :: dV
    Real( wp ) :: s
    
    Integer, Dimension( 1:3 ) :: n_grid
    Integer, Dimension( 1:3 ) :: i_atom_centre
    Integer, Dimension( 1:3 ) :: i_atom_grid
    Integer, Dimension( 1:3 ) :: i_point
    Integer, Dimension( 1:3 ) :: i_grid

    Integer :: n
    Integer :: i1, i2, i3
    Integer :: i, i_alpha, i_beta
    Integer :: error

    n      = Size( q )
    n_grid = Ubound( pot_grid ) + 1
    dV     = l%get_volume() / Product( n_grid )

    q_norm = ( ( ( alpha * alpha ) / pi ) ** 1.5_wp )

    stress = 0.0_wp

    ! Set up the pot_grid with a halo
    ! PROBLEM - Should range_gauss should be an array? do we need different ranges in each direction?
    ! PROBLEM - Why +1? I think it may be because while an atom is in the current domain the nearest
    !           grid point could belong to the neighbouring domain if the atom is very close to
    !           the upper edge of the domain. But need to think this through. An alternative
    !           would be to use Floor instead of Nint in finding the centre grid point below,
    !           but this would be less accurate
    Call pot_swapper%allocate( lb, ub, range_gauss( 1 ) + 1, pot_with_halo )
    Call pot_swapper%fill( range_gauss( 1 ) + 1, Lbound( pot_with_halo ), pot_grid, pot_with_halo, error )
    If( error /= 0 ) Then
       Error Stop "halo filler problem in forces"
    End If
    
    !$omp parallel shared( n, l, alpha, r, q, q_norm, n_grid, pot_grid, range_gauss, dV, ei, f, stress ) &
    !$omp                          private( i, i1, i2, i3, qi_norm, ri, fi, i_atom_centre,   &
    !$omp                                   i_atom_grid, i_point, f_point, r_point, grid_vec, g_val, i_grid, s )
    ! Loop over atoms
    ! Stress small array, size invariant, so reduction won't stack smash
    !$omp do reduction( +:stress )
    Do i = 1, n
       ! Loop over points associated with atoms
       ! Find point nearest to the atom, and call this the centre for the atom grid
       ! Assumes atom in fractional 0 < ri < 1
       ri = r( :, i )
       qi_norm = q_norm * q( i ) * dV
       Call l%to_fractional( ri, fi )
       i_atom_centre = Nint( fi * n_grid )
       !
       ei(    i ) = 0.0_wp
       f ( :, i ) = 0.0_wp
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
                ! Vector to the point of interest from the centre of the gaussin
                grid_vec = r_point - ri
                ! Gaussian at that point times normalisation times the volume element
                g_val = qi_norm * Exp( - alpha * alpha * Dot_product( grid_vec, grid_vec ) )
                i_grid = i_point
                ! Include the potential term
                g_val = g_val * pot_with_halo( i_grid( 1 ), i_grid( 2 ), i_grid( 3 ) )
                ! Add into the per particle energy and the force
                ei(    i ) = ei(    i ) + 0.5_wp *                            g_val
                f ( :, i ) = f ( :, i ) - 2.0_wp * alpha * alpha * grid_vec * g_val
                ! Stress term - TOTALLY UNTESTED AND PROBABLY INCOMPLETE
                ! Need short range term due to differentiation of coulomb operator,
                ! and SIC term
                Do i_beta = 1, 3
                   Do i_alpha = i_beta, 3
                      s = - 2.0_wp * g_val * alpha * alpha * grid_vec( i_alpha ) * grid_vec( i_beta )
                      stress( i_alpha, i_beta ) = stress( i_alpha, i_beta ) + s 
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do
    !$omp end do
    !$omp end parallel

    ! Stress - Untested!
    ! Symmetrize stress tensor
    Do i_beta = 1, 2
       Do i_alpha = i_beta + 1, 3
          stress( i_beta, i_alpha ) = stress( i_alpha, i_beta )
       End Do
    End Do
    
  End Subroutine charge_grid_forces

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
