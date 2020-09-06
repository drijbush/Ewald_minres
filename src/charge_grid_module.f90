Module charge_grid_module

  Use constants, Only : wp, pi
  Use lattice_module, Only : lattice

  Implicit None

  Public :: charge_grid_calculate, charge_grid_find_range, charge_grid_forces, charge_grid_get_n_grid

  Private

Contains

  Subroutine charge_grid_calculate( l, alpha, q, r, range_gauss, n_grid, lb, ub, comms, grid_integrator, q_grid,error )

    !$ Use omp_lib, only : omp_get_max_threads, omp_get_thread_num
    Use lattice_module,          Only : lattice
    Use comms_base_class_module, Only : comms_base_class
    Use quadrature_base_module,  Only : quadrature_base_class

    Implicit None

    Type( lattice ),                     Intent( In    ) :: l
    Real( wp ),                          Intent( In    ) :: alpha
    Real( wp ), Dimension( 1: ),         Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1: ),     Intent( In    ) :: r
    Integer,                             Intent( In    ) :: range_gauss
    Integer,    Dimension( 1:3        ), Intent( In    ) :: n_grid ! global size of grid
    Integer,    Dimension( 1:3        ), Intent( In    ) :: lb( 1:3 )
    Integer,    Dimension( 1:3        ), Intent( In    ) :: ub( 1:3 )
    Class( quadrature_base_class      ), Intent( In    ) :: grid_integrator
    Class( comms_base_class           ), Intent( In    ) :: comms
    Real( wp ), Dimension( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) ), Intent(   Out ) :: q_grid
    Integer,                             Intent(   Out ) :: error

    ! If set adjust the charge so it adds to zero. If we have zero charge there is
    ! gauranteed to be a solution of the linear equations, as the RHS vector
    ! has no component of the null space of the LHS matrix
    ! If we have to do this this really tells us that our representation of the charge density
    ! is not good enough, so print a warning
    Logical,    Parameter :: stabilise_q     = .True.

    Real( wp ), Parameter :: stabilise_q_tol = 1e-10_wp

    Real( wp ), Dimension( :, :, :, : ), Allocatable :: q_grid_red_hack

    Real( wp ), Dimension( 1:3 ) :: ri
    Real( wp ), Dimension( 1:3 ) :: fi
    Real( wp ), Dimension( 1:3 ) :: f_point
    Real( wp ), Dimension( 1:3 ) :: r_point
    Real( wp ), Dimension( 1:3 ) :: grid_vec

    Real( wp ), Dimension( 1:3 ) :: r_0
    Real( wp ), Dimension( 1:3, 1:3 ) :: dr
    Real( wp ) :: g_r
    Real( wp ), Dimension(3) :: g_r0, g_dr, f_rdr, f_rdr0, f_rdr00, f_drdr
    Real( wp ), Dimension(3,3) :: g_drdr

    Real( wp ) :: q_norm
    Real( wp ) :: qi_norm
    Real( wp ) :: q_val, p_val
    Real( wp ) :: q_tot, q_av

    Integer, Dimension( 1:3 ) :: domain_lo, n_domain, domain_hi

    Integer, Dimension( 1:3 ) :: i_atom_centre
    Integer, Dimension( 1:3 ) :: i_point
    Integer, Dimension( 1:3 ) :: i_grid
    Integer, Dimension( 1:3 ) :: i_g_lo, i_g_hi

    Integer :: n_th, iam
    Integer :: i1, i2, i3
    Integer :: i, i_th

    error = 0

    domain_lo = Lbound( q_grid )
    n_domain  = Ubound( q_grid ) - Lbound( q_grid ) + 1

    domain_hi = domain_lo + n_domain - 1

    q_norm = ( ( ( alpha * alpha ) / pi ) ** 1.5_wp )

    dr = l%get_direct_vectors()
    Do i1 = 1, 3
       dr( :, i1 ) = dr( :, i1 ) / n_grid( i1 )
       g_dr(i1) = g(dr(:, i1), alpha)
       f_drdr(i1) = f(dr(:, i1), dr(:, i1), alpha)
    end do

    do i2 = 1, 3
      do i1 = 1, 3
        g_drdr(i1, i2) = f(dr(:,i1), dr(:,i2), alpha)
      end do
    end do

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
    !$omp parallel shared( l, alpha, r, q, q_norm, n_grid, q_grid, range_gauss, q_grid_red_hack, n_th, &
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
    Do i = 1, Size( q )
       ! Loop over points associated with this atom which are in this domain
       ! Find point nearest to the atom, and call this the centre for the atom grid
       ! Assumes atom in fractional 0 < ri < 1
       ri = r( :, i )
       qi_norm = q_norm * q( i )
       Call l%to_fractional( ri, fi )
       i_atom_centre = Nint( fi * n_grid )
       i_g_lo = Max( i_atom_centre - range_gauss, domain_lo )
       i_g_hi = Min( i_atom_centre + range_gauss, domain_hi )

       ! Transform to fractional coordinates
       f_point = Real( i_atom_centre, wp ) / n_grid
       ! And fractional to real
       Call l%to_direct( f_point, r_point )
       ! Vector to the point of interest from the centre of the gaussian
       grid_vec = r_point - ri
       
       
!       r_0 = i_g_lo(1) * dr(:,1) + i_g_lo(2) * dr(:,2) + i_g_lo(3) * dr(:,3)
       r_0 = matmul(dr, i_atom_centre - i_g_lo)
       r_0 = grid_vec - r_0

       ! Gaussian at zero point
       g_r = g(r_0, alpha)
       g_r0 = g_r

       do i1 = 1,3
         f_rdr(i1) = f(r_0, dr(:, i1), alpha)
       end do
       f_rdr0 = f_rdr
       f_rdr00 = f_rdr

       r_0 = grid_vec

       Do i3 = i_g_lo( 3 ), i_g_hi( 3 )
          Do i2 = i_g_lo( 2 ), i_g_hi( 2 )
             Do i1 = i_g_lo( 1 ), i_g_hi( 1 )
                ! ! The indices of the point in space
                i_point = [ i1, i2, i3 ]
                ! Transform to fractional coordinates
                ! f_point = Real( i_point, wp ) / n_grid
                ! ! And fractional to real
                ! Call l%to_direct( f_point, r_point )
                ! ! Calculate the contribution to the total charge at the point
                ! ! R_POINT due to the charge distribution I
                ! grid_vec = r_point - ri
                ! q_val = qi_norm * Exp( - alpha * alpha * Dot_product( grid_vec, grid_vec ) )
                ! print*, "OLD", i_point, q_val

                ! And add in
                q_val = qi_norm * g_r
                ! print*, "NEW", i_point, q_val
                g_r = g_r * f_rdr(1) * g_dr(1)
                f_rdr = f_rdr * g_drdr(:, 1)

                i_grid = i_point

                q_grid_red_hack( i_grid( 1 ), i_grid( 2 ), i_grid( 3 ), iam ) = &
                     q_grid_red_hack( i_grid( 1 ), i_grid( 2 ), i_grid( 3 ), iam ) + q_val
              End Do

              g_r = g_r0(2) * f_rdr0(2) * g_dr(2)
              g_r0(2) = g_r

              f_rdr = f_rdr0 * g_drdr(:, 2)
              f_rdr0 = f_rdr

          End Do

          g_r = g_r0(1) * f_rdr00(3) * g_dr(3)
          g_r0 = g_r

          f_rdr = f_rdr00 * g_drdr(:, 3)
          f_rdr0 = f_rdr
          f_rdr00 = f_rdr

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
       q_tot = grid_integrator%integrate( comms, l, n_grid, q_grid ) * Product( n_grid ) / l%get_volume()
       If( Abs( q_tot ) > stabilise_q_tol ) Then
          error = -1
       End If
       q_av = q_tot / Product( n_grid )
       q_grid = q_grid - q_av
       q_tot = grid_integrator%integrate( comms, l, n_grid, q_grid ) * Product( n_grid ) / l%get_volume()
       If( Abs( q_tot ) > stabilise_q_tol ) Then
          error = -2
       End If
    End If

  End Subroutine charge_grid_calculate

  Subroutine charge_grid_forces( l, alpha, q, r, range_gauss, n_grid, pot_swapper, lb, ub, pot_grid, ei, force )

    Use lattice_module,           Only : lattice
    Use halo_setter_base_module,  Only : halo_setter_base_class

    Implicit none

    Type( lattice ),                     Intent( In    ) :: l
    Real( wp ),                          Intent( In    ) :: alpha
    Real( wp ), Dimension( 1:         ), Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1:     ), Intent( In    ) :: r
    Integer,                             Intent( In    ) :: range_gauss
    Integer,    Dimension( 1:3        ), Intent( In    ) :: n_grid ! global size of grid
    Class( halo_setter_base_class ),     Intent( InOut ) :: pot_swapper
    Integer,    Dimension( 1:3        ), Intent( In    ) :: lb( 1:3 )
    Integer,    Dimension( 1:3        ), Intent( In    ) :: ub( 1:3 )
    Real( wp ), Dimension( lb( 1 ):ub( 1 ), lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) ), Intent( In    ) :: pot_grid
    Real( wp ), Dimension( 1:         ), Intent(   Out ) :: ei
    Real( wp ), Dimension( 1:, 1:     ), Intent(   Out ) :: force

    Real( wp ), Dimension( :, :, : ), Allocatable :: pot_with_halo

    Real( wp ), Dimension( 1:3, 1:3 ) :: stress
    Real( wp ), Dimension( 1:3, 1:3 ) :: dr
    Real( wp ), Dimension( 1:3, 1:3 ) :: g_drdr

    Real( wp ), Dimension( 1:3 ) :: ri
    Real( wp ), Dimension( 1:3 ) :: fi
    Real( wp ), Dimension( 1:3 ) :: f_point
    Real( wp ), Dimension( 1:3 ) :: r_point
    Real( wp ), Dimension( 1:3 ) :: grid_vec

    Real( wp ), Dimension( 1:3 ) :: r_0
    Real( wp ), Dimension( 1:3 ) :: g_r0, g_dr, f_rdr, f_rdr0, f_rdr00, f_drdr

    Real( wp ) :: q_norm
    Real( wp ) :: g_r
    Real( wp ) :: qi_norm
    Real( wp ) :: fix, fiy, fiz
    Real( wp ) :: g_val
    Real( wp ) :: dV
    Real( wp ) :: s

    Integer, Dimension( 1:3 ) :: i_atom_centre
    Integer, Dimension( 1:3 ) :: i_grid

    Integer :: n
    Integer :: i1, i2, i3
    Integer :: i
    Integer :: error

    n      = Size( q )
    dV     = l%get_volume() / Product( n_grid )

    q_norm = ( ( ( alpha * alpha ) / pi ) ** 1.5_wp )

    dr = l%get_direct_vectors()
    Do i1 = 1, 3
       dr( :, i1 ) = dr( :, i1 ) / n_grid( i1 )
       g_dr( i1 ) = g( dr( :, i1 ), alpha )
       f_drdr( i1 ) = f( dr( :, i1 ), dr( :, i1 ), alpha )
    end do

    do i2 = 1, 3
       do i1 = 1, 3
          g_drdr( i1, i2 ) = f( dr( :, i1 ), dr( :, i2 ), alpha )
       end do
    end do

    stress = 0.0_wp

    ! Set up the pot_grid with a halo
    ! PROBLEM - Should range_gauss should be an array? do we need different ranges in each direction?
    ! PROBLEM - Why +1? I think it may be because while an atom is in the current domain the nearest
    !           grid point could belong to the neighbouring domain if the atom is very close to
    !           the upper edge of the domain. But need to think this through. An alternative
    !           would be to use Floor instead of Nint in finding the centre grid point below,
    !           but this would be less accurate
    Call pot_swapper%allocate( lb, ub, range_gauss + 1, pot_with_halo )
    Call pot_swapper%fill( range_gauss + 1, Lbound( pot_with_halo ), pot_grid, pot_with_halo, error )
    If( error /= 0 ) Then
       Error Stop "halo filler problem in forces"
    End If

    !$omp parallel shared( n, l, alpha, r, q, q_norm, n_grid, pot_grid, range_gauss, &
    !$omp                  dV, ei, force, stress ) &
    !$omp          private( i, i1, i2, i3, qi_norm, ri, fi, i_atom_centre,   &
    !$omp                   f_point, r_point, grid_vec, g_val, i_grid, s )
    ! Loop over atoms
    ! Stress small array, size invariant, so reduction won't stack smash
    !$omp do reduction( +:stress )
    particle_loop: Do i = 1, n
       ! Loop over points associated with atoms
       ! Find point nearest to the atom, and call this the centre for the atom grid
       ! Assumes atom in fractional 0 < ri < 1
       ri = r( :, i )
       qi_norm = q_norm * q( i ) * dV
       Call l%to_fractional( ri, fi )
       i_atom_centre = Nint( fi * n_grid )

       ! Transform to fractional coordinates
       f_point = Real( i_atom_centre, wp ) / n_grid
       ! And fractional to real
       Call l%to_direct( f_point, r_point )
       ! Vector to the point of interest from the centre of the gaussian
       grid_vec = r_point - ri

       r_0 = grid_vec - sum( range_gauss * dr, dim=2 )

       ! Gaussian at zero point
       g_r = g( r_0, alpha )
       g_r0 = g_r

       do i1 = 1,3
          f_rdr( i1 ) = f( r_0, dr( :, i1 ), alpha )
       end do
       f_rdr0 = f_rdr
       f_rdr00 = f_rdr

       r_0 = grid_vec

       !
       ei(    i ) = 0.0_wp
!!$       force ( :, i ) = 0.0_wp
       fix = 0.0_wp
       fiy = 0.0_wp
       fiz = 0.0_wp
       Do i3 = - range_gauss, range_gauss
          Do i2 = - range_gauss, range_gauss
             Do i1 = - range_gauss, range_gauss

                ! Gaussian at that point times normalisation times the volume element
!!$                g_val = qi_norm * g_r
                g_val = g_r

                g_r = g_r * f_rdr( 1 ) * g_dr( 1 )
                f_rdr = f_rdr * g_drdr( :, 1 )

                i_grid = i_atom_centre + [ i1, i2, i3 ]

                grid_vec = r_0 + i1 * dr( :, 1 ) + i2 * dr( :, 2 ) + i3 * dr( :, 3 )

                ! Include the potential term
                g_val = g_val * pot_with_halo( i_grid( 1 ), i_grid( 2 ), i_grid( 3 ) )
                ! Add into the per particle energy and the force
                ei    (    i ) = ei    (    i ) +            g_val
!!$                force ( :, i ) = force ( :, i ) - grid_vec * g_val
                fix = fix - grid_vec( 1 ) * g_val
                fiy = fiy - grid_vec( 2 ) * g_val
                fiz = fiz - grid_vec( 3 ) * g_val
                ! Stress term - TOTALLY UNTESTED AND PROBABLY INCOMPLETE
                ! Need short range term due to differentiation of coulomb operator,
                ! and SIC term
!!$                Do i_beta = 1, 3
!!$                   Do i_alpha = i_beta, 3
!!$                      s = - 2.0_wp * g_val * alpha * alpha * grid_vec( i_alpha ) * grid_vec( i_beta )
!!$                      stress( i_alpha, i_beta ) = stress( i_alpha, i_beta ) + s
!!$                   End Do
!!$                End Do
             End Do

             g_r = g_r0( 2 ) * f_rdr0( 2 ) * g_dr( 2 )
             g_r0( 2 ) = g_r

             f_rdr = f_rdr0 * g_drdr( :, 2 )
             f_rdr0 = f_rdr

          End Do

          g_r = g_r0( 1 ) * f_rdr00( 3 ) * g_dr( 3 )
          g_r0 = g_r

          f_rdr = f_rdr00 * g_drdr( :, 3 )
          f_rdr0 = f_rdr
          f_rdr00 = f_rdr

       End Do

       ! Apply appropriate scalings to force and energy
       ei   (    i ) = ei   (    i ) * 0.5_wp * qi_norm
!!$       force( :, i ) = force( :, i ) * 2.0_wp * alpha * alpha
       force( :, i ) = [ fix, fiy, fiz ] * 2.0_wp * alpha * alpha  * qi_norm
       
    End Do particle_loop
    !$omp end do
    !$omp end parallel

    ! Stress - Untested!
    ! Symmetrize stress tensor
!!$    Do i_beta = 1, 2
!!$       Do i_alpha = i_beta + 1, 3
!!$          stress( i_beta, i_alpha ) = stress( i_alpha, i_beta )
!!$       End Do
!!$    End Do

  End Subroutine charge_grid_forces

  Subroutine charge_grid_find_range( l, alpha, n_grid, range_gauss )

    ! DEPRECATED!!!!  -base all on range_gauss now, not n_grid
    ! i.e. invert the argument

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

    Use lattice_module, Only : lattice

    Implicit None

    Type( lattice ),              Intent( In    ) :: l
    Real( wp ),                   Intent( In    ) :: alpha
    Integer,    Dimension( 1:3 ), Intent( In    ) :: n_grid
    Integer,    Dimension( 1:3 ), Intent(   Out ) :: range_gauss

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

  Subroutine charge_grid_get_n_grid( l, alpha, range_gauss, gauss_tol, n_grid )

    Type( lattice ),              Intent( In    ) :: l
    Real( wp ),                   Intent( In    ) :: alpha
    Integer,                      Intent( In    ) :: range_gauss
    Real( wp ),                   Intent( In    ) :: gauss_tol
    Integer,    Dimension( 1:3 ), Intent(   Out ) :: n_grid

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

  End Subroutine charge_grid_get_n_grid

  function g(r, alpha)
    Real( wp ) :: g
    Real( wp ), Dimension(3), Intent( In    ) :: r
    Real( wp ), Intent( In    ) :: alpha

    g = exp(-alpha**2 * dot_product(r, r))

  end function g

  function f(r, dr, alpha)
    Real( wp ) :: f
    Real( wp ), Dimension(3), Intent( In    ) :: r
    Real( wp ), Dimension(3), Intent( In    ) :: dr
    Real( wp ), Intent( In    ) :: alpha

    f = exp(-2.0_wp * alpha**2 * dot_product(r, dr))

  end function f

End Module charge_grid_module
