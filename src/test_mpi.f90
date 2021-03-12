Program test_mpi

  !!$ Use omp_lib
  Use mpi_f08, Only : mpi_comm, mpi_init, mpi_finalize, mpi_comm_world, mpi_bcast, mpi_double_precision, mpi_integer, &
       mpi_dims_create, mpi_cart_create, mpi_cart_coords, mpi_allreduce, mpi_sum, mpi_in_place, mpi_character, mpi_comm_rank, &
       mpi_comm_size

  Use constants,                                 Only : wp
  Use lattice_module,                            Only : lattice
  Use charge_grid_module,                        Only : charge_grid_get_n_grid
  Use fast_fourier_poisson_module,               Only : ffp_long_range, ffp_sic
!!$  Use symetrically_screened_poisson_module,      Only : ssp_long_range, ssp_sic
  Use grid_io_module,                            Only : grid_io_save
  Use domains_module,                            Only : domain_build, domain_halo_build, domain_get_params
  Use comms_parallel_module,                     Only : comms_parallel
  Use halo_parallel_module,                      Only : halo_parallel_setter
  Use quadrature_trapezium_rule_module,          Only : quadrature_trapezium_rule
  Use FD_Laplacian_3d_module,                    Only : FD_Laplacian_3D
  Use equation_solver_base_class_module,         Only : equation_solver_base_class
  Use equation_solver_minres_module,             Only : equation_solver_minres
  Use equation_solver_conjugate_gradient_module, Only : equation_solver_conjugate_gradient
  Use Ewald_3d_module,                           Only : Ewald_3d_recipe, Ewald_3d_status

  Implicit None

!  Class( equation_solver_base_class ), Allocatable :: solver

  Type( lattice  ) :: l
  Type( mpi_comm ) :: cart_comm
!  Type( halo_parallel_setter          ) :: fd_swapper
!  Type( halo_parallel_setter          ) :: pot_swapper
!  Type( quadrature_trapezium_rule     ) :: grid_integrator
!  Type( comms_parallel                ) :: comms
!  Type( FD_Laplacian_3D               ) :: FD
  Type( Ewald_3d_recipe               ) :: ewald_recipe
  Type( Ewald_3d_status               ) :: status

  Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid_ffp
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_ffp
  Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid_ssp
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_ssp
  Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid_full

  Real( wp ), Dimension( :, : ), Allocatable :: r
  Real( wp ), Dimension( :, : ), Allocatable :: velocity
  Real( wp ), Dimension( :, : ), Allocatable :: dlp_force
  Real( wp ), Dimension( :, : ), Allocatable :: r_domain
  Real( wp ), Dimension( :, : ), Allocatable :: r_halo
  Real( wp ), Dimension( :, : ), Allocatable :: force_ffp
  Real( wp ), Dimension( :, : ), Allocatable :: force_ssp
  Real( wp ), Dimension( :, : ), Allocatable :: force_full
  Real( wp ), Dimension( :, : ), Allocatable :: dr

  Real( wp ), Dimension( 1:3, 1:3 ) :: a
!  Real( wp ), Dimension( 1:3, 1:3 ) :: dGrid_vecs
  Real( wp ), Dimension( 1:3, 1:3 ) :: stress

  Real( wp ), Dimension( : ), Allocatable :: q
  Real( wp ), Dimension( : ), Allocatable :: mass
  Real( wp ), Dimension( : ), Allocatable :: q_domain
  Real( wp ), Dimension( : ), Allocatable :: q_halo
  Real( wp ), Dimension( : ), Allocatable :: ei_ffp
  Real( wp ), Dimension( : ), Allocatable :: ei_ssp

  Real( wp ), Dimension( 1:3 ) :: dg
  Real( wp ), Dimension( 1:3 ) :: t
  Real( wp ), Dimension( 1:3 ) :: f_ssp

  Real( wp ) :: alpha, rtol
  Real( wp ) :: gauss_tol
  Real( wp ) :: xshift
  Real( wp ) :: recip_E_ffp
!  Real( wp ) :: rnorm
  Real( wp ) :: t_grid_ffp, t_pot_solve_ffp
  Real( wp ) :: recip_E_ssp
  Real( wp ) :: ei_full
  !  Real( wp ) :: t_grid_ssp, t_pot_solve_ssp, t_forces_ssp
  Real( wp ) :: dt = 0.001_wp

  Integer, Dimension( : ), Allocatable :: id, id_domain

  Integer, Dimension( 1:3 ) :: n_grid
  Integer, Dimension( 1:3 ) :: np_grid, np_grid_serial
  Integer, Dimension( 1:3 ) :: p_coords
  Integer, Dimension( 1:3 ) :: n_grid_domain
  Integer, Dimension( 1:3 ) :: domain_base_coords
  Integer, Dimension( 1:3 ) :: domain_end_coords

  Integer :: range_gauss, FD_order
  Integer :: nproc, me
  Integer :: nproc_cart, me_cart
  Integer :: level
  Integer :: n
  Integer :: n_at_loc, n_at
!  Integer :: itn, istop
  Integer :: which_solver
  Integer :: error
  Integer :: i

  Integer :: config_unit

!  Character( Len = 132 ) :: istop_message
  Character( Len = 132 ) :: what

  Call mpi_init( error )
  Call mpi_comm_size( mpi_comm_world, nproc )
  Call mpi_comm_rank( mpi_comm_world, me    )
  Write( *, * ) 'Hello from ', me, ' of ', nproc

  Read_params_on_proc_0: If( me == 0 ) Then

     Write( *, * ) 'Ewald param (alpha)?'
     Read ( *, * ) alpha
     Write( *, * ) alpha
     Write( *, * ) 'Number of points for gaussians?'
     Read ( *, * ) range_gauss
     Write( *, * ) range_gauss
     Write( *, * ) 'Cut off tolerance for gaussians?'
     Read ( *, * ) gauss_tol
     Write( *, * ) gauss_tol
     Write( *, * ) 'Residual tolerance for solver?'
     Read ( *, * ) rtol
     Write( *, * ) rtol
     Write( *, * ) 'FD_order?'
     Read ( *, * ) FD_order
     Write( *, * ) FD_order
     Write( *, * ) 'xshift?'
     Read ( *, * ) xshift
     Write( *, * ) xshift
     Write( *, * ) 'which solver?'
     Read ( *, * ) what
     Write( *, * ) what

     Select Case( Trim( Adjustl( what ) ) )
     Case( 'CG' )
        which_solver = 1
     Case( 'MINRES' )
        which_solver = 2
     Case( 'WJAC' )
        which_solver = 3
     Case( 'PFMG' )
        which_solver = 4
     Case Default
        Error Stop 'Unrecongnized solver'
     End Select

     Open( newunit=config_unit, File = 'CONFIG', status='OLD' )
     Call read_header( config_unit, n, level, a )
     Write( *, * ) 'n = ', n

  End If Read_params_on_proc_0

  Call mpi_bcast( alpha,                 1, mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( range_gauss,           1, mpi_integer,          0, mpi_comm_world )
  Call mpi_bcast( gauss_tol,             1, mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( rtol,                  1, mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( FD_order,              1, mpi_integer,          0, mpi_comm_world )
  Call mpi_bcast( n,                     1, mpi_integer,          0, mpi_comm_world )
  Call mpi_bcast( a,             Size( a ), mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( which_solver,          1, mpi_integer,          0, mpi_comm_world )

  Call mpi_bcast( what, Len( what ), mpi_character, 0, mpi_comm_world )

!!$  Select Case( which_solver )
!!$  Case( 1 )
!!$     Allocate( equation_solver_conjugate_gradient :: solver )
!!$  Case( 2 )
!!$     Allocate( equation_solver_minres :: solver )
!!$  End Select

  Call l%initialise( 3, a, alpha )

  Allocate( q   ( 1:n ) )
  Allocate( mass( 1:n ) )
  
  Allocate( r        ( 1:3, 1:n ) )
  Allocate( velocity ( 1:3, 1:n ) )
  Allocate( dlp_force( 1:3, 1:n ) )

  Read_CONFIG_on_proc_0: If( me == 0 ) Then

     Call read_config( level, config_unit, q, r, mass, velocity, dlp_force )
     Close( config_unit )

     ! Shift and move into reference cell
     Write( *, * ) 'r(1) before shift and reference = ', r( :, 1 )
     Do i = 1, n
        t = r( :, i )
        t( 1 ) = t( 1 ) + xshift
        Call l%to_reference( t, r( :, i ) )
     End Do
     Write( *, * ) 'r(1) after  shift and reference = ', r( :, 1 )

  End If Read_CONFIG_on_proc_0

  Call mpi_bcast( q        , Size( q         ), mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( mass     , Size( mass      ), mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( r        , Size( r         ), mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( velocity , Size( velocity  ), mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( dlp_force, Size( dlp_force ), mpi_double_precision, 0, mpi_comm_world )

  Allocate( id( 1:n ) )
  id = [ ( i, i = 1, n ) ]

  If( me == 0 ) Then
     Write( *, '( a, 3( t25, 3( f10.5, 1x ) / ) )' ) &
          'Lattice vectors: ', l%get_direct_vectors()
  End If

  ! Find the size of the grid
  ! Find the minimum size of grid commensurate with the range of gauss required
  Call charge_grid_get_n_grid( l, alpha, range_gauss, gauss_tol, n_grid )
  ! Now make the charge/pot grid commensurate with the process grid
  ! First factor the procs into a grid
  np_grid = 0
  Call mpi_dims_create( nproc, 3, np_grid )
  ! Now increase the size of the grid where it is not a multiple of the number of procs
  Where( Mod( n_grid, np_grid ) /= 0 )
     n_grid = ( n_grid / np_grid + 1 ) * np_grid
  End Where
  If( me == 0 ) Then
     Write( *, * ) 'N_grid = ', n_grid
  End If

  Do i = 1, 3
     dG( i ) = Sqrt( Dot_product( a( :, i ), a( :, i ) ) ) / n_grid( i )
  End Do

  ! First do the calculation in serial on proc 0
  Reference_serial: If( me == 0 ) Then
     ! Build the arrays of particles in the domain and halo
     np_grid_serial = [ 1, 1, 1 ]
     Call domain_build( l, q, r, n_grid, np_grid_serial, [ 0, 0, 0 ], q_domain, r_domain )
     Call domain_halo_build( l, q, r, n_grid, np_grid_serial, [ 0, 0, 0 ], [ range_gauss, range_gauss, range_gauss ], &
          q_halo, r_halo )

     ! Calculate the long range term by fourier space methods
     Allocate( q_grid_ffp( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
     Allocate( pot_grid_ffp( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
     Allocate( ei_ffp( 1:n ) )
     Allocate( force_ffp( 1:3, 1:n ) )
     Call ffp_long_range( l, q, r, alpha, FD_order, q_halo, r_halo, &
          recip_E_ffp, q_grid_ffp, pot_grid_ffp, ei_ffp, force_ffp, t_grid_ffp, t_pot_solve_ffp, error )
     ! Report the energy
     Write( *, '( "Serial FFP reciprocal energy: ", g24.14 )' ) recip_E_ffp
     ! Save the charge grid
     Call grid_io_save( 11, 'q_grid_FFP_serial.dat', l, q_grid_ffp )
     ! Save the FFP potential
     Call grid_io_save( 11, 'pot_grid_FFP_serial.dat', l, pot_grid_ffp )
     ! Save the forces
     Open( 11, file = 'forces_ffp_serial.dat' )
     Write( 11, * ) n, '     #number of particles'
     Write( 11, * ) Sum( force_ffp( 1, : ) ), Sum( force_ffp( 2, : ) ), Sum( force_ffp( 3, : ) ), '     #Nett force'
     Do i = 1, n
        Write( 11, * ) i, force_ffp( :, i )
     End Do
     Close( 11 )
     ! Report some potential  data on some potential problems
     Write( *, '( "FFP: Sum of charge over grid: ", g24.14 )' ) Sum( q_grid_ffp )
     Write( *, '( "FFP: Sum over pot grid      : ", g24.14 )' ) Sum( pot_grid_ffp )
     Write( *, '( "FFP: Nett force: ", 3( g20.14, 1x ) )' ) &
          Sum( force_ffp( 1, : ) ), Sum( force_ffp( 2, : ) ), Sum( force_ffp( 3, : ) )
     ! timing data
     Write( *, '( "FFP grid  time: ", f7.3 )' ) t_grid_ffp
     Write( *, '( "FFP solve time: ", f7.3 )' ) t_pot_solve_ffp
     ! Tidy up
     Deallocate( force_ffp )
     Deallocate( ei_ffp )
     Deallocate( pot_grid_ffp )
     Deallocate( q_grid_ffp )
  End If Reference_serial

  ! Now try and do it in parallel
  ! Factor the procs into a grid
  np_grid = 0
  Call mpi_dims_create( nproc, 3, np_grid )

  ! Create a cartesian communicator
  Call mpi_cart_create( mpi_comm_world, 3, np_grid, [ .True., .True., .True. ], .True., &
       cart_comm )
  Call mpi_comm_size( cart_comm, nproc_cart, error )
  Call mpi_comm_rank( cart_comm, me_cart,    error )
  If( me_cart == 0 ) Then
     Write( *, '( "Process grid dimensions: ", 3( i3, 1x ) )' ) np_grid
  End If

  ! Get the coordinates for this processor in the process grid
  Call mpi_cart_coords( cart_comm, me_cart, 3, p_coords )

  ! Get which bit of the grid I own
  Call domain_get_params( n_grid, np_grid, p_coords, n_grid_domain, domain_base_coords )
  domain_end_coords = domain_base_coords + n_grid_domain - 1

  ! Build the domain and halo
  Call domain_build( l, q, r, n_grid, np_grid, p_coords, q_domain, r_domain, id, id_domain )
  Call domain_halo_build( l, q, r, n_grid, np_grid, p_coords, [ range_gauss, range_gauss, range_gauss ], &
       q_halo, r_halo )
  ! Check total number of atoms is conserved
  n_at_loc = Size( q_domain )
  Call mpi_allreduce( n_at_loc, n_at, 1, mpi_integer, mpi_sum, cart_comm )
  If( me_cart == 0 ) Then
     Write( *, * ) 'n at = ', n_at
  End If

  Allocate( ei_ssp( 1:n_at_loc ) )
  Allocate( force_ssp( 1:3, 1:n_at_loc ) )

  ! START USE OF FULL INTERFACE !!!!!!!!!!!!!!
  ! ------------------------------------------

  ! Need to set these up earlier so use consistent sizes throughout calculation
  ! Set the recipe for the calculation
  Call ewald_recipe%mix( l, alpha, error, communicator = cart_comm%mpi_val, equation_solver = what, &
       residual_tol = rtol, FD_order = FD_order, range_gauss = range_gauss, gauss_tol = gauss_tol )
  ! May have changed some dimnsions due to constrains of parallel implementation
  ! Need these for checking potential grids etc.
  Call ewald_recipe%get_ingredients( n_grid = n_grid, range_gauss = range_gauss, &
       domain_base_coords = domain_base_coords, &
       domain_end_coords = domain_end_coords )
  If( me_cart == 0 ) Then
     Write( *, * )
     Write( *, * ) 'N_grid in full interface ', n_grid
  End If

  recip_E_ssp   = 0.0_wp
  Call ewald_recipe%consume(  q_domain, r_domain, q_halo, r_halo,  &
       recip_E_ssp, force_ssp, stress, error,          &
       q_grid = q_grid_ssp, pot_grid = pot_grid_ssp, ei = ei_ssp, status = status )
  Do i = 1, 3
     f_ssp( i ) = Sum( force_ssp( i, : ) )
  End Do
  Call mpi_allreduce( mpi_in_place, f_ssp, Size( f_ssp ), mpi_double_precision, mpi_sum, cart_comm )
  If( me_cart == 0 ) Then
     Write( *, * )
     Write( *, * ) 'NEW Energy: ', recip_E_ssp, recip_E_ssp - recip_E_ffp
     Write( *, * ) 'iterations = ', status%solver_iterations
     Write( *, * ) 'terminated because ', Trim( Adjustl( status%solver_stop_message ) )
     Write( *, * ) 'norm of error vector ', status%solver_residual_norm
     Write( *, '( "SSP grid   time: ", f7.3 )' ) status%t_grid
     Write( *, '( "SSP solve  time: ", f7.3 )' ) status%t_pot_solve
     Write( *, '( "SSP forces time: ", f7.3 )' ) status%t_forces
     Write( *, * ) 'Nett force: ', f_ssp
  End If
  Allocate( q_grid_full( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  q_grid_full = 0.0_wp
  q_grid_full( domain_base_coords( 1 ):domain_end_coords( 1 ), &
       domain_base_coords( 2 ):domain_end_coords( 2 ), &
       domain_base_coords( 3 ):domain_end_coords( 3 ) ) = q_grid_ssp
  Call mpi_allreduce( mpi_in_place, q_grid_full, Size( q_grid_full ), mpi_double_precision, mpi_sum, cart_comm )
  If( me_cart == 0 ) Then
     Call grid_io_save( 11, 'q_grid_ssp_parallel_new.dat', l, q_grid_full )
  End If
  q_grid_full = 0.0_wp
  q_grid_full( domain_base_coords( 1 ):domain_end_coords( 1 ), &
       domain_base_coords( 2 ):domain_end_coords( 2 ), &
       domain_base_coords( 3 ):domain_end_coords( 3 ) ) = pot_grid_ssp
  Call mpi_allreduce( mpi_in_place, q_grid_full, Size( q_grid_full ), mpi_double_precision, mpi_sum, cart_comm )
  If( me_cart == 0 ) Then
     Call grid_io_save( 11, 'pot_grid_ssp_parallel_new.dat', l, q_grid_full )
  End If

  Allocate( force_full( 1:3, 1:n ) )
  force_full = 0.0_wp
  Do i = 1, Size( force_ssp, Dim = 2 )
     force_full( :, id_domain( i ) ) = force_ssp( :, i )
  End Do
  Call mpi_allreduce( mpi_in_place, force_full, Size( force_full ), mpi_double_precision, mpi_sum, cart_comm )
  ei_full = Sum( ei_ssp )
  Call mpi_allreduce( mpi_in_place, ei_full, 1, mpi_double_precision, mpi_sum, cart_comm )
  If( me_cart == 0 ) Then
     Open( 11, file = 'forces_ssp_parallel_new.dat' )
     Write( 11, * ) n, '     #number of particles'
     Write( 11, * ) Sum( force_full( 1, : ) ), Sum( force_full( 2, : ) ), Sum( force_full( 3, : ) ), '     #Nett force'
     Do i = 1, n
        Write( 11, * ) i, force_full( :, i )
     End Do
     Close( 11 )
     Write( *, * ) 'Energy from summing per particle contributions ', ei_full
  End If
  
  ! See if we get the same answer with a restart
  If( me == 0 ) Then
     Write( *, * )
     Write( *, * ) 'Restart without moving the particles'
     Write( *, * )
  End If
  recip_E_ssp   = 0.0_wp
  Call ewald_recipe%consume(  q_domain, r_domain, q_halo, r_halo,  &
       recip_E_ssp, force_ssp, stress, error, &
       q_grid_old = q_grid_ssp, pot_grid_old = pot_grid_ssp, ei = ei_ssp, status = status )
  Do i = 1, 3
     f_ssp( i ) = Sum( force_ssp( i, : ) )
  End Do
  Call mpi_allreduce( mpi_in_place, f_ssp, Size( f_ssp ), mpi_double_precision, mpi_sum, cart_comm )
  If( me_cart == 0 ) Then
     Write( *, * )
     Write( *, * ) 'OLD Energy: ', recip_E_ssp, recip_E_ssp - recip_E_ffp
     Write( *, * ) 'iterations = ', status%solver_iterations
     Write( *, * ) 'terminated because ', Trim( Adjustl( status%solver_stop_message ) )
     Write( *, * ) 'norm of error vector ', status%solver_residual_norm
     Write( *, '( "SSP grid   time: ", f7.3 )' ) status%t_grid
     Write( *, '( "SSP solve  time: ", f7.3 )' ) status%t_pot_solve
     Write( *, '( "SSP forces time: ", f7.3 )' ) status%t_forces
     Write( *, * ) 'Nett force: ', f_ssp
  End If

  ! Now move the atoms slightly ...
  ! Can use the velocities and forces in the DLP CONFIG file to
  ! approximate a time step
  If( me == 0 ) Then
     Write( *, * )
     Write( *, * )
     Write( *, * ) 'Moving the particles! Using previous result as initial guess'
     Write( *, * )
  End If
  
  Allocate( dr, mold = r )
  ! Approximate time step
  dr = dt * velocity
  Do  i = 1, n
     dr( :, i ) = dr( :, i ) + dt * dt * dlp_force( :, i ) / ( 2.0_wp * mass( i ) )
  End Do
  r = r + dr

  ! PBC's
  Do i = 1, n
     t = r( :, i )
     Call l%to_reference( t, r( :, i ) )
  End Do

  ! Generate a new refernece energy for comparison
  Moved_serial: If( me == 0 ) Then
     ! Build the arrays of particles in the domain and halo
     np_grid_serial = [ 1, 1, 1 ]
     Call domain_build( l, q, r, n_grid, np_grid_serial, [ 0, 0, 0 ], q_domain, r_domain )
     Call domain_halo_build( l, q, r, n_grid, np_grid_serial, [ 0, 0, 0 ], [ range_gauss, range_gauss, range_gauss ], &
          q_halo, r_halo )

     ! Calculate the long range term by fourier space methods
     Allocate( q_grid_ffp( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
     Allocate( pot_grid_ffp( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
     Allocate( ei_ffp( 1:n ) )
     Allocate( force_ffp( 1:3, 1:n ) )
     Call ffp_long_range( l, q, r, alpha, FD_order, q_halo, r_halo, &
          recip_E_ffp, q_grid_ffp, pot_grid_ffp, ei_ffp, force_ffp, t_grid_ffp, t_pot_solve_ffp, error )
     ! Report the energy
     Write( *, '( "Serial FFP reciprocal energy: ", g24.14 )' ) recip_E_ffp
     ! Save the charge grid
     Call grid_io_save( 11, 'q_grid_FFP_serial_moved.dat', l, q_grid_ffp )
     ! Save the FFP potential
     Call grid_io_save( 11, 'pot_grid_FFP_serial_moved.dat', l, pot_grid_ffp )
     ! Save the forces
     Open( 11, file = 'forces_ffp_serial_moved.dat' )
     Write( 11, * ) n, '     #number of particles'
     Write( 11, * ) Sum( force_ffp( 1, : ) ), Sum( force_ffp( 2, : ) ), Sum( force_ffp( 3, : ) ), '     #Nett force'
     Do i = 1, n
        Write( 11, * ) i, force_ffp( :, i )
     End Do
     Close( 11 )
     ! Report some potential  data on some potential problems
     Write( *, '( "FFP: Sum of charge over grid: ", g24.14 )' ) Sum( q_grid_ffp )
     Write( *, '( "FFP: Sum over pot grid      : ", g24.14 )' ) Sum( pot_grid_ffp )
     Write( *, '( "FFP: Nett force: ", 3( g20.14, 1x ) )' ) &
          Sum( force_ffp( 1, : ) ), Sum( force_ffp( 2, : ) ), Sum( force_ffp( 3, : ) )
     ! timing data
     Write( *, '( "FFP grid  time: ", f7.3 )' ) t_grid_ffp
     Write( *, '( "FFP solve time: ", f7.3 )' ) t_pot_solve_ffp
     ! Tidy up
     Deallocate( force_ffp )
     Deallocate( ei_ffp )
     Deallocate( pot_grid_ffp )
     Deallocate( q_grid_ffp )
  End If Moved_serial

  ! Build the domain and halo
  Call domain_build( l, q, r, n_grid, np_grid, p_coords, q_domain, r_domain, id, id_domain )
  Call domain_halo_build( l, q, r, n_grid, np_grid, p_coords, [ range_gauss, range_gauss, range_gauss ], &
       q_halo, r_halo )
  Deallocate( ei_ssp )
  Deallocate( force_ssp )
  n_at_loc = Size( q_domain )
  Allocate( ei_ssp( 1:n_at_loc ) )
  Allocate( force_ssp( 1:3, 1:n_at_loc ) )
  recip_E_ssp   = 0.0_wp
  Call ewald_recipe%consume(  q_domain, r_domain, q_halo, r_halo,  &
       recip_E_ssp, force_ssp, stress, error, &
       q_grid_old = q_grid_ssp, pot_grid_old = pot_grid_ssp, ei = ei_ssp, status = status )
  Do i = 1, 3
     f_ssp( i ) = Sum( force_ssp( i, : ) )
  End Do
  Call mpi_allreduce( mpi_in_place, f_ssp, Size( f_ssp ), mpi_double_precision, mpi_sum, cart_comm )
  If( me_cart == 0 ) Then
     Write( *, * )
     Write( *, * ) 'SHIFTED Energy: ', recip_E_ssp, recip_E_ssp - recip_E_ffp
     Write( *, * ) 'iterations = ', status%solver_iterations
     Write( *, * ) 'terminated because ', Trim( Adjustl( status%solver_stop_message ) )
     Write( *, * ) 'norm of error vector ', status%solver_residual_norm
     Write( *, '( "SSP grid   time: ", f7.3 )' ) status%t_grid
     Write( *, '( "SSP solve  time: ", f7.3 )' ) status%t_pot_solve
     Write( *, '( "SSP forces time: ", f7.3 )' ) status%t_forces
     Write( *, * ) 'Nett force: ', f_ssp
  End If

  Call mpi_finalize( error )

Contains

  Subroutine read_header( unit, n, level, vecs )

    Integer,                       Intent( In    ) :: unit
    Integer,                       Intent(   Out ) :: n
    Integer,                       Intent(   Out ) :: level
    Real( wp ), Dimension( :, : ), Intent(   Out ) :: vecs

    Integer :: junk

    Read( unit, * )
    Read( unit, * ) level, junk, n
    Read( unit, * ) vecs

  End Subroutine read_header

  Subroutine read_config( level, unit, q, r, mass, velocity, force )

    Use, Intrinsic :: ieee_arithmetic, Only : ieee_value, ieee_signaling_nan
    
    Integer,                       Intent( In    ) :: level
    Integer,                       Intent( In    ) :: unit
    Real( wp ), Dimension( :    ), Intent(   Out ) :: q
    Real( wp ), Dimension( :, : ), Intent(   Out ) :: r
    Real( wp ), Dimension(    : ), Intent(   Out ) :: mass
    Real( wp ), Dimension( :, : ), Intent(   Out ) :: velocity
    Real( wp ), Dimension( :, : ), Intent(   Out ) :: force

    Integer :: n
    Integer :: i

    Character( Len = 3 ) :: what

    velocity = ieee_value( velocity, ieee_signaling_nan )
    force    = ieee_value( force   , ieee_signaling_nan )
    
    n = Size( q )

    Do i = 1, n
       Read( unit, * ) what
       If( what == 'Na+' ) Then
          q( i )    = +1.0_wp
          mass( i ) = 22.9898_wp 
       Else If( what == 'Cl-' ) Then
          q( i )    = -1.0_wp
          mass( i ) = 35.4530_wp
       Else
          Write( *, * ) 'what is a ', what, ' ?'
          Stop 'Error in CONFIG reading'
       End If
       Read( unit, * ) r( :, i )
       If( level >= 1 ) Then
          Read( unit, * ) velocity( :, i )
       End If
       If( level >= 2 ) Then
          Read( unit, * ) force( :, i )
       End If
    End Do

  End Subroutine read_config

End Program test_mpi
