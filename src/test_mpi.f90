Program test_mpi

  !$ Use omp_lib
  Use mpi_f08, Only : mpi_comm, mpi_init, mpi_finalize, mpi_comm_world, mpi_bcast, mpi_double_precision, mpi_integer, &
       mpi_dims_create, mpi_cart_create, mpi_cart_coords, mpi_allreduce, mpi_sum, mpi_in_place
  
  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

  Use lattice_module                      , Only : lattice
  Use charge_grid_module                  , Only : charge_grid_get_n_grid
  Use fast_fourier_poisson_module         , Only : ffp_long_range, ffp_sic
  Use symetrically_screened_poisson_module, Only : ssp_long_range, ssp_sic
!!$  Use real_space_module                   , Only : real_space_energy
  Use grid_io_module                      , Only : grid_io_save
  Use domains_module                      , Only : domain_build, domain_halo_build, domain_get_params
  Use comms_parallel_module               , Only : comms_parallel
  Use halo_parallel_module                , Only : halo_parallel_setter
!!$  Use halo_serial_module                , Only : halo_serial_setter
  Use quadrature_trapezium_rule_module    , Only : quadrature_trapezium_rule
  Use FD_Laplacian_3d_module              , Only : FD_Laplacian_3D
  
  Implicit None

  Type( lattice  ) :: l
  Type( mpi_comm ) :: cart_comm
  Type( halo_parallel_setter          ) :: fd_swapper
!!$  Type( halo_serial_setter          ) :: fd_swapper
  Type( halo_parallel_setter          ) :: pot_swapper
  Type( quadrature_trapezium_rule     ) :: grid_integrator
  Type( comms_parallel                ) :: comms
  Type( FD_Laplacian_3D               ) :: FD
  
  Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid_ffp
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_ffp
  Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid_ssp
  Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_ssp
  Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid_full
  
  Real( wp ), Dimension( :, : ), Allocatable :: r
  Real( wp ), Dimension( :, : ), Allocatable :: r_domain
  Real( wp ), Dimension( :, : ), Allocatable :: r_halo
  Real( wp ), Dimension( :, : ), Allocatable :: force_ffp
  Real( wp ), Dimension( :, : ), Allocatable :: force_ssp
  Real( wp ), Dimension( :, : ), Allocatable :: force_full

  Real( wp ), Dimension( 1:3, 1:3 ) :: a
  Real( wp ), Dimension( 1:3, 1:3 ) :: dGrid_vecs
  
  Real( wp ), Dimension( : ), Allocatable :: q
  Real( wp ), Dimension( : ), Allocatable :: q_domain
  Real( wp ), Dimension( : ), Allocatable :: q_halo
  Real( wp ), Dimension( : ), Allocatable :: ei_ffp
  Real( wp ), Dimension( : ), Allocatable :: ei_ssp

  Real( wp ), Dimension( 1:3 ) :: dg
  Real( wp ), Dimension( 1:3 ) :: t

  Real( wp ) :: alpha
  Real( wp ) :: gauss_tol
  Real( wp ) :: xshift
  Real( wp ) :: recip_E_ffp
  Real( wp ) :: rnorm
  Real( wp ) :: t_grid_ffp, t_pot_solve_ffp
  Real( wp ) :: recip_E_ssp
  Real( wp ) :: ei_full
  Real( wp ) :: t_grid_ssp, t_pot_solve_ssp, t_forces_ssp

  Integer, Dimension( : ), Allocatable :: id, id_domain
  
  Integer, Dimension( 1:3 ) :: n_grid
  Integer, Dimension( 1:3 ) :: np_grid
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
  Integer :: itn, istop
  Integer :: error
  Integer :: i

  Character( Len = 132 ) :: istop_message

  Call mpi_init( error )
  Call mpi_comm_size( mpi_comm_world, nproc )
  Call mpi_comm_rank( mpi_comm_world, me    )
  Write( *, * ) 'Hello from ', me, ' of ', nproc

  Read_params_on_proc_0: If( me == 0 ) Then

     Write( *, * ) 'Ewald param ?'
     Read ( *, * ) alpha
     Write( *, * ) 'Number of ppoints for gaussians?'
     Read ( *, * ) range_gauss
     Write( *, * ) 'Cut off tolerance for gaussians ?'
     Read ( *, * ) gauss_tol
     Write( *, * ) 'FD_order?'
     Read ( *, * ) FD_order
     Write( *, * ) 'xshift?'
     Read ( *, * ) xshift
     
     Open( 10, File = 'CONFIG' )
     Call read_header( 10, n, level, a )
     Write( *, * ) 'n = ', n

  End If Read_params_on_proc_0
     
  Call mpi_bcast( alpha      ,         1, mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( range_gauss,         1, mpi_integer         , 0, mpi_comm_world )
  Call mpi_bcast( gauss_tol  ,         1, mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( FD_order   ,         1, mpi_integer         , 0, mpi_comm_world )
  Call mpi_bcast( n          ,         1, mpi_integer         , 0, mpi_comm_world )
  Call mpi_bcast( a          , Size( a ), mpi_double_precision, 0, mpi_comm_world )

  Call l%initialise( 3, a, alpha )
  
  Allocate( q( 1:n ) )
  Allocate( r( 1:3, 1:n ) )

  Read_CONFIG_on_proc_0: If( me == 0 ) Then  
  
     Call read_config( level, 10, q, r )
     Close( 10 )
     
     ! Shift and move into reference cell
     Write( *, * ) 'r(1) before shift and reference = ', r( :, 1 )
     Do i = 1, n
        t = r( :, i )
        t( 1 ) = t( 1 ) + xshift
        Call l%to_reference( t, r( :, i ) ) 
     End Do
     Write( *, * ) 'r(1) after  shift and reference = ', r( :, 1 )

  End If Read_CONFIG_on_proc_0

  Call mpi_bcast( q, Size( q ), mpi_double_precision, 0, mpi_comm_world )
  Call mpi_bcast( r, Size( r ), mpi_double_precision, 0, mpi_comm_world )

  Allocate( id( 1:n ) )
  id = [ ( i, i = 1, n ) ]

  If( me == 0 ) Then
     Write( *, '( a, 3( t25, 3( f10.5, 1x ) / ) )' ) &
          'Lattice vectors: ', l%get_direct_vectors()
  End If
  
  Call charge_grid_get_n_grid( l, alpha, range_gauss, gauss_tol, n_grid )
  If( me == 0 ) Then
     Write( *, * ) 'N_grid = ', n_grid
  End If
  
  Do i = 1, 3
     dG( i ) = Sqrt( Dot_product( a( :, i ), a( :, i ) ) ) / n_grid( i )
  End Do

  ! First do the calculation in serial on proc 0
  Reference_serial: If( me == 0 ) Then
     ! Build the arrays of particles in the domain and halo
     np_grid = [ 1, 1, 1 ]
     Call domain_build( l, q, r, n_grid, np_grid, [ 0, 0, 0 ], q_domain, r_domain )
     Call domain_halo_build( l, q, r, n_grid, np_grid, [ 0, 0, 0 ], [ range_gauss, range_gauss, range_gauss ], &
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
  Call mpi_comm_rank( cart_comm, me_cart   , error )
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

  Allocate( q_grid_ssp( domain_base_coords( 1 ):domain_end_coords( 1 ), &
       domain_base_coords( 2 ):domain_end_coords( 2 ), &
       domain_base_coords( 3 ):domain_end_coords( 3 ) ) )
  Allocate( pot_grid_ssp( domain_base_coords( 1 ):domain_end_coords( 1 ), &
       domain_base_coords( 2 ):domain_end_coords( 2 ), &
       domain_base_coords( 3 ):domain_end_coords( 3 ) ) )
  Allocate( ei_ssp( 1:n_at_loc ) )
  Allocate( force_ssp( 1:3, 1:n_at_loc ) )
  dGrid_vecs = l%get_direct_vectors()
  Do i = 1, 3
     dGrid_vecs( :, i ) = dGrid_vecs( :, i ) / n_grid( i )
  End Do
  ! Set up comms
  Call comms%set_comm( cart_comm%mpi_val, error )
  If( error /= 0 ) Then
     Error stop "Error in setting comms communicator"
  End If
  ! set up finite difference operator
  Call FD%init( FD_order, dGrid_vecs )
!!$  ! Set up halo swapper for finit difference operator
  Call fd_swapper%init ( n_grid_domain, FD%get_order(), cart_comm, error )
!!$  Call fd_swapper%init ( error )
  If( error /= 0 ) Then
     Error stop "Error in setting comms communicator"
  End If
  ! Set up potential halo swapper
  ! Need to track down +2
  Call pot_swapper%init( n_grid_domain, range_gauss + 1, cart_comm, error )
!!$  Call pot_swapper%init( n_grid_domain, range_gauss, cart_comm, error )
  ! Solve for the long range using finiste difference
  Call ssp_long_range( l, q_domain, r_domain, alpha, FD, q_halo, r_halo, range_gauss, n_grid, Lbound( q_grid_ssp ), &
       recip_E_ssp, q_grid_ssp, pot_grid_ssp, comms, fd_swapper, pot_swapper, grid_integrator, &
       ei_ssp, force_ssp, t_grid_ssp, t_pot_solve_ssp, t_forces_ssp, itn, istop, istop_message, rnorm, error )
  ! Check charge grid
  If( me_cart == 0 ) Then
     Write( *, * ) 'SSP Energy: ', recip_E_ssp, recip_E_ssp - recip_E_ffp
     Write( *, * ) itn, istop, istop_message, rnorm
     Write( *, '( "SSP grid   time: ", f7.3 )' ) t_grid_ssp
     Write( *, '( "SSP solve  time: ", f7.3 )' ) t_pot_solve_ssp
     Write( *, '( "SSP forces time: ", f7.3 )' ) t_forces_ssp
  End If
  Allocate( q_grid_full( 0:n_grid( 1 ) - 1, 0:n_grid( 2 ) - 1, 0:n_grid( 3 ) - 1 ) )
  q_grid_full = 0.0_wp
  q_grid_full( domain_base_coords( 1 ):domain_end_coords( 1 ), &
       domain_base_coords( 2 ):domain_end_coords( 2 ), &
       domain_base_coords( 3 ):domain_end_coords( 3 ) ) = q_grid_ssp
  Call mpi_allreduce( mpi_in_place, q_grid_full, Size( q_grid_full ), mpi_double_precision, mpi_sum, cart_comm )
  If( me_cart == 0 ) Then
     Call grid_io_save( 11, 'q_grid_ssp_parallel.dat', l, q_grid_full )
  End If
  q_grid_full = 0.0_wp
  q_grid_full( domain_base_coords( 1 ):domain_end_coords( 1 ), &
       domain_base_coords( 2 ):domain_end_coords( 2 ), &
       domain_base_coords( 3 ):domain_end_coords( 3 ) ) = pot_grid_ssp
  Call mpi_allreduce( mpi_in_place, q_grid_full, Size( q_grid_full ), mpi_double_precision, mpi_sum, cart_comm )
  If( me_cart == 0 ) Then
     Call grid_io_save( 11, 'pot_grid_ssp_parallel.dat', l, q_grid_full )
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
     Open( 11, file = 'forces_ssp_parallel.dat' )
     Write( 11, * ) n, '     #number of particles'
     Write( 11, * ) Sum( force_full( 1, : ) ), Sum( force_full( 2, : ) ), Sum( force_full( 3, : ) ), '     #Nett force'
     Do i = 1, n
        Write( 11, * ) i, force_full( :, i )
     End Do
     Close( 11 )
     Write( *, * ) 'Energy from summing per particle contributions ', ei_full
  End If
  
  Call mpi_finalize( error )
  
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

End Program test_mpi
