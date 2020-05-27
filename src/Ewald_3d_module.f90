Module Ewald_3d_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

  Implicit None
  
  Public :: Ewald_3d

  Type, Public :: Ewald_3d_status
!!$     Private
     Real( wp )              :: t_grid                  = - Huge( 1.0_wp )
     Real( wp )              :: t_pot_solve             = - Huge( 1.0_wp )
     Real( wp )              :: t_forces                = - Huge( 1.0_wp )
     Integer                 :: solver_iterations       = - Huge( 1      )
     Integer                 :: solver_stop_code        = - Huge( 1      )
     Character( Len = 132 )  :: solver_stop_message     = "How should I know?"
     Real( wp )              :: solver_residual_norm    = - Huge( 1.0_wp )
  End type Ewald_3d_status
  
  Private

  ! Default values for the methods parameters
  Real( wp )          , Parameter :: gauss_tol_default    = 1e-15_wp ! Value for cutting of the gaussian
  Integer             , Parameter :: range_gauss_default  = 14       ! Number of points to grid one side of a gaussian
  Integer             , Parameter :: FD_order_default     = 14       ! Default order of the finite difference approximation
!!$  Character( Len = * ), Parameter :: default_solver       = "minres" ! Default equation solver 
  Character( Len = * ), Parameter :: default_solver       = "CG" ! Default equation solver 
  Real( wp )          , Parameter :: residual_tol_default = 1e-10_wp ! Tolerance on residual in equation solver
  
Contains

  Subroutine Ewald_3d(  l, alpha, q, r, q_halo, r_halo, &
       recip_E, forces, stress, error,                  &
       communicator, method, q_grid_old, pot_grid_old, ei, q_grid, pot_grid, status ) 

    Use lattice_module                   , Only : lattice

    Use comms_base_class_module          , Only : comms_base_class
    Use comms_serial_module              , Only : comms_serial
    Use comms_parallel_module            , Only : comms_parallel

    Use charge_grid_module               , Only : charge_grid_calculate, charge_grid_get_n_grid, charge_grid_forces

    Use FD_template_module               , Only : FD_template
    Use FD_Laplacian_3d_module           , Only : FD_Laplacian_3D

    Use halo_setter_base_module          , Only : halo_setter_base_class
    Use halo_serial_module               , Only : halo_serial_setter
    Use halo_parallel_module             , Only : halo_parallel_setter

    Use quadrature_base_module           , Only : quadrature_base_class
    Use quadrature_trapezium_rule_module , Only : quadrature_trapezium_rule
    
    Use equation_solver_base_class_module        , Only : equation_solver_base_class
    Use equation_solver_minres_module            , Only : equation_solver_minres
    Use equation_solver_conjugate_gradient_module, Only : equation_solver_conjugate_gradient

    Use symetrically_screened_poisson_module, Only : ssp_long_range
    
    ! Should eventually get rid of this one once HALO_HACK is sorted
    ! But will need a more general solution to find base coords of domain
    Use domains_module, Only : domain_get_params

    Implicit None
    
    Type( lattice )                    ,              Intent( In    )           :: l
    Real( wp )                         ,              Intent( In    )           :: alpha
    Real( wp ), Dimension( 1:         ),              Intent( In    )           :: q
    Real( wp ), Dimension( 1:, 1:     ),              Intent( In    )           :: r
    Real( wp ), Dimension( 1:         ),              Intent( In    )           :: q_halo
    Real( wp ), Dimension( 1:, 1:     ),              Intent( In    )           :: r_halo
    Real( wp )                         ,              Intent(   Out )           :: recip_E
    Real( wp ), Dimension( 1:, 1:     ),              Intent(   Out )           :: forces
    Real( wp ), Dimension( 1:, 1:     ),              Intent(   Out )           :: stress
    Integer                            ,              Intent(   Out )           :: error
    Integer                            ,              Intent( In    ), Optional :: communicator
    Character( Len = * )               ,              Intent( In    ), Optional :: method
    Real( wp ), Dimension(  :,  :,  : ),              Intent( In    ), Optional :: q_grid_old
    Real( wp ), Dimension(  :,  :,  : ),              Intent( In    ), Optional :: pot_grid_old
    Real( wp ), Dimension(  :         ),              Intent(   Out ), Optional :: ei
    Real( wp ), Dimension(  :,  :,  : ), Allocatable, Intent(   Out ), Optional :: q_grid
    Real( wp ), Dimension(  :,  :,  : ), Allocatable, Intent(   Out ), Optional :: pot_grid
    Type( Ewald_3d_status )                         , Intent(   Out ), Optional :: status

    Class( comms_base_class           ), Allocatable :: comms
    Class( FD_template                ), Allocatable :: FD
    Class( halo_setter_base_class     ), Allocatable :: FD_swapper
    Class( halo_setter_base_class     ), Allocatable :: pot_swapper
    Class( equation_solver_base_class ), Allocatable :: solver
    Class( quadrature_base_class      ), Allocatable :: grid_integrator

    Type( Ewald_3d_status ) :: scope_status

    Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid_ssp
    Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_ssp

    Real( wp ), Dimension( 1:3, 1:3 ) :: Grid_vecs, dGrid_vecs
    
    Real( wp ), Dimension( : ), Allocatable :: ei_ssp

    Real( wp ) :: gauss_tol
    Real( wp ) :: residual_tol

    Integer, Dimension( 1:3 ) :: n_grid
    Integer, Dimension( 1:3 ) :: n_proc_grid
    Integer, Dimension( 1:3 ) :: me_proc_coords
    Integer, Dimension( 1:3 ) :: n_grid_domain
    Integer, Dimension( 1:3 ) :: domain_base_coords
    Integer, Dimension( 1:3 ) :: domain_end_coords
    
    Integer :: range_gauss
    Integer :: FD_order
    Integer :: max_swap
    Integer :: scope_communicator
    Integer :: i

    Character( Len = : ), Allocatable :: scope_method
    
    ! Stress not calculated yet
    stress = 0.0_wp

    ! Set up communication layer. Provide a communicator to indicate the run is MPI parallel
    If( Present( communicator ) ) Then
       Allocate( comms_parallel :: comms )
       scope_communicator = communicator
    Else
       Allocate( comms_serial :: comms )
       ! Any rubbish for the communicator will do if in serial
       scope_communicator = - Huge( scope_communicator )
    End If
    Call comms%set_comm( scope_communicator, error )
    If( error /= 0 ) Then
       ! Will rationalise error codes later
       Return
    End If

    ! Get the size of the proc grid, and where I am
    Call comms%get_proc_grid(  n_proc_grid )
    Call comms%get_proc_coords( me_proc_coords )

    ! Select the order for the finite difference approximation
    ! Select now as need now to set the grid size due to the hacks below involved in setting the grid size
    ! But can't init the FD operator itself yet as need to know the grid size for that.
    FD_order = FD_order_default

    ! Set the size of the grid
    ! First find the minimum consistent with input alpha and the current gauss sizing parameters
    range_gauss = range_gauss_default
    gauss_tol   = gauss_tol_default    
    Call charge_grid_get_n_grid( l, alpha, range_gauss, gauss_tol, n_grid )

    ! Due to the current domain mapping routines working on the grid and not the actual volumes of the
    ! domains make the grid consistent with the volumes. Simply making sure the grid size is a multiple
    ! of the procs down each dimension is enough
    ! WANT TO FIX THIS SO THE NEXT BLOCK GOES AWAY - but note next hack depends on this one!
    Domain_hack: Block
      Where( Mod( n_grid, n_proc_grid ) /= 0 )
         n_grid = ( n_grid / n_proc_grid + 1 ) * n_proc_grid
      End Where
    End Block Domain_hack
    
    ! Due to the halo swapping routines only coping with halos of maximum size one domain
    ! make the grid size consistent with this. An issue for small grids compared to proc count
    ! WANT TO FIX THIS SO THE NEXT BLOCK GOES AWAY
    Halo_hack: Block
      Do
         ! Find the minimum size of the domain with the current parameters
         ! We know all the domains are the same size due to the previous hack - NOTE DEPENDENCY -
         ! so don't need comms to find the minimum number of grdi points asscoiated with each proc
         Call domain_get_params( n_grid, n_proc_grid, me_proc_coords, n_grid_domain, domain_base_coords )
         ! Get the maximum size for halo swapping
         ! Note bug in FD%get_order returns half the order, which is actually the amount we need to swap
         max_swap = Max( ( 2 * FD_order ) / 2, range_gauss + 1 )
         ! Check max_swap is at most the size of the domain
         If( All( max_swap <= n_grid_domain ) ) Exit
         ! Domain size not big enough. Increase n_grid and try again. Keep n_grid a multiple of the
         ! proc grid size
         n_grid = n_grid + n_proc_grid
      End Do
    End Block Halo_hack

    ! Now have size of grid get size of domain and where our domain starts
    Call domain_get_params( n_grid, n_proc_grid, me_proc_coords, n_grid_domain, domain_base_coords )

    ! Set up Finite difference operator
    Grid_vecs = l%get_direct_vectors()
    Do i = 1, 3
       dGrid_vecs( :, i ) = Grid_vecs( :, i ) / n_grid( i )
    End Do
    Allocate( FD_Laplacian_3D :: FD )
    ! Hack ... All init methods need a rethink
    Select Type( FD )
    Class is ( FD_Laplacian_3D )
       Call FD%init( FD_order, dGrid_vecs )
    End Select

    ! Set up the halo swapper for the FD operator
    ! Is this condition better done by using the intrinsic extends_type_of of on comms? Probably ...
    If( Present( communicator ) ) Then
       Allocate( halo_parallel_setter :: FD_swapper )
    Else
       Allocate( halo_serial_setter :: FD_swapper )
    End If
    Select Type( FD_swapper )
    Class is ( halo_parallel_setter )
       Call FD_swapper%init ( n_grid_domain, FD%get_order(), scope_communicator, error )
    Class is ( halo_serial_setter )
       Call FD_swapper%init( error )
    End Select
    If( error /= 0 ) Then
       Return
    End If

    ! Set up potential halo swapper
    ! The +1 is to do with atoms near the edge of the grid where
    ! the nearest grid point to the centre of the gaussian actually lies in the next domain, but the atom
    ! is in this domain
    If( Present( communicator ) ) Then
       Allocate( halo_parallel_setter :: pot_swapper )
    Else
       Allocate( halo_serial_setter :: pot_swapper )
    End If
    Select Type( pot_swapper )
    Class is ( halo_parallel_setter )
       Call pot_swapper%init ( n_grid_domain, range_gauss + 1, scope_communicator, error )
    Class is ( halo_serial_setter )
       Call pot_swapper%init( error )
    End Select
    If( error /= 0 ) Then
       Return
    End If

    ! Set up the equation solver
    ! First decide on the accuracy
    residual_tol = residual_tol_default
    If( Present( method ) ) Then
       scope_method = method 
    Else
       scope_method = default_solver
    End If
    Select Case( Trim( Adjustl( scope_method ) ) )
    Case( "cg", "CG" )
       Allocate( equation_solver_conjugate_gradient :: solver )
    Case( "minres", "MINRES" )
       Allocate( equation_solver_minres :: solver )
    Case Default
       error = 3
       Return
    End Select

    ! Set up the quadrature method
    Allocate( quadrature_trapezium_rule :: grid_integrator )
    
    ! Allocate the charge and potential grids
    domain_end_coords = domain_base_coords + n_grid_domain - 1
    Allocate( q_grid_ssp( domain_base_coords( 1 ):domain_end_coords( 1 ), &
         domain_base_coords( 2 ):domain_end_coords( 2 ), &
         domain_base_coords( 3 ):domain_end_coords( 3 ) ) )
    Allocate( pot_grid_ssp( domain_base_coords( 1 ):domain_end_coords( 1 ), &
         domain_base_coords( 2 ):domain_end_coords( 2 ), &
         domain_base_coords( 3 ):domain_end_coords( 3 ) ) )

    ! Allocate the per particle energy
    Allocate( ei_ssp( 1:Size( q ) ) )

    ! If provided with old data can try and solve for a delta problem which should be much quicker
    ! for iterative solvers at least
!!$    If( Present( q_grid_old ) .And. Present(  pot_grid_old ) ) Then
    
    ! Solve the problem
    Call ssp_long_range( l, q, r, alpha, FD, q_halo, r_halo, range_gauss, n_grid, Lbound( q_grid_ssp ), residual_tol, &
       recip_E, q_grid_ssp, pot_grid_ssp, solver, comms, FD_swapper, pot_swapper, grid_integrator      , &
       ei_ssp, forces, scope_status%t_grid, scope_status%t_pot_solve, scope_status%t_forces,             &
       scope_status%solver_iterations, scope_status%solver_stop_code, scope_status%solver_stop_message , &
       scope_status%solver_residual_norm, error, q_grid_old, pot_grid_old )
!!$    If( error /= 0 ) Then
!!$       Return
!!$    End If

    ! Return optional arguments as required
    If( Present( q_grid ) ) Then
       Call Move_alloc( q_grid_ssp, q_grid )
    End If

    If( Present( pot_grid ) ) Then
       Call Move_alloc( pot_grid_ssp, pot_grid )
    End If

    If( Present( ei ) ) Then
       ei = ei_ssp
    End If

    If( Present( status ) ) Then
       status = scope_status
    End If

  End Subroutine Ewald_3d
    
End Module Ewald_3d_module
