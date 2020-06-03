Module Ewald_3d_module

  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Use lattice_module                          , Only : lattice
  Use comms_base_class_module                 , Only : comms_base_class
  Use FD_template_module                      , Only : FD_template
  Use halo_setter_base_module                 , Only : halo_setter_base_class
  Use quadrature_base_module                  , Only : quadrature_base_class
  Use equation_solver_precon_base_class_module, Only : equation_solver_precon_base_class

  Implicit None

  Integer, Parameter :: EWALD_3D_ISSUED_WARNING                = -1
  Integer, Parameter :: EWALD_3D_SUCCESS                       = 0
  Integer, Parameter :: EWALD_3D_INIT_COMMS_SETUP_FAILED       = 1
  Integer, Parameter :: EWALD_3D_INIT_FD_SWAPPER_SETUP_FAILED  = 2
  Integer, Parameter :: EWALD_3D_INIT_POT_SWAPPER_SETUP_FAILED = 3
  Integer, Parameter :: EWALD_3D_INIT_UNKNOWN_EQUATION_SOLVER  = 4
  Integer, Parameter :: EWALD_3D_SOLVER_FAILED                 = 5
  
  Type, Public :: Ewald_3d_recipe
     Private
     Type( lattice )                                         :: l
     Real( wp )                                              :: alpha
     Integer, Dimension( : )                   , Allocatable :: n_grid
     Integer, Dimension( : )                   , Allocatable :: domain_base_coords
     Integer, Dimension( : )                   , Allocatable :: domain_end_coords
     Integer                                                 :: range_gauss
     Real( wp )                                              :: gauss_tol
     Real( wp )                                              :: residual_tol
     Class( comms_base_class                  ), Allocatable :: comms
     Class( FD_template                       ), Allocatable :: FD
     Class( halo_setter_base_class            ), Allocatable :: FD_swapper
     Class( halo_setter_base_class            ), Allocatable :: pot_swapper
     Class( equation_solver_precon_base_class ), Allocatable :: solver
     Class( quadrature_base_class             ), Allocatable :: grid_integrator
   Contains
     Procedure, Public :: mix             => Ewald_3d_init
     Procedure, Public :: consume         => Ewald_3d 
     Procedure, Public :: get_ingredients => Ewald_3d_get_ingredients
  End type Ewald_3d_recipe

  Type, Public :: Ewald_3d_status
!!$     Private  ! Need to add inquiry routines
     Real( wp )              :: t_grid                  = - Huge( 1.0_wp )
     Real( wp )              :: t_pot_solve             = - Huge( 1.0_wp )
     Real( wp )              :: t_forces                = - Huge( 1.0_wp )
     Integer                 :: solver_iterations       = - Huge( 1      )
     Integer                 :: solver_stop_code        = - Huge( 1      )
     Character( Len = 132  ) :: solver_stop_message     = "How should I know?"
     Real( wp )              :: solver_residual_norm    = - Huge( 1.0_wp )
  End type Ewald_3d_status
  
  Private

  ! Default values for the methods parameters
  Real( wp )          , Parameter :: gauss_tol_default    = 1e-15_wp ! Value for cutting of the gaussian
  Integer             , Parameter :: range_gauss_default  = 12       ! Number of points to grid one side of a gaussian
  Integer             , Parameter :: FD_order_default     = 12       ! Default order of the finite difference approximation
  Character( Len = * ), Parameter :: default_solver       = "minres" ! Default equation solver 
!!$  Character( Len = * ), Parameter :: default_solver       = "CG" 
!!$  Character( Len = * ), Parameter :: default_solver       = "wjac" 
  Real( wp )          , Parameter :: residual_tol_default = 1e-7_wp ! Tolerance on residual in equation solver
  
Contains

  Subroutine Ewald_3d( recipe, q, r, q_halo, r_halo, &
       recip_E, forces, stress, error, q_grid_old, pot_grid_old, ei, q_grid, pot_grid, status ) 

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64
    
    Use symetrically_screened_poisson_module, Only : ssp_long_range

    Implicit None

    Class( Ewald_3d_recipe )           ,              Intent( InOut )           :: recipe ! Need to fix intent eventualy
    Real( wp ), Dimension( 1:         ),              Intent( In    )           :: q
    Real( wp ), Dimension( 1:, 1:     ),              Intent( In    )           :: r
    Real( wp ), Dimension( 1:         ),              Intent( In    )           :: q_halo
    Real( wp ), Dimension( 1:, 1:     ),              Intent( In    )           :: r_halo
    Real( wp )                         ,              Intent(   Out )           :: recip_E
    Real( wp ), Dimension( 1:, 1:     ),              Intent(   Out )           :: forces
    Real( wp ), Dimension( 1:, 1:     ),              Intent(   Out )           :: stress
    Integer                            ,              Intent(   Out )           :: error
    Real( wp ), Dimension(  :,  :,  : ),              Intent( In    ), Optional :: q_grid_old
    Real( wp ), Dimension(  :,  :,  : ),              Intent( In    ), Optional :: pot_grid_old
    Real( wp ), Dimension(  :         ),              Intent(   Out ), Optional :: ei
    Real( wp ), Dimension(  :,  :,  : ), Allocatable, Intent(   Out ), Optional :: q_grid
    Real( wp ), Dimension(  :,  :,  : ), Allocatable, Intent(   Out ), Optional :: pot_grid
    Type( Ewald_3d_status )                         , Intent(   Out ), Optional :: status

    Type( Ewald_3d_status ) :: loc_status

    Real( wp ), Dimension( :, :, : ), Allocatable :: q_grid_loc
    Real( wp ), Dimension( :, :, : ), Allocatable :: pot_grid_loc

    Real( wp ), Dimension( : ), Allocatable :: ei_loc
    
    ! Stress not calculated yet
    stress = 0.0_wp

    ! Allocate the charge and potential grids
    Allocate( q_grid_loc( recipe%domain_base_coords( 1 ):recipe%domain_end_coords( 1 ), &
         recipe%domain_base_coords( 2 ):recipe%domain_end_coords( 2 ), &
         recipe%domain_base_coords( 3 ):recipe%domain_end_coords( 3 ) ) )
    Allocate( pot_grid_loc( recipe%domain_base_coords( 1 ):recipe%domain_end_coords( 1 ), &
         recipe%domain_base_coords( 2 ):recipe%domain_end_coords( 2 ), &
         recipe%domain_base_coords( 3 ):recipe%domain_end_coords( 3 ) ) )

    ! Allocate the per particle energy
    Allocate( ei_loc( 1:Size( q ) ) )

    ! Solve the problem
    Call ssp_long_range( recipe%l, q, r, recipe%alpha, q_halo, r_halo,                                 &
         recipe%range_gauss, recipe%n_grid, Lbound( q_grid_loc ), recipe%residual_tol,                 &
         recip_E, q_grid_loc, pot_grid_loc, recipe%solver, recipe%pot_swapper, recipe%grid_integrator, &
         ei_loc, forces, loc_status%t_grid, loc_status%t_pot_solve, loc_status%t_forces,               &
         loc_status%solver_iterations, loc_status%solver_stop_code, loc_status%solver_stop_message ,   &
         loc_status%solver_residual_norm, error, q_grid_old, pot_grid_old )

    ! Return optional arguments as required
    If( Present( q_grid ) ) Then
       Call Move_alloc( q_grid_loc, q_grid )
    End If

    If( Present( pot_grid ) ) Then
       Call Move_alloc( pot_grid_loc, pot_grid )
    End If

    If( Present( ei ) ) Then
       ei = ei_loc
    End If

    If( Present( status ) ) Then
       status = loc_status
    End If

    If( error > 0 ) Then
       error = EWALD_3D_SOLVER_FAILED
    Else If( error < 0 ) Then
       error = EWALD_3D_ISSUED_WARNING
    Else
       error = EWALD_3D_SUCCESS
    End If

  End Subroutine Ewald_3d

  Subroutine Ewald_3d_init( recipe, l, alpha, error, communicator, equation_solver )

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
    
    Use equation_solver_precon_base_class_module , Only : equation_solver_precon_base_class
    Use equation_solver_minres_module            , Only : equation_solver_minres
    Use equation_solver_conjugate_gradient_module, Only : equation_solver_conjugate_gradient
    Use equation_solver_weighted_jacobi_module   , Only : equation_solver_weighted_jacobi

    Use symetrically_screened_poisson_module, Only : ssp_long_range
    
    ! Should eventually get rid of this one once HALO_HACK is sorted
    ! But will need a more general solution to find base coords of domain
    Use domains_module, Only : domain_get_params

    Implicit None
    

    Class( Ewald_3d_recipe ), Intent(   Out )           :: recipe
    Type( lattice )         , Intent( In    )           :: l
    Real( wp )              , Intent( In    )           :: alpha
    Integer                 , Intent(   Out )           :: error
    Integer                 , Intent( In    ), Optional :: communicator
    Character( Len = * )    , Intent( In    ), Optional :: equation_solver

    Class( comms_base_class                  ), Allocatable :: comms
    Class( FD_template                       ), Allocatable :: FD
    Class( halo_setter_base_class            ), Allocatable :: FD_swapper
    Class( halo_setter_base_class            ), Allocatable :: pot_swapper
    Class( equation_solver_precon_base_class ), Allocatable :: solver
    Class( equation_solver_precon_base_class ), Allocatable :: precon
    Class( quadrature_base_class             ), Allocatable :: grid_integrator

    Real( wp ), Dimension( 1:3, 1:3 ) :: Grid_vecs, dGrid_vecs
    
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
    Integer :: loc_communicator
    Integer :: max_iter
    Integer :: i

    Character( Len = : ), Allocatable :: loc_equation_solver

    ! Set up communication layer. Provide a communicator to indicate the run is MPI parallel
    If( Present( communicator ) ) Then
       Allocate( comms_parallel :: comms )
       loc_communicator = communicator
    Else
       Allocate( comms_serial :: comms )
       ! Any rubbish for the communicator will do if in serial
       loc_communicator = - Huge( loc_communicator )
    End If
    Call comms%set_comm( loc_communicator, error )
    If( error /= 0 ) Then
       error = EWALD_3D_INIT_COMMS_SETUP_FAILED
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
    domain_end_coords = domain_base_coords + n_grid_domain - 1

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
       Call FD_swapper%init ( n_grid_domain, FD%get_order(), loc_communicator, error )
    Class is ( halo_serial_setter )
       Call FD_swapper%init( error )
    End Select
    If( error /= 0 ) Then
       error = EWALD_3D_INIT_FD_SWAPPER_SETUP_FAILED
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
       Call pot_swapper%init ( n_grid_domain, range_gauss + 1, loc_communicator, error )
    Class is ( halo_serial_setter )
       Call pot_swapper%init( error )
    End Select
    If( error /= 0 ) Then
       error = EWALD_3D_INIT_POT_SWAPPER_SETUP_FAILED
       Return
    End If

    ! Set up the equation solver
    ! First decide on the accuracy
    residual_tol = residual_tol_default
    ! TESTING PRECONDITIONER
    Allocate( equation_solver_weighted_jacobi :: precon )
    Call precon%init( max_iter = 3, comms = comms, FD_operator = FD, halo_swapper = FD_swapper  )
    If( Present( equation_solver ) ) Then
       loc_equation_solver = equation_solver 
    Else
       loc_equation_solver = default_solver
    End If
    Select Case( Trim( Adjustl( loc_equation_solver ) ) )
    Case( "cg", "CG" )
       max_iter = 1000
       Allocate( equation_solver_conjugate_gradient :: solver )
    Case( "minres", "MINRES" )
       max_iter = 1000
       Allocate( equation_solver_minres :: solver )
    Case( "WJAC", "wjac" )
       max_iter = 80
       Allocate( equation_solver_weighted_jacobi :: solver )
    Case Default
       error = EWALD_3D_INIT_UNKNOWN_EQUATION_SOLVER
       Return
    End Select
    call solver%init( max_iter = max_iter, comms = comms, FD_operator = FD, halo_swapper = FD_swapper  )
!!$    call solver%init( max_iter = max_iter, comms = comms, FD_operator = FD, halo_swapper = FD_swapper, precon = precon  )


    ! Set up the quadrature method
    Allocate( quadrature_trapezium_rule :: grid_integrator )

    ! Finally set up the recipe
    recipe%l                  = l
    recipe%alpha              = alpha
    recipe%n_grid             = n_grid
    recipe%domain_base_coords = domain_base_coords
    recipe%domain_end_coords  = domain_end_coords
    recipe%range_gauss        = range_gauss
    recipe%gauss_tol          = gauss_tol
    recipe%residual_tol       = residual_tol
    recipe%comms              = comms
    recipe%FD                 = FD
    recipe%FD_swapper         = FD_swapper
    recipe%pot_swapper        = pot_swapper
    recipe%solver             = solver
    recipe%grid_integrator    = grid_integrator

    ! And all is well with the world
    error = EWALD_3D_SUCCESS

  End Subroutine Ewald_3d_init

  Subroutine Ewald_3d_get_ingredients( recipe, l, alpha, n_grid, &
       domain_base_coords, domain_end_coords, range_gauss, gauss_tol, residual_tol, &
       comms, FD, FD_swapper, pot_swapper, solver, grid_integrator )

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use lattice_module                          , Only : lattice
    Use comms_base_class_module                 , Only : comms_base_class
    Use FD_template_module                      , Only : FD_template
    Use halo_setter_base_module                 , Only : halo_setter_base_class
    Use quadrature_base_module                  , Only : quadrature_base_class
    Use equation_solver_precon_base_class_module, Only : equation_solver_precon_base_class

    Implicit None
    
    Class( Ewald_3d_recipe )                               , Intent( In    )           :: recipe
    Type( lattice )                                        , Intent(   Out ), Optional :: l
    Real( wp )                                             , Intent(   Out ), Optional :: alpha
    Integer, Dimension( 1:3 )                              , Intent(   Out ), Optional :: n_grid
    Integer, Dimension( 1:3 )                              , Intent(   Out ), Optional :: domain_base_coords
    Integer, Dimension( 1:3 )                              , Intent(   Out ), Optional :: domain_end_coords
    Integer                                                , Intent(   Out ), Optional :: range_gauss
    Real( wp )                                             , Intent(   Out ), Optional :: gauss_tol
    Real( wp )                                             , Intent(   Out ), Optional :: residual_tol
    Class( comms_base_class                  ), Allocatable, Intent(   Out ), Optional :: comms
    Class( FD_template                       ), Allocatable, Intent(   Out ), Optional :: FD
    Class( halo_setter_base_class            ), Allocatable, Intent(   Out ), Optional :: FD_swapper
    Class( halo_setter_base_class            ), Allocatable, Intent(   Out ), Optional :: pot_swapper
    Class( equation_solver_precon_base_class ), Allocatable, Intent(   Out ), Optional :: solver
    Class( quadrature_base_class             ), Allocatable, Intent(   Out ), Optional :: grid_integrator

    If( Present( l ) ) Then
       l = recipe%l
    End If

    If( Present( alpha ) ) Then
       alpha = recipe%alpha
    End If

    If( Present( n_grid ) ) Then
       n_grid = recipe%n_grid
    End If
       
    If( Present( domain_base_coords ) ) Then
       domain_base_coords = recipe%domain_base_coords
    End If

    If( Present( domain_end_coords ) ) Then
       domain_end_coords = recipe%domain_end_coords
    End If

    If( Present( range_gauss ) ) Then
       range_gauss = recipe%range_gauss
    End If
       
    If( Present( gauss_tol ) ) Then
       gauss_tol = recipe%gauss_tol
    End If

    If( Present( residual_tol ) ) Then
       residual_tol = recipe%residual_tol
    End If

    If( Present( comms ) ) Then
       comms = recipe%comms
    End If

    If( Present( FD ) ) Then
       FD = recipe%FD
    End If

    If( Present( FD_swapper ) ) Then
       FD_swapper = recipe%FD_swapper
    End If

    If( Present( pot_swapper ) ) Then
       pot_swapper = recipe%pot_swapper
    End If

    If( Present( solver ) ) Then
       solver = recipe%solver
    End If

    If( Present( grid_integrator ) ) Then
       grid_integrator = recipe%grid_integrator
    End If

  End Subroutine Ewald_3d_get_ingredients
  
End Module Ewald_3d_module
