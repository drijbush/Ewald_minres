Module equation_solver_precon_base_class_module

  Use equation_solver_base_class_module, Only : equation_solver_base_class

  Type, Public, Abstract, Extends( equation_solver_base_class ) :: equation_solver_precon_base_class
     Class( equation_solver_base_class ), Allocatable :: precon
   Contains
     Generic  , Public  :: init => precon_init
     Procedure, Private :: precon_init
  End type equation_solver_precon_base_class

Contains

  Subroutine precon_init( method, max_iter, comms, FD_operator, halo_swapper, precon )

    Use comms_base_class_module          , Only : comms_base_class
    Use halo_setter_base_module          , Only : halo_setter_base_class
    Use FD_template_module               , Only : FD_template
    
    Implicit None

    Class( equation_solver_precon_base_class ), Intent( InOut )           :: method
    Integer                                   , Intent( In    ), Optional :: max_iter
    Class( comms_base_class       )           , Intent( In    ), Optional :: comms
    Class( FD_template )                      , Intent( In    ), Optional :: FD_operator
    Class( halo_setter_base_class )           , Intent( In    ), Optional :: halo_swapper
    Class( equation_solver_base_class )       , Intent( In    )           :: precon

    method%precon = precon

    Call method%init( max_iter, comms, FD_operator, halo_swapper )

  End Subroutine precon_init

  
End Module equation_solver_precon_base_class_module
