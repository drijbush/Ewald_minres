Module FD_template_module
  !!-----------------------------------------------------------
  !! Module containing type to construct finite difference weights
  !! for arbitrary order derivatives, stencils and arbitrary grids
  !!
  !! Written by I.J. Bush
  !!-----------------------------------------------------------

  Use numbers_module     , Only : wp
  Use grid_vectors_module, Only : grid_vectors

  Implicit None

  Private

  Type, Abstract, Extends( grid_vectors ), Public :: FD_template
     !!-----------------------------------------------------------
     !! Low-level type for handling weights for arbitray stencils
     !! and derivatives, contains grid_vector information through
     !! the extension of the grid_vectors type.
     !!
     !! Do not use this type directly.
     !!
     !! Written by I.J. Bush
     !!-----------------------------------------------------------
     Integer                                   , Private :: max_deriv
     Integer                                   , Private :: order
     Real( wp ), Dimension( :, : ), Allocatable, Private :: weights
     Real( wp )                                , Public  :: diag_inv ! HACK AS PUBLIC - NEED TO FIX EVENTUALLY
   Contains
     Procedure                    , Public           :: FD_init
     Procedure                    , Public           :: set_order
     Procedure                    , Public           :: get_order
     Procedure                    , Public           :: get_max_deriv
     Procedure                    , Public           :: get_weight
     Procedure                    , Public           :: get_diag_inv
     Procedure(  apply_interface ), Public, Deferred :: apply
     Procedure( jacobi_interface ), Public, Deferred :: jacobi_sweep
  End type FD_template

  Interface

     Subroutine apply_interface( FD, grid_lb, lap_lb, start, final, grid, laplacian )
       Import :: wp, FD_template
       Implicit None
       Class( FD_template )                                          , Intent( In    ) :: FD
       Integer, Dimension( 1:3 )                                     , Intent( In    ) :: grid_lb !! lower bounds grid
       Integer, Dimension( 1:3 )                                     , Intent( In    ) :: lap_lb  !! lower bounds laplacian
       Integer, Dimension( 1:3 )                                     , Intent( In    ) :: start   !! start point for calculation
       Integer, Dimension( 1:3 )                                     , Intent( In    ) :: final   !! final point for calculation
       Real( wp ), Dimension( grid_lb(1):, grid_lb(2):, grid_lb(3): ), Intent( In    ) :: grid    !! Thing to be differentiated
       Real( wp ), Dimension(  lap_lb(1):,  lap_lb(2):,  lap_lb(3): ), Intent(   Out ) :: laplacian
     End Subroutine apply_interface

     Subroutine jacobi_interface( FD, grid_lb, lap_lb, start, final, jac_weight, grid, soln_in, soln_out )
       Import :: wp, FD_template
       Class( FD_template )                                          , Intent( In    ) :: FD
       Integer, Dimension( 1:3 )                                     , Intent( In    ) :: grid_lb    !! lower bounds grid
       Integer, Dimension( 1:3 )                                     , Intent( In    ) :: lap_lb     !! lower bounds laplacian
       Integer, Dimension( 1:3 )                                     , Intent( In    ) :: start      !! start point for calculation
       Integer, Dimension( 1:3 )                                     , Intent( In    ) :: final      !! final point for calculation
       Real( wp )                                                    , Intent( In    ) :: jac_weight !! The weight for the weighted jacobi
       Real( wp ), Dimension( grid_lb(1):, grid_lb(2):, grid_lb(3): ), Intent( In    ) :: grid       !! Thing to be differentiated
       Real( wp ), Dimension(  lap_lb(1):,  lap_lb(2):,  lap_lb(3): ), Intent( In    ) :: soln_in    !! Solution at start of sweep 
       Real( wp ), Dimension( grid_lb(1):, grid_lb(2):, grid_lb(3): ), Intent(   Out ) :: soln_out   !! Solution at end of sweep

     End Subroutine jacobi_interface

  End Interface

Contains

  Subroutine FD_init( FD, max_deriv, order, vecs )
    !!-----------------------------------------------------------
    !! Initialise finite difference object for arbitrary grids
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------

    Class( FD_template )         , Intent(   Out ) :: FD
    Integer                      , Intent( In    ) :: max_deriv
    Integer                      , Intent( In    ) :: order
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs

    Call FD%set_order( max_deriv, order )
    Call FD%set_dir_vecs( vecs )

  End Subroutine FD_init

  Subroutine set_order( FD, max_deriv, order )
    !!-----------------------------------------------------------
    !! Set derivative and stencil and calculate respective weights
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------

    Class( FD_template ), Intent( InOut ) :: FD
    Integer             , Intent( In    ) :: max_deriv
    Integer             , Intent( In    ) :: order

    Integer :: i

    FD%max_deriv = max_deriv
    FD%order     = order

    Call weights( 0.0_wp, [ ( Real( i, wp ), i = -FD%order, FD%order ) ], &
         FD%max_deriv, FD%weights )

  End Subroutine set_order

  Pure Function get_order( FD ) Result( order )
    !!-----------------------------------------------------------
    !! Return stencil order
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------

    Integer :: order

    Class( FD_template ), Intent( In ) :: FD

    order = FD%order

  End Function get_order

  Pure Function get_max_deriv( FD ) Result( max_deriv )
    !!-----------------------------------------------------------
    !! Return maximum possible derivative order
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------

    Integer :: max_deriv

    Class( FD_template ), Intent( In ) :: FD

    max_deriv = FD%max_deriv

  End Function get_max_deriv

  Pure Function get_weight( FD, deriv  ) Result( weight )
    !!-----------------------------------------------------------
    !! Return FD weights for given derivative order
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------
    Real( wp ), Allocatable, Dimension( : ) :: weight

    Class( FD_template ), Intent( In ) :: FD
    Integer             , Intent( In ) :: deriv

    Allocate( weight( -FD%order:FD%order ) )
    weight = FD%weights( :, deriv )

  End Function get_weight

  Pure Function get_diag_inv( FD ) Result( diag_inv )
    !!-----------------------------------------------------------
    !! Return maximum possible derivative order
    !!
    !! Written by I.J. Bush
    !!-----------------------------------------------------------

    Real( wp ) :: diag_inv

    Class( FD_template ), Intent( In ) :: FD

    diag_inv = FD%diag_inv

  End Function get_diag_inv

  Subroutine weights (centre,sample_points,deriv,coeffs)
    !!-----------------------------------------------------------
    !! Input Arguments
    !! centre                      -- location where approximations are to be accurate,
    !! sample_points(0:n_samples)  -- grid point locations
    !! deriv                       -- highest derivative for which weights are sought,
    !!
    !! Output Arguments
    !! coeffs(0:n_samples,0:deriv) -- weights at sample_points for derivatives order deriv
    !!
    !! Taken from https://pdfs.semanticscholar.org/8bf5/912bde884f6bd4cfb4991ba3d077cace94c0.pdf
    !! By Bengt Fornberg
    !!-----------------------------------------------------------

    Implicit None
    Real( wp ), Intent( In ) :: centre
    Real( wp ), Dimension( 0: ), Intent( In ) :: sample_points
    Integer, Intent( In ) :: deriv
    Real( wp ), Dimension( :, : ), Allocatable, Intent( Out ) :: coeffs

    Integer :: n
    Real( wp ) :: c1, c2, c3, c4, c5
    Integer :: mn
    Integer :: i, j, k
    Integer :: ierr

    n = size(sample_points) - 1

    Allocate( coeffs( 0:n, 0:deriv ), stat = ierr )
    If ( ierr /= 0 ) Stop "Error in allocation of coeffs in weights"

    c1 = 1.0_wp
    c4 = sample_points(0) - centre
    coeffs = 0.0_wp
    coeffs(0,0) = 1.0_wp

    Do i=1,n

       mn = Min(i,deriv)
       c2 = 1.0_wp
       c5 = c4
       c4 = sample_points(i) - centre

       Do j=0,i-1

          c3 = sample_points(i)-sample_points(j)
          c2 = c2*c3

          If (j.Eq.i-1) Then ! If final iteration

             Do k=mn,1,-1
                coeffs(i,k) = c1*(k*coeffs(i-1,k-1)-c5*coeffs(i-1,k))/c2
             End Do

             coeffs(i,0) = -c1*c5*coeffs(i-1,0)/c2
          Endif

          Do k=mn,1,-1
             coeffs(j,k) = (c4*coeffs(j,k)-k*coeffs(j,k-1))/c3
          End Do

          coeffs(j,0) = c4*coeffs(j,0)/c3
       End Do

       c1 = c2
    End Do

  End Subroutine weights

End Module FD_template_module
