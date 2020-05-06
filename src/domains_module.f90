Module domains_module

  ! NEED to think how better to map domains onto DL_POLY domains
  ! Here we map purely using the grid
  ! DLP (effectively) maps using the fractional coordinates
  ! Use this for debugging for the moment
  
  Implicit None

  Public :: domain_get_params
  Public :: domain_is_in_grid_volume
  Public :: domain_build
  Public :: domain_halo_build
  
  Private

Contains  

  Pure Subroutine domain_build( l, q, r, n_grid, n_proc, domain_coords, q_domain, r_domain )

    ! GRIDS START AT ZERO

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use lattice_module, Only : lattice

    Implicit none

    Type( lattice )                ,              Intent( In    ) :: l
    Real( wp ), Dimension( 1:     ),              Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1: ),              Intent( In    ) :: r
    Integer   , Dimension( 1:3    ),              Intent( In    ) :: n_grid
    Integer   , Dimension( 1:3    ),              Intent( In    ) :: n_proc
    Integer   , Dimension( 1:3    ),              Intent( In    ) :: domain_coords
    Real( wp ), Dimension(      : ), Allocatable, Intent(   Out ) :: q_domain
    Real( wp ), Dimension(  :,  : ), Allocatable, Intent(   Out ) :: r_domain

    Integer, Dimension( : ), Allocatable :: i_domain

    Integer, Dimension( 1:3 ) :: domain_base_coords
    Integer, Dimension( 1:3 ) :: n_grid_domain
    
    Integer :: n
    Integer :: i
    
    n = Size( q )

    Call domain_get_params( n_grid, n_proc, domain_coords, n_grid_domain, domain_base_coords )
    
    Allocate( i_domain( 1:0 ) )

    Do i = 1, n
       If( domain_is_in_grid_volume( l, r( :, i ), n_grid, domain_base_coords, &
            domain_base_coords + n_grid_domain ) ) Then
          i_domain = [ i_domain, i ]
       End If
    End Do

    q_domain = q(    i_domain )
    r_domain = r( :, i_domain )
    
  End Subroutine domain_build

  Pure Subroutine domain_halo_build( l, q, r, n_grid, n_proc, domain_coords, halo_width, q_halo, r_halo )

    ! GRIDS START AT ZERO

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use lattice_module, Only : lattice

    Implicit none

    Type( lattice )                ,              Intent( In    ) :: l
    Real( wp ), Dimension( 1:     ),              Intent( In    ) :: q
    Real( wp ), Dimension( 1:, 1: ),              Intent( In    ) :: r
    Integer   , Dimension( 1:3    ),              Intent( In    ) :: n_grid
    Integer   , Dimension( 1:3    ),              Intent( In    ) :: n_proc
    Integer   , Dimension( 1:3    ),              Intent( In    ) :: domain_coords
    Integer   , Dimension( 1:3    ),              Intent( In    ) :: halo_width
    Real( wp ), Dimension(      : ), Allocatable, Intent(   Out ) :: q_halo
    Real( wp ), Dimension(  :,  : ), Allocatable, Intent(   Out ) :: r_halo

    Real( wp ), Dimension( :, : ), Allocatable :: rtmp
    
    Real( wp ), Dimension( 1:3 ) :: G, ri, riG
    
    Integer, Dimension( 1:3 ) :: domain_base_coords
    Integer, Dimension( 1:3 ) :: n_grid_domain
    
    Integer :: n
    Integer :: i
    Integer :: iGx, iGy, iGz
    
    n = Size( q )

    Call domain_get_params( n_grid, n_proc, domain_coords, n_grid_domain, domain_base_coords )

    Allocate( q_halo( 1:0 ) )
    Allocate( r_halo( 1:3, 1:0 ) )
    
    Do i = 1, n
       ! For the halo we have to consider periodic images
       ri = r( :, i )
       Do iGz = -1, 1
          Do iGy = -1, 1
             Do iGx = -1, 1
                Call l%get_dir_vec( [ iGz, iGy, iGz ], G )
                riG = ri + G
                ! HACK assume we only need look at most one lattice vector away
                ! Is it in the volume including the domain and the halo
                If( domain_is_in_grid_volume( l, riG, n_grid, &
                     domain_base_coords - halo_width, domain_base_coords + n_grid_domain + halo_width ) ) Then
                   ! If it is and it is not in the domain it is in the halo
                   If( .Not. domain_is_in_grid_volume( l, riG, n_grid, &
                        domain_base_coords, domain_base_coords + n_grid_domain ) ) Then
                      ! Need to store shifted position of atom
                      q_halo = [ q_halo, q( i ) ]
                      Allocate( rtmp( 1:3, 1:Size( r_halo, Dim = 2 ) + 1 ) )
                      rtmp( :, 1:Size( r_halo, Dim = 2 ) ) = r_halo
                      rtmp( :, Size( r_halo, Dim = 2 ) + 1 ) = riG
                      Call move_alloc( rtmp, r_halo )
                   End If
                End If
             End Do
          End Do
       End Do
    End Do
    
  End Subroutine domain_halo_build

  Pure Function domain_is_in_grid_volume( l, ri, n_grid, lo, hi ) Result( in_volume )

    ! GRIDS START AT ZERO

    Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

    Use lattice_module, Only : lattice

    Logical :: in_volume

    Type( lattice ),                   Intent( In ) :: l
    Real( wp )     , Dimension( 1:3 ), Intent( In ) :: ri
    Integer        , Dimension( 1:3 ), Intent( In ) :: n_grid
    Integer        , Dimension( 1:3 ), Intent( In ) :: lo
    Integer        , Dimension( 1:3 ), Intent( In ) :: hi

    Real( wp ), Dimension( 1:3 ) :: fi

    Integer, Dimension( 1:3 ) :: i_grid

    Call l%to_fractional( ri, fi )
    i_grid = Int( fi * n_grid )
    in_volume = All( i_grid >= lo ) .And. All( i_grid < hi )
    
  End Function domain_is_in_grid_volume

  Elemental Subroutine domain_get_params( n_grid, n_proc, domain_coords, n_grid_domain, domain_base_coords )

    ! GRIDS START AT ZERO

    Integer, Intent( In    ) :: n_grid
    Integer, Intent( In    ) :: n_proc
    Integer, Intent( In    ) :: domain_coords
    Integer, Intent(   Out ) :: n_grid_domain
    Integer, Intent(   Out ) :: domain_base_coords

    Integer :: n_per_proc, n_left, n_me, n_first

    n_per_proc = n_grid / n_proc
    n_left     = n_grid - n_proc * n_per_proc

    If( domain_coords < n_left ) Then
       n_me = n_per_proc + 1
    Else
       n_me = n_per_proc
    End If

    If( domain_coords < n_left ) Then
       n_first = ( n_per_proc + 1 ) * domain_coords
    Else
       n_first = ( n_per_proc + 1 ) * n_left + ( domain_coords - n_left ) * n_per_proc
    End If

    n_grid_domain      = n_me
    domain_base_coords = n_first
    
  End Subroutine domain_get_params
  
End Module domains_module
