Module domains_module

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
       If( domain_is_in_grid_volume( l, r( :, i ), n_grid, domain_base_coords, domain_base_coords + n_grid_domain ) ) Then
          i_domain = [ i_domain, i ]
       End If
    End Do

    q_domain = q(    i_domain )
    r_domain = r( :, i_domain )
    
  End Subroutine domain_build

  Pure Subroutine domain_halo_build( l, q, r, n_grid, n_proc, domain_coords, halo_width, q_domain, r_domain )

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
    Real( wp ), Dimension(      : ), Allocatable, Intent(   Out ) :: q_domain
    Real( wp ), Dimension(  :,  : ), Allocatable, Intent(   Out ) :: r_domain

    Integer, Dimension( : ), Allocatable :: i_halo

    Integer, Dimension( 1:3 ) :: domain_base_coords
    Integer, Dimension( 1:3 ) :: n_grid_domain
    
    Integer :: n
    Integer :: i
    
    n = Size( q )

    Call domain_get_params( n_grid, n_proc, domain_coords, n_grid_domain, domain_base_coords )
    
    Allocate( i_halo( 1:0 ) )

    Do i = 1, n
       ! Is it in the volume including the domain and the halo
       If( domain_is_in_grid_volume( l, r( :, i ), n_grid, &
            domain_base_coords - halo_width, domain_base_coords + n_grid_domain + halo_width ) ) Then
          ! If it is and it is not in the domain it is in the halo
          If( .Not. domain_is_in_grid_volume( l, r( :, i ), n_grid, &
               domain_base_coords, domain_base_coords + n_grid_domain ) ) Then
             i_halo = [ i_halo, i ]
          End If
       End If
    End Do

    q_domain = q(    i_halo )
    r_domain = r( :, i_halo )
    
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
