Module grid_vectors_module
  !!----------------------------------------------------
  !! Module containing type and functions pertaining to
  !! calculating finite-difference grid contributions
  !! for non-orthonormal grids.
  !!
  !! Written by I.J. Bush
  !!----------------------------------------------------



  Use numbers_module, Only : wp

  Implicit None

  Type, Public :: grid_vectors
    !!----------------------------------------------------
    !! Type to handle grids of arbitrary size and orientation
    !! and the calculation of relevant parameters thereof
    !!
    !! Written by I.J. Bush
    !!----------------------------------------------------
    Real( wp ), Dimension( :, : ), Allocatable, Private :: dir_vecs
    Real( wp ),                                 Private :: volume
    Real( wp ), Dimension( :, : ), Allocatable, Private :: inv_vecs
  Contains
    Procedure, Public :: set_dir_vecs
    Procedure, Public :: get_dir_vec
    Procedure, Public :: get_volume
    Procedure, Public :: get_inv_vec
  End type grid_vectors

  Private

Contains

  Subroutine set_dir_vecs( g, vecs )
    !!----------------------------------------------------
    !! Set grid vectors as an NxN array (vecs) and thereby
    !! calculate the various properties for use with the
    !! rest of the solvers
    !!
    !! Written by I.J. Bush
    !!----------------------------------------------------

    Class( grid_vectors )        , Intent( InOut ) :: g
    Real( wp ), Dimension( :, : ), Intent( In    ) :: vecs

    Real( wp ) :: V

    g%dir_vecs = vecs

    Call givens_invert( g%dir_vecs, g%inv_vecs, V )
    g%volume = Abs( V )

  End Subroutine set_dir_vecs

  Pure Function get_dir_vec( g, which ) Result( v )
    !!----------------------------------------------------
    !! Get the direction of a given (which) vector from
    !! the direction vecs
    !!
    !! Written by I.J. Bush
    !!----------------------------------------------------

    Real( wp ), Dimension( : ), Allocatable :: v

    Class( grid_vectors ), Intent( In ) :: g
    Integer              , Intent( In ) :: which

    v = g%dir_vecs( :, which )

  End Function get_dir_vec

  Pure Function get_inv_vec( g, which ) Result( v )
    !!----------------------------------------------------
    !! Get the inverse direction vectors
    !!
    !! Written by I.J. Bush
    !!----------------------------------------------------

    Real( wp ), Dimension( : ), Allocatable :: v

    Class( grid_vectors ), Intent( In ) :: g
    Integer              , Intent( In ) :: which

    v = g%inv_vecs( which, : )

  End Function get_inv_vec

  Pure Function get_volume( g ) Result( V )
    !!----------------------------------------------------
    !! Get the volume element from the vectors
    !!
    !! Written by I.J. Bush
    !!----------------------------------------------------

    Real( wp ) :: V

    Class( grid_vectors ), Intent( In ) :: g

    V = g%volume

  End Function get_volume

  Subroutine givens_invert( A, B, det )
    !!----------------------------------------------------
    !! Invert the matrix A, returning the inverse in B
    !! and the determinant of A in det
    !!
    !! NOT FOR USE ON LARGE MATRICES
    !!
    !! Really designed as a robust, not totally disastorous, way
    !! to find inverse lattice vectors and the volume of the unit cell
    !! ( i.e. the determinant )
    !!
    !! The method used is to QR factorise ( in fact LQ factorise ),
    !! invert the triangular matrix and the from that form the
    !! inverse of the original matrix. The determinant is simply the
    !! product of the diagonal elements of the triangular matrix.
    !! Method chosen as easy to implement and nicely numerically stable,
    !! especially for the small matrices need here.
    !!
    !! Written by I.J. Bush
    !!----------------------------------------------------

    Real( wp ), Dimension( 1:, 1: ),               Intent( In    )           :: A
    Real( wp ), Dimension(  :,  : ), Allocatable,  Intent(   Out )           :: B
    Real( wp )                     ,               Intent(   Out ), Optional :: det

    Real( wp ), Dimension( :, : ), Allocatable :: Q, L

    Real( wp ) :: theta, c, s
    Real( wp ) :: ann, anm
    Real( wp ) :: Lim, Lin
    Real( wp ) :: Qmi, Qni
    Real( wp ) :: mv

    Integer :: dim
    Integer :: m, n
    Integer :: i

    Allocate( B, Mold = A )

    ! Size of the problem
    dim = Size( A, dim = 1 )

    ! The factors of the matrix A
    Allocate( L( 1:dim, 1:dim ) )
    Allocate( Q( 1:dim, 1:dim ) )

    ! Factorise A by Givens rotations
    L = A
    Q = 0.0_wp
    Do i = 1, dim
      Q( i, i ) = 1.0_wp
    End Do
    LQ_factorise: Do m = 2, dim
      Do n = 1, m - 1
        ! For element mn generate the rotation which
        ! will zero that element
        ann = L( n, n )
        anm = L( n, m )
        theta = atan2( - anm, ann )
        c = Cos( theta )
        s = Sin( theta )
        ! And apply that rotation to the original matrix
        Do i = 1, Dim
          Lim = L( i, m )
          Lin = L( i, n )
          L( i, n ) = c * Lin - s * Lim
          L( i, m ) = c * Lim + s * Lin
        End Do
        ! And then apply the rotation to the orthognal matrix
        Do i = 1, Dim
          Qmi = Q( m, i )
          Qni = Q( n, i )
          Q( n, i ) = c * Qni - s * Qmi
          Q( m, i ) = c * Qmi + s * Qni
        End Do
      End Do
    End Do LQ_factorise

    ! Now invert the triangular matrix, returning the inverse in B
    B(  1, 1  ) = 1.0_wp / L( 1, 1 )
    B(  1, 2: ) = 0.0_wp
    Do m = 2, dim
      Do n = 1, m
        mv = 0.0_wp
        Do i = 1, m - 1
          mv = mv - L( m, i ) * B( i, n )
        End Do
        If( m == n ) Then
          mv = mv + 1.0_wp
        End If
        mv = mv / L( m, m )
        B( m, n ) = mv
      End Do
      B( m, n: ) = 0.0_wp
    End Do

    ! The inverse of A is now trivially formed
    B = Matmul( Transpose( Q ), B )

    ! And generate the determinant
    If( Present( det ) ) Then
      det = 1.0_wp
      Do i = 1, dim
        det = det * L( i, i )
      End Do
    End If

    Deallocate( Q )
    Deallocate( L )

  End Subroutine givens_invert

End Module grid_vectors_module
