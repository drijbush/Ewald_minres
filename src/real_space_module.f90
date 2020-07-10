Module real_space_module

  Implicit None

  Public :: real_space_energy

  Private

Contains

  Subroutine real_space_energy( l, q, r, alpha, n_G_shells, E )

    Use constants,      Only : wp
    Use lattice_module, Only : lattice

    Type( lattice ),               Intent( In    ) :: l
    Real( wp ), Dimension( :    ), Intent( In    ) :: q
    Real( wp ), Dimension( :, : ), Intent( In    ) :: r
    Real( wp ),                    Intent( In    ) :: alpha
    Integer,                       Intent( In    ) :: n_G_shells
    Real( wp ),                    Intent(   Out ) :: E

    Real( wp ), Dimension( 1:3 ) :: G
    Real( wp ), Dimension( 1:3 ) :: ri, rjG
    Real( wp ), Dimension( 1:3 ) :: rj, rij

    Real( wp ) :: rijG, rijG_sq
    Real( wp ) :: c, y, t

    Integer :: n
    Integer :: side, max_G_vector
    Integer :: i, j, iG

    n = Size( q )

    ! Implicit barrier at end of single makes sure shared E is inited by all threds
    !$omp single
    E = 0.0_wp
    !$omp end single

    Select Case( n_G_shells )

    Case( 1: )
       ! Include lattice sum over the number of shells of G vectors
       side =  2 * n_G_shells + 1
       max_G_vector = n_G_shells * side * side + ( ( side + 1 ) * side ) / 2 - 1
       c = 0.0_wp
       !$omp do reduction( +:E )
       Do i = 1, n
          ri = r( :, i )
          Do j = 1, n
             If( i /= j ) Then
                ! Zero G vector
                rjG = r( :, j )
                rijG_sq = Dot_product( ri - rjG, ri - rjG )
                rijG = Sqrt( rijG_sq )
                y = q( i ) * q( j ) * Erfc( alpha * rijG ) / rijG - c
                t = E + y
                c = ( t - E ) - y
                E = t
             End If
             Do iG = 1, max_G_vector
                !+ve G vector
                Call l%get_nth_dir_vec( ig + 1, G )
                rjG = r( :, j ) + G
                rijG_sq = Dot_product( ri - rjG, ri - rjG )
                rijG = Sqrt( rijG_sq )
                y = q( i ) * q( j ) * Erfc( alpha * rijG ) / rijG - c
                t = E + y
                c = ( t - E ) - y
                E = t
                !-ve G vector
                rjG = r( :, j ) - G
                rijG_sq = Dot_product( ri - rjG, ri - rjG )
                rijG = Sqrt( rijG_sq )
                y = q( i ) * q( j ) * Erfc( alpha * rijG ) / rijG - c
                t = E + y
                c = ( t - E ) - y
                E = t
             End Do
          End Do
       End Do
       !$omp end do

    Case( 0 )
       ! Just minimum image convention
       c = 0.0_wp
       !$omp do reduction( +:E )
       Do i = 1, n
          ri = r( :, i )
          Do j = 1, n
             If( i /= j ) Then
                rj = r( :, j )
                ! Apply minimum image convention
                Call l%minimum_image( ri - rj, rij )
                rijG_sq = Dot_product( rij, rij )
                rijG = Sqrt( rijG_sq )
                ! Should implement rcut at some point
!!$                If( rijG <= 11.0_wp ) Then
                y = q( i ) * q( j ) * Erfc( alpha * rijG ) / rijG - c
                t = E + y
                c = ( t - E ) - y
                E = t
!!$                End If
             End If
          End Do
       End Do
       !$omp end do

    Case Default
       Stop "Negative number of shells in real space lattice summation"

    End Select

    ! Make sure scaling is just applied by one thread
    !$omp single
    E = E * 0.5_wp
    !$omp end single

  End Subroutine real_space_energy

End Module real_space_module
