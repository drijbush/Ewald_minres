Module lattice_module
  
  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Type, Public :: lattice
     Private
     Integer   ,                                     Private :: dim           !! Dimensionality
     Real( wp ),                                     Private :: V             !! Volume
     Real( wp ), Dimension( :, :     ), Allocatable, Private :: given_vecs    !! Vectors initally given
     Real( wp ), Dimension( 1:3, 1:3 ),              Private :: given_vecs_3d !! Vectors Used
     Real( wp ), Dimension( 1:3, 1:3 ),              Private :: rotation      !! Rotation given -> used
     Real( wp ), Dimension( 1:3, 1:3 ),              Private :: dir_vecs      !! Direct vectors Used
     Real( wp ), Dimension( 1:3, 1:3 ),              Private :: rec_vecs      !! Recip vectors used
     Integer   , Dimension(  :,   :  ), Allocatable, Private :: dir_vec_list
     Integer   , Dimension(  :,   :  ), Allocatable, Private :: rec_vec_list
   Contains
     Procedure :: initialise
     Procedure :: print
     Procedure :: get_volume
     Procedure :: get_direct_vectors
     Procedure :: get_reciprocal_vectors
     Procedure :: get_n_dir_vecs
     Procedure :: get_n_rec_vecs
     Procedure :: get_nth_dir_vec
     Procedure :: get_nth_rec_vec
     Procedure :: to_reference
     Procedure :: minimum_image
     Procedure :: to_direct
     Procedure :: to_fractional
     Procedure :: get_dir_vec
     Procedure :: get_rec_vec
     Procedure :: rotate_given_to_used
  End type lattice
  
  Private

  Real( wp ), Parameter :: pi = 3.141592653589793238462643383279502884197_wp
  
Contains

  Pure Subroutine initialise( l, nd, vecs, alpha )

    Use sort_module          , Only : sort_index, SORT_ASCEND
    Use vec3_utilities_module, Only : Operator( .vec. ), vec3_invert
    
    Class( lattice ),                    Intent( InOut ) :: l
    Integer         ,                    Intent( In    ) :: nd
    Real( wp )      , Dimension( :, : ), Intent( In    ) :: vecs
    Real( wp )      ,                    Intent( In    ) :: alpha

    Real( wp ), Parameter :: default_dir_vec_cutoff = 500.0_wp
!!$    Real( wp ), Parameter :: default_dir_vec_cutoff = 5000.0_wp
    Real( wp ), Parameter :: default_rec_vec_cutoff = 100.0_wp

    Real( wp ), Dimension( 1:3, 1:3 ) :: Q12, Q23, Q

    Real( wp ), Dimension( 1:3 ) :: a, a23

    Real( wp ) :: s, c, theta
    Real( wp ) :: l1, l2
    Real( wp ) :: vec_cutoff

    l%dim = nd
    l%given_vecs = vecs( 1:3, 1:nd )

    ! Work out the "volume" for each dimensionality
    Select Case( l%dim )
    Case( 0 )
       l%V = Huge( l%V )
    Case( 1 )
       l%V = Sqrt( Dot_product( l%given_vecs( :, 1 ), l%given_vecs( :, 1 ) ) )
    Case( 2 )
       a = l%given_vecs( :, 1 ) .vec. l%given_vecs( :, 2 )
       l%V = Sqrt( Dot_product( a, a ) )
    Case( 3 )
       l%V = Abs( Dot_product( l%given_vecs( :, 1 ), l%given_vecs( :, 2 ) .vec. l%given_vecs( :, 3 ) ) )
    End Select

    ! For cases where the system is not periodic in all 3 dimensions
    ! Generate vectors orthogonal to the given vectors for the other (non-periodic) directions
    ! and make them very long
    l%given_vecs_3d( :, 1:nd ) = l%given_vecs
    Select Case( l%dim )

    Case( 0 )
       l%given_vecs_3d = 0.0_wp
       l%given_vecs_3d( 1, 1 ) = 1.0_wp
       l%given_vecs_3d( 2, 2 ) = 1.0_wp
       l%given_vecs_3d( 3, 3 ) = 1.0_wp

    Case( 1 )
       ! Generate two vectors orthogonal to the original one
       ! Do this by working out what rotation takes the original vector to the x axis,
       ! and then apply that to the y and z vectors to produce the orthogonal vectors
       l1 = Sqrt( Dot_product( l%given_vecs_3d( :, 1 ), l%given_vecs_3d( :, 1 ) ) )

       ! First rotate the z component to zero
       theta = Atan2( - l%given_vecs_3d( 3, 1 ), l%given_vecs_3d( 2, 1 ) )
       s = Sin( theta )
       c = Cos( theta )
       Q23 = 0.0_wp
       Q23( 1, 1 ) = 1.0_wp
       Q23( 2, 2 ) =   c
       Q23( 3, 2 ) = + s
       Q23( 2, 3 ) = - s
       Q23( 3, 3 ) =   c
       a23 = Matmul( Q23, l%given_vecs_3d( :, 1 ) )

       ! Now rotate the y component to zero
       theta = Atan2( - a23( 2 ), a23( 1 ) )
       s = Sin( theta )
       c = Cos( theta )
       Q12 = 0.0_wp
       Q12( 3, 3 ) = 1.0_wp
       Q12( 1, 1 ) =   c
       Q12( 2, 1 ) = + s
       Q12( 1, 2 ) = - s
       Q12( 2, 2 ) =   c

       ! Generate the total rotation
       Q = Matmul( Q12, Q23 )

       ! Rotate the y and z axes to produce the new vectors
       l%given_vecs_3d( :, 2 ) = Matmul( Transpose( Q ), [  0.0_wp, l1    , 0.0_wp ] )
       l%given_vecs_3d( :, 3 ) = Matmul( Transpose( Q ), [  0.0_wp, 0.0_wp, l1     ] )
       
    Case( 2 )
       ! Generate a vector orthogonal to the two exisiting vectors
       l1 = Sqrt( Dot_product( l%given_vecs_3d( :, 1 ), l%given_vecs_3d( :, 1 ) ) )
       l2 = Sqrt( Dot_product( l%given_vecs_3d( :, 2 ), l%given_vecs_3d( :, 2 ) ) )
       l%given_vecs_3d( :, 3 ) = l%given_vecs_3d( :, 1 ) .vec. l%given_vecs_3d( :, 2 ) / ( Sqrt( l1 * l2 ) )

    Case( 3 )
       ! Nothing to do in 3d

    End Select

    ! This allows a possible rotation to a standardised orientation. Don't use it for now
    ! i.e. set to unit matrix
    l%rotation = 0.0_wp
    l%rotation( 1, 1 ) = 1.0_wp
    l%rotation( 2, 2 ) = 1.0_wp
    l%rotation( 3, 3 ) = 1.0_wp

    ! And the vectors we shall use
    l%dir_vecs = Matmul( l%rotation, l%given_vecs_3d )

    ! Now generate a list of vectors in a sphere of the given size, and then sort them by increasing length
    ! Only use the "z" is positive part of space as the nexgative part can be accessed by transpose symmetry
    vec_cutoff = default_dir_vec_cutoff
    Call generate_vec_list( l%dim, vec_cutoff, l%dir_vecs, l%dir_vec_list )

    ! Generate the reciprocal of the vectors we are using
    Call vec3_invert( l%dir_vecs, l%rec_vecs )
    l%rec_vecs = Transpose( l%rec_vecs )

    ! For large cells recip vectors take up too much memory
    ! Generate a list of reciprocal ectors in order of increasin length
    ! Only use the "z" is positive part of space as the nexgative part can be accessed by transpose symmetry
    vec_cutoff = default_rec_vec_cutoff
    ! NEED TO THIK ABOUT STORAGE OF RECIP VECS
    Do While( 4.0_wp * pi * Exp( - 4.0_wp * pi * vec_cutoff * vec_cutoff / ( 4.0_wp * alpha * alpha ) ) / &
         ( 4.0_wp * pi * vec_cutoff * vec_cutoff * l%V ) < 1e-12_wp )
!!$         ( 4.0_wp * pi * vec_cutoff * vec_cutoff * l%V ) < 1e-15_wp )
       vec_cutoff = vec_cutoff - 0.1_wp
    End Do
    Call generate_vec_list( l%dim, vec_cutoff, l%rec_vecs, l%rec_vec_list )

  Contains

    Pure Subroutine generate_vec_list( dim, vec_cutoff, vecs, vectors )

      Integer   ,                                 Intent( In    ) :: dim
      Real( wp ),                                 Intent( In    ) :: vec_cutoff
      Real( wp ), Dimension( :, : ),              Intent( In    ) :: vecs
      Integer   , Dimension( :, : ), Allocatable, Intent(   Out ) :: vectors

      Real( wp ), Dimension( : ), Allocatable :: vec_lengths

      Real( wp ), Dimension( 1:3 ) :: lattice_vec

      Integer, Dimension( :, : ), Allocatable :: vec_list

      Integer, Dimension( : ), Allocatable :: index

      Real( wp ) :: l1, l2, l3
      Real( wp ) :: length

      Integer :: n_vecs, max_vecs
      Integer :: g1_max, g2_max, g3_max
      Integer :: g1, g2, g3
      
      n_vecs = 0
      Select Case( dim )

      Case( 0 )
         n_vecs = 0
         max_vecs = 0
         Allocate( vec_lengths( 1:max_vecs ) )
         Allocate( vec_list( 1:0, 1:max_vecs ) )

      Case( 1 )
         l1 = Sqrt( Dot_product( vecs( :, 1 ), vecs( :, 1 ) ) )
         g1_max = Int( vec_cutoff / l1 ) + 1
         max_vecs = g1_max + 1
         Allocate( vec_lengths( 1:max_vecs ) )
         Allocate( vec_list( 1:1, 1:max_vecs ) )
         n_vecs = 0
         Do g1 = 0, g1_max
            length = g1 * l1
            If( length < vec_cutoff ) Then
               n_vecs = n_vecs + 1
               vec_list( 1, n_vecs ) = g1
               vec_lengths( n_vecs ) = length
            End If
         End Do

      Case( 2 )
         l1 = Sqrt( Dot_product( vecs( :, 1 ), vecs( :, 1 ) ) )
         l2 = Sqrt( Dot_product( vecs( :, 2 ), vecs( :, 2 ) ) )
         g1_max = Int( vec_cutoff / l1 ) + 1
         g2_max = Int( vec_cutoff / l2 ) + 1
         max_vecs = ( g1_max + 1 ) * ( 2 * g2_max + 1 )
         Allocate( vec_lengths( 1:max_vecs ) )
         Allocate( vec_list( 1:2, 1:max_vecs ) )
         n_vecs = 0
         g1 = 0
         Do g2 = 0, g2_max
            lattice_vec = g1 * vecs( :, 1 ) + g2 * vecs( :, 2 )
            length = Sqrt( Dot_product( lattice_vec, lattice_vec ) )
            If( length < vec_cutoff ) Then
               n_vecs = n_vecs + 1
               vec_list( :, n_vecs ) = [ g1, g2 ]
               vec_lengths( n_vecs ) = length
            End If
         End Do
         Do g1 = 1, g1_max
            Do g2 = - g2_max, g2_max
               lattice_vec = g1 * vecs( :, 1 ) + g2 * vecs( :, 2 )
               length = Sqrt( Dot_product( lattice_vec, lattice_vec ) )
               If( length < vec_cutoff ) Then
                  n_vecs = n_vecs + 1
                  vec_list( :, n_vecs ) = [ g1, g2 ]
                  vec_lengths( n_vecs ) = length
               End If
            End Do
         End Do

      Case( 3 )
         l1 = Sqrt( Dot_product( vecs( :, 1 ), vecs( :, 1 ) ) )
         l2 = Sqrt( Dot_product( vecs( :, 2 ), vecs( :, 2 ) ) )
         l3 = Sqrt( Dot_product( vecs( :, 3 ), vecs( :, 3 ) ) )
         g1_max = Int( vec_cutoff / l1 ) + 1
         g2_max = Int( vec_cutoff / l2 ) + 1
         g3_max = Int( vec_cutoff / l3 ) + 1
         max_vecs = ( g1_max + 1 ) * ( 2 * g2_max + 1 )* ( 2 * g3_max + 1 )
         Allocate( vec_lengths( 1:max_vecs ) )
         Allocate( vec_list( 1:3, 1:max_vecs ) )
         n_vecs = 0
         g1 = 0
         g2 = 0
         Do g3 = 0, g3_max
            lattice_vec = g1 * vecs( :, 1 ) + g2 * vecs( :, 2 ) + g3 * vecs( :, 3 )
            length = Sqrt( Dot_product( lattice_vec, lattice_vec ) )
            If( length < vec_cutoff ) Then
               n_vecs = n_vecs + 1
               vec_list( :, n_vecs ) = [ g1, g2, g3 ]
               vec_lengths( n_vecs ) = length
            End If
         End Do
         g1 = 0
         Do g2 = 1, g2_max
            Do g3 = - g3_max, g3_max
               lattice_vec = g1 * vecs( :, 1 ) + g2 * vecs( :, 2 ) + g3 * vecs( :, 3 )
               length = Sqrt( Dot_product( lattice_vec, lattice_vec ) )
               If( length < vec_cutoff ) Then
                  n_vecs = n_vecs + 1
                  vec_list( :, n_vecs ) = [ g1, g2, g3 ]
                  vec_lengths( n_vecs ) = length
               End If
            End Do
         End Do
         Do g1 = 1, g1_max
            Do g2 = - g2_max, g2_max
               Do g3 = - g3_max, g3_max
                  lattice_vec = g1 * vecs( :, 1 ) + g2 * vecs( :, 2 ) + g3 * vecs( :, 3 )
                  length = Sqrt( Dot_product( lattice_vec, lattice_vec ) )
                  If( length < vec_cutoff ) Then
                     n_vecs = n_vecs + 1
                     vec_list( :, n_vecs ) = [ g1, g2, g3 ]
                     vec_lengths( n_vecs ) = length
                  End If
               End Do
            End Do
         End Do

      Case Default
         ! Need to implement error checking - do this to avoid stupid warnings from gfortran-8
         n_vecs = 0
         max_vecs = 0
         Allocate( vec_lengths( 1:max_vecs ) )
         Allocate( vec_list( 1:0, 1:max_vecs ) )

      End Select

      If( dim /= 0 ) Then
         Allocate( index( 1:n_vecs ) )
         Call sort_index( vec_lengths( 1:n_vecs ), SORT_ASCEND, index )
         vectors = vec_list( :, [ index ] )
      Else
         Allocate( vectors( 1:0, 1:0 ) )
      End If

    End Subroutine generate_vec_list

  End Subroutine initialise

  Subroutine print( l )

    Class( lattice ), Intent( In ) :: l

    Real( wp ), Dimension( 1:3 ) :: vec 
    
    Real( wp ) :: l1, l2, l3, llat
    Real( wp ) :: alpha, beta, gamma

    Integer :: i

    l1 = Sqrt( Dot_product( l%dir_vecs( :, 1 ), l%dir_vecs( :, 1 ) ) )
    l2 = Sqrt( Dot_product( l%dir_vecs( :, 2 ), l%dir_vecs( :, 2 ) ) )
    l3 = Sqrt( Dot_product( l%dir_vecs( :, 3 ), l%dir_vecs( :, 3 ) ) )
    
    alpha = Acos( Dot_product( l%dir_vecs( :, 1 ), l%dir_vecs( :, 2 ) ) / ( l1 * l2 ) )
    beta  = Acos( Dot_product( l%dir_vecs( :, 1 ), l%dir_vecs( :, 3 ) ) / ( l1 * l3 ) )
    gamma = Acos( Dot_product( l%dir_vecs( :, 2 ), l%dir_vecs( :, 3 ) ) / ( l2 * l3 ) ) 

    alpha = alpha * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) )
    beta  = beta  * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) )
    gamma = gamma * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) )

    Write( *, '( "Dimensionality of lattice:", t30, i0   )' ) l%dim
    If( l%dim /= 0 ) then
       Write( *, '( "Volume                   :", t30, f0.6 )' ) l%V
    End If
    Write( *, '( "Given vectors            :", t30, 3(  3( f9.3, 1x ) : / t30 ) )' ) l%given_vecs
    Write( *, '( "Used vectors             :", t30, 3(  3( f9.3, 1x ) : / t30 ) )' ) l%dir_vecs
    Write( *, '( "Angles                   :", t30, 3(     f9.3, 1x )           )' ) alpha, beta, gamma
    Write( *, '( "N dir vectors          :", t30, i0 )' ) Size( l%dir_vec_list, Dim = 2 )
    If( l%dim /= 0 ) then
       Write( *, '( "First vecs               :", t30, 20( i7, 1x ) )' ) &
            l%dir_vec_list( 1, 1:Min( 20, Size( l%dir_vec_list, Dim = 2 ) ) )
    End If
    Do i = 2, l%dim
       Write( *, '( t30, 1000( i7, 1x ) )' ) &
            l%dir_vec_list( i, 1:Min( 20, Size( l%dir_vec_list, Dim = 2 ) ) )
    End Do
    If( l%dim /= 0 ) Then
       Write( *, '( "Vector lengths             :", t30 )', Advance = 'No' )
       Do i = 1, Min( 20, Size( l%dir_vec_list, Dim = 2 ) )
          vec = Matmul( l%given_vecs( :, 1:l%dim ), Real( l%dir_vec_list( 1:l%dim, i ), wp ) )
          vec = Matmul( l%dir_vecs( :, 1:l%dim ), Real( l%dir_vec_list( 1:l%dim, i ), wp ) )
          llat = Sqrt( Dot_product( vec, vec ) )
          Write( *, '( f7.2, 1x )', Advance = 'no' ) llat
       End Do
    End If
    Write( *, * )
    Write( *, '( "Recip vectors            :", t30, 3(  3( f9.3, 1x ) : / t30 ) )' ) l%rec_vecs
    Write( *, '( "Recip vectors * dir_vecs :", t30, 3(  3( f9.3, 1x ) : / t30 ) )' ) &
         Matmul( Transpose( l%rec_vecs ), l%dir_vecs )
    Write( *, '( "N recip vectors          :", t30, i0 )' ) Size( l%rec_vec_list, Dim = 2 )
    If( l%dim /= 0 ) then
       Write( *, '( "First recip vecs      :", t30, 20( i7, 1x ) )' ) &
            l%rec_vec_list( 1, 1:Min( 20, Size( l%rec_vec_list, Dim = 2 ) ) )
    End If
    Do i = 2, l%dim
       Write( *, '( t30, 20( i7, 1x ) )' ) &
            l%rec_vec_list( i, 1:Min( 20, Size( l%rec_vec_list, Dim = 2 ) ) )
    End Do
    If( l%dim /= 0 ) Then
       Write( *, '( "Vector lengths             :", t30 )', Advance = 'No' )
       Do i = 1, Min( 20, Size( l%rec_vec_list, Dim = 2 ) )
          vec = Matmul( l%rec_vecs( :, 1:l%dim ), Real( l%rec_vec_list( 1:l%dim, i ), wp ) )
          llat = Sqrt( Dot_product( vec, vec ) )
          Write( *, '( f7.2, 1x )', Advance = 'no' ) llat
       End Do
    End If
    Write( *, * )
    
  End Subroutine print

  Pure Function get_volume( l ) Result( V )

    Implicit None
    
    Real( wp ) :: V

    Class( lattice ), Intent( In ) :: l

    V = l%V
    
  End Function get_volume

  Pure Function get_direct_vectors( l ) Result( G )

    Class( lattice ), Intent( In ) :: l

    Real( wp ), Dimension( :, : ), Allocatable :: G

    G = l%dir_vecs

  End Function get_direct_vectors

  Pure Function get_reciprocal_vectors( l ) Result( G )

    Class( lattice ), Intent( In ) :: l

    Real( wp ), Dimension( :, : ), Allocatable :: G

    G = l%rec_vecs

  End Function get_reciprocal_vectors

  Pure Function get_n_dir_vecs( l ) Result( n )

    Implicit None
    
    Integer :: n

    Class( lattice ), Intent( In ) :: l

    n = Size( l%dir_vec_list, Dim = 2 )
    
  End Function get_n_dir_vecs
  
  Pure Function get_n_rec_vecs( l ) Result( n )

    Implicit None
    
    Integer :: n

    Class( lattice ), Intent( In ) :: l

    n = Size( l%rec_vec_list, Dim = 2 )
    
  End Function get_n_rec_vecs
  
  Pure Subroutine get_nth_dir_vec( l, n, vec ) 

    Implicit None

    Class( lattice )          , Intent( In    ) :: l
    Integer                   , Intent( In    ) :: n
    Real( wp ), Dimension( : ), Intent(   Out ) :: vec

    vec( 1:l%dim ) = Matmul( l%dir_vecs, Real( l%dir_vec_list( :, n ), wp ) )
    
  End Subroutine get_nth_dir_vec
  
  Pure Subroutine get_nth_rec_vec( l, n, vec ) 

    Implicit None
    
    Class( lattice )          , Intent( In    ) :: l
    Integer                   , Intent( In    ) :: n
    Real( wp ), Dimension( : ), Intent(   Out ) :: vec

    vec( 1:l%dim ) = Matmul( l%rec_vecs, Real( l%rec_vec_list( :, n ), wp ) )
    
  End Subroutine get_nth_rec_vec

  Pure Subroutine rotate_given_to_used( l, r, vec )

    Implicit None

    Class( lattice )          , Intent( In    ) :: l
    Real( wp ), Dimension( : ), Intent( In    ) :: r
    Real( wp ), Dimension( : ), Intent(   Out ) :: vec
    
    vec( 1:l%dim )  = Matmul( l%rotation, r( 1:l%dim ) )
    
  End Subroutine rotate_given_to_used

  Pure Subroutine to_reference( l, r, t ) 

    Implicit None

    Class( lattice )          , Intent( In    ) :: l
    Real( wp ), Dimension( : ), Intent( In    ) :: r
    Real( wp ), Dimension( : ), Intent(   Out ) :: t
    
    t( 1:l%dim ) = Matmul( Transpose( l%rec_vecs ), r( 1:l%dim ) )
    t( 1:l%dim ) = t( 1:l%dim ) - Floor( t( 1:l%dim ) )
    t( 1:l%dim ) = Matmul( l%dir_vecs, t( 1:l%dim ) )
    
  End Subroutine to_reference

  Subroutine minimum_image( l, r, t ) 

    Implicit None

    Class( lattice )          , Intent( In    ) :: l
    Real( wp ), Dimension( : ), Intent( In    ) :: r
    Real( wp ), Dimension( : ), Intent(   Out ) :: t

    ! Find the fractional coordinates 
    t( 1:l%dim ) = Matmul( Transpose( l%rec_vecs ), r( 1:l%dim ) )
    ! Shift into reference cell in range 0 < f < 1
    t( 1:l%dim ) = t( 1:l%dim ) - Floor( t( 1:l%dim ) )
    Where( t >= 0.5_wp )
       t = t - 1.0_wp
    End Where
    ! And back to direct
    t( 1:l%dim ) = Matmul( l%dir_vecs, t( 1:l%dim ) )
    
  End Subroutine minimum_image

  Pure Subroutine to_direct( l, f, r ) 

    Implicit None

    Class( lattice )          , Intent( In    ) :: l
    Real( wp ), Dimension( : ), Intent( In    ) :: f
    Real( wp ), Dimension( : ), Intent(   Out ) :: r

    r( 1:l%dim ) = Matmul( l%dir_vecs, f( 1:l%dim ) )
    
  End Subroutine to_direct

  Pure Subroutine to_fractional( l, r, f ) 

    Implicit None

    Class( lattice )          , Intent( In    ) :: l
    Real( wp ), Dimension( : ), Intent( In    ) :: r
    Real( wp ), Dimension( : ), Intent(   Out ) :: f

    f( 1:l%dim ) = Matmul( Transpose( l%rec_vecs ), r( 1:l%dim ) )
    
  End Subroutine to_fractional

  Pure Subroutine get_dir_vec( l, n, G )

    Implicit None
    Class( lattice )          , Intent( In    ) :: l
    Integer   , Dimension( : ), Intent( In    ) :: n
    Real( wp ), Dimension( : ), Intent(   Out ) :: G

    G( 1:l%dim ) = Matmul( l%dir_vecs, Real( n( 1:l%dim ), wp ) )
    
  End Subroutine get_dir_vec
  
  Pure Subroutine get_rec_vec( l, n, G )

    Implicit None
    Class( lattice )          , Intent( In    ) :: l
    Integer   , Dimension( : ), Intent( In    ) :: n
    Real( wp ), Dimension( : ), Intent(   Out ) :: G

    G( 1:l%dim ) = Matmul( l%rec_vecs, Real( n( 1:l%dim ), wp ) )
    
  End Subroutine get_rec_vec
  
End Module lattice_module
