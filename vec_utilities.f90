Module vec3_utilities_module
    
  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64

  Implicit None

  Interface Operator( .vec. )
     Module Procedure cross
  End Interface Operator( .vec. )

  Public :: Operator( .vec. )
  Public :: vec3_invert
  
  Private

Contains
  
  Pure Function cross( a, b )

    ! perform a vector product

    Real( wp ), Dimension( 1:3 ) :: cross

    Real( wp ), Dimension( 1:3 ), Intent( In ) :: a
    Real( wp ), Dimension( 1:3 ), Intent( In ) :: b

    cross( 1 ) =   ( a( 2 ) * b( 3 ) - a( 3 ) * b( 2 ) )
    cross( 2 ) = - ( a( 1 ) * b( 3 ) - a( 3 ) * b( 1 ) )
    cross( 3 ) =   ( a( 1 ) * b( 2 ) - a( 2 ) * b( 1 ) )

  End Function cross

  Pure Subroutine vec3_invert( A, B )

    Real( wp ), Dimension( 1:3, 1:3 ), Intent( In    ) :: A
    Real( wp ), Dimension( 1:3, 1:3 ), Intent(   Out ) :: B

    Real( wp ) :: V

    V = Dot_product( A( :, 1 ), A( :, 2 ) .vec. A( :, 3 ) )
    
    B( 1, : ) = A( :, 2 ) .vec. A( :, 3 )
    B( 2, : ) = A( :, 3 ) .vec. A( :, 1 )
    B( 3, : ) = A( :, 1 ) .vec. A( :, 2 )

    B = B / V
    
  End Subroutine vec3_invert
   
End Module vec3_utilities_module
