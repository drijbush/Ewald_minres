Program gen_conf

  Real, Parameter :: a = 5.6402

  Real, Dimension( 1:3 ), Parameter :: rna_1 = [ 0.0, 0.0, 0.0 ]
  Real, Dimension( 1:3 ), Parameter :: rna_2 = [ 0.5, 0.5, 0.0 ]
  Real, Dimension( 1:3 ), Parameter :: rna_3 = [ 0.5, 0.0, 0.5 ]
  Real, Dimension( 1:3 ), Parameter :: rna_4 = [ 0.0, 0.5, 0.5 ]
  
  Real, Dimension( 1:3 ), Parameter :: rcl_1 = [ 0.5, 0.0, 0.0 ]
  Real, Dimension( 1:3 ), Parameter :: rcl_2 = [ 0.0, 0.5, 0.0 ]
  Real, Dimension( 1:3 ), Parameter :: rcl_3 = [ 0.0, 0.0, 0.5 ]
  Real, Dimension( 1:3 ), Parameter :: rcl_4 = [ 0.5, 0.5, 0.5 ]

  Real, Dimension( 1:3 ) :: r
  
  Real :: a_super
  
  Integer :: ix, iy, iz, i
  Integer :: ncell, n

  Write( *, * ) 'ncell ?'
  Read ( *, * ) ncell
  n = 8 * ncell ** 3

  a_super = a * ncell

  Open( 10, file = 'CONFIG.small' )
  Write( 10, * ) 'Small NaCl'
  Write( 10, * ) 0, 1, n
  Write( 10, * ) a_super,     0.0, 0.0
  Write( 10, * ) 0.0    , a_super, 0.0
  Write( 10, * ) 0.0    ,     0.0, a_super

  i = 0
  Do iz = 0, ncell - 1
     Do iy = 0, ncell - 1
        Do ix = 0, ncell - 1

           i = i + 1
           Write( 10, * ) 'Na+', i
           r = a * ( rna_1 + [ ix, iy, iz ] )
           Write( 10, * ) r
           i = i + 1
           Write( 10, * ) 'Na+', i
           r = a * ( rna_2 + [ ix, iy, iz ] )
           Write( 10, * ) r
           i = i + 1
           Write( 10, * ) 'Na+', i
           r = a * ( rna_3 + [ ix, iy, iz ] )
           Write( 10, * ) r
           i = i + 1
           Write( 10, * ) 'Na+', i
           r = a * ( rna_4 + [ ix, iy, iz ] )
           Write( 10, * ) r
           i = i + 1

           Write( 10, * ) 'Cl-', i
           r = a * ( rcl_1 + [ ix, iy, iz ] )
           Write( 10, * ) r
           i = i + 1
           Write( 10, * ) 'Cl-', i
           r = a * ( rcl_2 + [ ix, iy, iz ] )
           Write( 10, * ) r
           i = i + 1
           Write( 10, * ) 'Cl-', i
           r = a * ( rcl_3 + [ ix, iy, iz ] )
           Write( 10, * ) r
           i = i + 1
           Write( 10, * ) 'Cl-', i
           r = a * ( rcl_4 + [ ix, iy, iz ] )
           Write( 10, * ) r

        End Do
     End Do
  End Do

  Write( *, * ) 'n = ', n
    
End Program gen_conf
