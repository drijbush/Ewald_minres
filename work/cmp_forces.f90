Program cmp_forces
  
  Use, Intrinsic :: iso_fortran_env, Only :  wp => real64, li => int64

  Implicit None

  Real( wp ), Dimension( 1:3 ) :: fa, fb

  Real( wp ) :: f_diff_max, f_diff
  Real( wp ) :: f_diff_rel_max
  
  Integer :: iostat
  Integer :: i
  
  Character( Len = 132 ) :: file_a, file_b

  Write( *, * ) 'first file?'
  Read( * , * ) file_a
  Write( *, * ) 'second file?'
  Read( * , * ) file_b
  
  Open( 11, file = file_a )
  Open( 12, file = file_b )

  Do i = 1, 2
     Read( 11, * )
     Read( 12, * )
  End Do

  f_diff_max = - Huge( f_diff_max )
  f_diff_rel_max = - Huge( f_diff_rel_max )
  Do
     Read( 11, *, iostat = iostat ) i, fa
     If( iostat < 0 ) Exit
     Read( 12, *, iostat = iostat ) i, fb
     f_diff = Maxval( Abs( fa - fb ) )
     f_diff_max = Max( f_diff_max, f_diff )
     f_diff_rel_max = Max( f_diff_rel_max, f_diff / Maxval( Abs( fa ) ) )
  End Do

  Write( *, * ) 'Max difference ', f_diff_max
  Write( *, * ) 'Max relative difference ', f_diff_rel_max
  
End Program cmp_forces
