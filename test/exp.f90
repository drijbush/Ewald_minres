program test_exp

  Use iso_fortran_env, Only : wp => real64
  Implicit None

  Real( wp ) :: r = 0.0_wp
  Real( wp ) :: g_2dr, g_dr, g_r, f_rdr, f_drdr
  Real( wp ), Parameter :: dr = 0.01_wp
  Real( wp ), Parameter :: alpha = 0.4_wp
  Integer,    Parameter :: nSamp = 100
  Integer :: i
  Real( wp ), Parameter :: gauss_exact(0:nSamp) = [(exp(alpha*alpha*(i*dr)**2), i = 0, nSamp)]

  g_dr = g(dr, alpha)
  g_2dr = g( Sqrt( 2.0_wp ) * dr, alpha)
  g_r = g(r, alpha)
  f_rdr = f(r, dr, alpha)
  f_drdr = f(dr, dr, alpha)

  print*, 'i    r         apprx                  exact                 diff'
  do i = 0, nSamp
    print('(i3.1,1x,f5.3,3(1x, g21.15))'), i, r, g_r, gauss_exact(i), g_r - gauss_exact(i)
    g_r = g_r * f_rdr * g_dr
    f_rdr = f_rdr * g_2dr
    r = r + dr
  end do

contains

  Real(wp) function g(r, alpha)
    Real( wp ), Intent( In    ) :: r
    Real( wp ), Intent( In    ) :: alpha

    g = exp(alpha**2 * r * r)

  end function g

  Real(wp) function f(r, dr, alpha)
    Real( wp ), Intent( In    ) :: r
    Real( wp ), Intent( In    ) :: dr
    Real( wp ), Intent( In    ) :: alpha

    f = exp(alpha**2 * r * dr)

  end function f

end program test_exp
